function [solution, number_rxns, percentage_changed, vector_distribution_pFBA, vector_distribution_FBA] = OptimizedGrowth (gr, enz_rxns, ub, lb, optimization)
% OPtimization problem to identify 1. Minimum number of reactions to be relaxed in order to achieve a higher growth rate.
% 2. Minimum total magnitud added to the bounds of the reactions in order to achieve a higher growth rate

% Parameters:
%   gr = vector of growth rates (vector of 'mu')
%   enz_rxns = vector with the rxns_ids of the enzymatic rxns
%   ub = Upper bound (Vmax * GeX) used to constrain the original model, the original prediction of Bacillus growing on tested medium (Basline -> Medium tested)
%   lb = Lower bound used to constrain the original model (Vmin * GeX)
%   optimization = {1, any other number} depends on what optimization we want to run 1 = Minimizing number of Rxns , any other number = Minimizing the change on the bounds of Rxns. 


%% Estas son las variables que nuestro LP nos puede regresar dependiendo de la optimizacion que estemos corriendo:

% OPTIMIZACION 1 : Identify the minimum number of reactions that need to be relaxed in order to achieve a higher growth rate (Min sum(alphas), when alpha = {0,1})
%solution
%number_rxns

%OPTIMIZACION 2 : Identify the minimum magnitude of change on the bounds of certain rxns in order to achieve a higher growth rate (Min sum(alphas), when alpha = R (real number))
%solution
%number_rxns
%rxns 
%percentage_changed
%vector_distribution_pFBA
%vector_distribution_FBA


model = readCbModel('gb-2009-10-6-r69-s4', 10000);

for i = 1:1258  % We set the enzymatic reaction bounds to an arbitrarily large number
index = find(strcmp(model.rxns, enz_rxns(i)) == 1); 
if model.rev(index) ~= 0
	model.ub(index) = 1e+06;
	model.lb(index) = -1e+06;
else
	model.ub(index) = 1e+06;
	model.lb(index) = 0;
end
end

modelMedium = readCbModel('gb-2009-10-6-r69-s4', 10000);  

for i = 1:1258   
	upper = ub(i); % Vmax * GeX ratio
	lower = lb(i); % Vmin * GeX ratio
	modelMedium.ub(find(strcmp(model.rxns, enz_rxns(i)) == 1)) = upper;
	modelMedium.lb(find(strcmp(model.rxns, enz_rxns(i)) == 1)) = lower;
end


solution_lb = optimizeCbModel(modelMedium);

%LP.A
I = eye(1685);
IS = [model.S;I;I];
zero = zeros(size(model.S));

if optimization == 1
	I(find(I == 1))= 1000;
	I_mil = I;
else
	I(find(I == 1))= 1;
	I_mil = I;
end

I = eye(1685);

if optimization == 1
	I(find(I == 1))= -1000;
	I_menosmil = I;
else
	I(find(I == 1))= -1;
	I_menosmil = I;

end


mat_zero = [zero;I_mil;I_menosmil];
new_matrix =  [IS mat_zero];

%LP.c
% Primera parte agregar el metabolito artificial que tome en cuenta todas las alfas
new_matrix(end+1,:) = zeros(size(new_matrix(1,:)));
alfas = 1686:3370;
new_matrix(end,alfas) = +1;

% Segunda parte agregar la reaccion encargada de consumir el metabolito anterior
new_matrix(:,end + 1) = zeros(size(new_matrix(:,1)));
new_matrix(end, end) = -1;

%Tercera parte crear el vector "c"
c = [model.c;model.c;1];
c(427) = 0;
c(2112) = 0;
%c = ones(size(new_matrix(1,:)));
%c = c';


%LP.lb
alfas_zeros = zeros(size(model.lb));
LB = [model.lb; alfas_zeros];
% Agregar el lower bound de la ultima reaccion que agregue
LB(end+1) = 0;
%LB(end+1) = 0;


%LP.ub
alfas_ones= ones(size(model.lb));
if optimization ~= 1
	alfas_ones(1:1685,1) = 1000;
end
UP = [model.ub; alfas_ones];
%Agregar el upper bound de la ultima reaccion que agregue
UP(end+1) = 10000;
%UP(end+1) = 10000;


%LP.csense
[nMets,nRxns] = size(model.S);
csense_equality(1:nMets,1) = 'E';
csense_greater(1:nRxns,1) = 'G';
csense_less(1:nRxns,1) = 'L';
csense =[csense_equality; csense_greater; csense_less];
%Agrego el csense de la reaccion artificial
%csense(end+1) = 'E';
csense(end+1) = 'E';

%Lp.b
b_stoimatrix(1:nMets,1) = 0;
b_upper = modelMedium.ub;
b_lower = modelMedium.lb;
B =[b_stoimatrix; b_lower; b_upper];
%B(end+1) = 0;
B(end+1) = 0;

%####################################################### Primera optimizacion: Encontrar el Minimo #Rxns a cambiar in order to achieve the fix Biomass value ######################

if optimization == 1
	%LP.vtype   % Variable CONTINUA para las reaccions y BINARIA para los alfas
	[nMets,nRxns] = size(model.S);
	vtype_rxns = repmat('C', nRxns, 1);
	vtype_alphas = repmat('B', nRxns, 1);
	vtype =[vtype_rxns; vtype_alphas];
	%Agrego el vtype de la ultima reaccion (artficial)
	vtype(end+1) = 'C';

	number_rxns=[];

	for i = 1:length(gr)
    		%Se fija el growth_rate
    		LB(427) = gr(i);
    		UP(427) = gr(i);
    		% Se minimiza la funcion objetivo
    		LP.osense = +1;
    		LP.A = new_matrix;
    		LP.c = c;
    		LP.lb = LB;
    		LP.ub = UP;
    		LP.csense = csense;
    		LP.b = B;
    		LP.vtype = vtype;
    		solution = solveCobraLP(LP);
    		n_rxns = solution.obj;
    		number_rxns = [number_rxns n_rxns];
	end

	rxns = [];
	percentage_changed = [];
	vector_distribution_pFBA = [];
	vector_distribution_FBA = [];

	LP =[]; % Para preparar el LP para la proxima optimizacion


    	s = 'number_rxns_opt1_SMMcutoff.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);


%######################### Segunda optimizacion: Encontrar el Minimo de CAMBIO (magnitud) que le tengo que hacer a las reacciones a cambiar in order to achieve the fix Biomass value ######################
else
	model_original = modelMedium;
	model_original.lb(find(model_original.c)) = solution_lb.f-0.0001;
	model_original.ub(find(model_original.c)) = solution_lb.f-0.0001;
	[minFlux] = pFBA_simple(model_original);
	sum_fluxes_original = minFlux.f;
	disp(sum_fluxes_original)


	number_rxns(1:length(gr),1) = 0;
	rxns = cell(length(gr),1);
	percentage_changed=[];
	vector_distribution_pFBA = [minFlux.x(1:1685)];
	vector_distribution_FBA = [solution_lb.x];

	for i = 1:length(gr)
	    %Se fija el growth_rate
	    LB(427) = gr(i);
	    UP(427) = gr(i);
	    % Se minimiza la funcion objetivo
	    LP.osense = +1;
	    LP.A = new_matrix;
	    LP.c = c;
	    LP.lb = LB;
	    LP.ub = UP;
	    LP.csense = csense;
	    LP.b = B;
	    solution = solveCobraLP(LP);

	    x = length(find(solution.full(1686:3370) ~=0));
	    number_rxns(i) = x;

	    rxns{i} = model.rxns(find(solution.full(1686:3370) ~=0));
	    
	    ind = (find(solution.full(1686:3370) ~=0));
	    alphas = solution.full(find(solution.full(1686:3370) ~=0) + 1685);
	    model_test = modelMedium;
	    for i = 1:length(ind)
		model_test.lb(ind(i)) = model_test.lb(ind(i)) - alphas(i);
		model_test.ub(ind(i)) = model_test.ub(ind(i)) + alphas(i);
	    end
	    fba = optimizeCbModel(model_test);
	    model_test.lb(427) = fba.f-0.0001;
	    model_test.ub(427) = fba.f-0.0001;
	    sol = pFBA_simple(model_test);

	    percentage = (sol.f/sum_fluxes_original);
	    percentage_changed =[percentage_changed percentage];

	    vector_distribution_pFBA = [vector_distribution_pFBA sol.x(1:1685)];
	    vector_distribution_FBA = [vector_distribution_FBA solution.full(1:1685)];

	end

	LP =[];

    	s = 'number_rxns_opt2_SMMtoCH_protein.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);

    	s = 'sumflux_change_SMMtoCH_protein.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(percentage_changed)
        	fprintf(fileID,  '%.2f\n', percentage_changed(row));
    	end
    	fclose(fileID);

	fileID = fopen('name_rxns_SMMtoCH_protein.txt','w');
	for row = 1:length(rxns{length(gr)})
		fprintf(fileID,  '%s\n', rxns{length(gr)}{row,:});
	end
	fclose(fileID);
   
    



	%dlmwrite('fluxes_pFBA.txt', vector_distribution_pFBA);
	%dlmwrite('fluxes_FBA.txt', vector_distribution_FBA);

	%% Buscando los indices de donde ocurren las reacciones verticales (saltos en el numero de reacciones, cuellos de botella que el modelo tiene que relajar)

	%indice_superior = [];
	%indice_inferior =[];
	%for i = 1:length(gr)-1
	%	if abs(number_rxns_opt2(i+1) - number_rxns_opt2(i)) > 8
	%		indice_superior = [indice_superior i+1];
	%		indice_inferior = [indice_inferior i];
	%	end
	%end


	%% Guardando el nombre de las reacciones en pasos verticales en archivos 'txt'

	%a = indice_inferior;
	%b = indice_superior;
	%z = 0;
	%dif =[];
	%for p = 1:length(a)
	%	rxns_vertical = cell((length(rxns{b(p)}) - length(rxns{a(p)})), 1);
	%	dif = [dif (length(rxns{b(p)}) - length(rxns{a(p)}))];
	%	for i = 1:length(rxns{b(p)})
	%		if find(strcmp(rxns{b(p)}{i,:}, rxns{a(p)}) == 1)
	%			c =1;
	%		else
	%			z= z + 1;
	%			rxns_vertical{z} = (rxns{b(p)}{i,:});
	%	end
	%	end
	%	z=0;
	%	s1 = 'rxns_vertical_';
	%	s2 = int2str((length(rxns{b(p)}) - length(rxns{a(p)})));
	%	s2_1 = '_'
	%        s2_2 = int2str(p);
	%	s3 = '.txt';
	%	s = strcat(s1,s2,s2_1, s2_2,s3);
	%	fileID = fopen(s,'w');
	%	for row = 1:length(rxns_vertical)
	%		fprintf(fileID,  '%s\n', rxns_vertical{row,:});
	%	end
	%	fclose(fileID);
	%
	%end
end
end
