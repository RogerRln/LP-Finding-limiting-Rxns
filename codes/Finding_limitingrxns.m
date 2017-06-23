function [solution, number_rxns, percentage_change, vector_distribution_pFBA, vector_distribution_FBA] = Finding_limitingrxns(gr, enz_rxns, ub, lb, optimization)
%
% We designed an LP problem to identify the reactions that we need to relax
% in order to achieve an experimental growth rate. We implemented this problem 
% in an iterative way, i.e. we used a vector of ascendant growth rates with the 
% last element of this vector as the experimental growth rate we wanted to reach.
% We iterated over this vector and so the LP calculated the minimum set of reactions 
% whose bounds are being expanded in order to reach that growth rate.
%
% OPTIMIZATION 1:
% Minimum number of reactions to be relaxed in order to achieve a higher growth rate.
%	Minimize Sum of Alphas
%	   s.t.  S*v = 0
%		 v_biomass = [growth rates vector]
%                lb - alpha*1000 <= v_enzymatic <= ub + alpha*1000
%		 alpha = {0,1}
%
% OPTIMIZATION 2:
% Minimum Change on the bounds (quantity added to the bounds) of the reactions in order to achieve a higher growth rate.
%	Minimize Sum of Alphas
%	   s.t.  S*v = 0
%		 v_biomass = [growth rates vector]
%                lb - alpha <= v_enzymatic <= ub + alpha
%		 alpha = ℝ
% Input:
%
%   gr = Vector of Growth Rates (vector of 'mu'), the LP problem iterates over this vector
%   enz_rxns = Vector with the Identifiers of the enzymatic reactions that are part of the metabolic model
%   ub = Vector of Upper bounds (Vmax * GeX ratio) used to constrain the upper bounds of the enzymatic reactions.
%   lb = Vector of Lower bounds (Vmin * GeX ratio) used to constrain the lower bounds of the enzymatic reactions.
%   optimization = {1, any other number} What Optimization we want to run 1: Minimizing number of Relaxed Rxns 
%									  any other number: Minimizing the change on the bounds of Rxns. 
%
% Output:
%
%  OPTIMIZATION 1 :
%  solution = LP solution returned by Gurobi
%  number_rxns = Number of relaxed reactions per iteration over the growth rates vector
%
%  OPTIMIZATION 2 :
%  solution = LP solution returned by Gurobi
%  number_rxns = Number of relaxed reactions per iteration over the growth rates vector
%  percentage_change = Ratio of the Sum of Fluxes at fixed Growth rate divided by the Sum of Fluxes at initial growth rate 
%  vector_distribution_pFBA = pFBA flux distribution per iteration over the growth rates vector
%  vector_distribution_FBA = FBA flux distribution per iteration over the growth rates vector

% LP problem Structure
%  A      LHS matrix
%  b      RHS vector
%  c      Objective coeff vector
%  lb     Lower bound vector
%  ub     Upper bound vector
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).


model = readCbModel('gb-2009-10-6-r69-s4', 10000);

for i = 1:1258  % We set the enzymatic reaction bounds to an arbitrarily large number (Part of the LP structure)
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
% We save the lb (Vmin * GeX ratio) and ub (Vmax * GeX ratio) vectors as the new bounds of the enzymatic reactions (RHS vector of the LP structure)
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
% First, we add an artificial metabolite that consists of putting stoichiometric coefficients of 1 in all the alpha columns.
new_matrix(end+1,:) = zeros(size(new_matrix(1,:)));
alfas = 1686:3370;
new_matrix(end,alfas) = +1;

% Second, we add an artificial reaction that "consumes" the above artificial metabolite (mimicking the consume of alphas)
new_matrix(:,end + 1) = zeros(size(new_matrix(:,1)));
new_matrix(end, end) = -1;

%Third, as we want to minimize the sum of alphas, we select the above reaction as our new objective function (we later minimize this objective function during the optimization).
c = [model.c;model.c;1];
c(427) = 0;
c(2112) = 0;


%LP.lb
alfas_zeros = zeros(size(model.lb));
LB = [model.lb; alfas_zeros];
LB(end+1) = 0;


%LP.ub
alfas_ones= ones(size(model.lb));
if optimization ~= 1
	alfas_ones(1:1685,1) = 1000;
end
UP = [model.ub; alfas_ones];
UP(end+1) = 10000;


%LP.csense
[nMets,nRxns] = size(model.S);
csense_equality(1:nMets,1) = 'E';
csense_greater(1:nRxns,1) = 'G';
csense_less(1:nRxns,1) = 'L';
csense =[csense_equality; csense_greater; csense_less];
csense(end+1) = 'E';

%Lp.b
b_stoimatrix(1:nMets,1) = 0;
b_upper = modelMedium.ub;
b_lower = modelMedium.lb;
B =[b_stoimatrix; b_lower; b_upper];
B(end+1) = 0;

%####################### Optimization 1 ######################

if optimization == 1
	%LP.vtype   % Reactions are treated as continuous Variables and Alphas as Binary Variables
	[nMets,nRxns] = size(model.S);
	vtype_rxns = repmat('C', nRxns, 1);
	vtype_alphas = repmat('B', nRxns, 1);
	vtype =[vtype_rxns; vtype_alphas];
	vtype(end+1) = 'C';

	number_rxns=[];

	for i = 1:length(gr)
    		% We fix the growth rate
    		LB(427) = gr(i);
    		UP(427) = gr(i);
    		% Minimize the objective function (Sum of Alphas)
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
	percentage_change = [];
	vector_distribution_pFBA = [];
	vector_distribution_FBA = [];

	LP =[];


    	s = 'number_rxns_opt1_SMMcutoff.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);


%####################### Optimization 2 ######################
else
	model_original = modelMedium;
	model_original.lb(find(model_original.c)) = solution_lb.f-0.0001;
	model_original.ub(find(model_original.c)) = solution_lb.f-0.0001;
	% Run pFBA on the original model (the one when optimized we get a growth rate equal to the first element of the Growth Rates Vector)
	[minFlux] = pFBA_simple(model_original);
	sum_fluxes_original = minFlux.f;


	number_rxns(1:length(gr),1) = 0;
	rxns = cell(length(gr),1);
	percentage_change=[];
	vector_distribution_pFBA = [minFlux.x(1:1685)];
	vector_distribution_FBA = [solution_lb.x];

	for i = 1:length(gr)
	    % We fix the growth rate
	    LB(427) = gr(i);
	    UP(427) = gr(i);
	    % Minimize the objective function (Sum of Alphas)
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
	    percentage_change =[percentage_change percentage];

	    vector_distribution_pFBA = [vector_distribution_pFBA sol.x(1:1685)];
	    vector_distribution_FBA = [vector_distribution_FBA solution.full(1:1685)];

	end

	LP =[];

    	s = 'number_rxns_opt2_SMMtoCH.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);

    	s = 'sumflux_change_SMMtoCH.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(percentage_change)
        	fprintf(fileID,  '%.2f\n', percentage_change(row));
    	end
    	fclose(fileID);

	fileID = fopen('name_rxns_SMMtoCH.txt','w');
	for row = 1:length(rxns{length(gr)})
		fprintf(fileID,  '%s\n', rxns{length(gr)}{row,:});
	end
	fclose(fileID);
   

	%dlmwrite('fluxes_pFBA.txt', vector_distribution_pFBA);
	%dlmwrite('fluxes_FBA.txt', vector_distribution_FBA);

end
end
