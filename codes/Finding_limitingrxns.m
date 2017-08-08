function [solution, number_rxns, percentage_change, vector_distribution_pFBA, vector_distribution_FBA] = Finding_limitingrxns (tested_condition, baseline_condition, optimization)

% We designed this LP problem with the aim of identifying the minimum number of changes on the bounds
% of the enzymatic reactions that allow a Biomasss flux prediction equal to the experimental growth rate
% of any of the 11 environmental conditions for which we have experimental data.
% We implement this LP problem on a Gene expression-bounded (using ) metabolic model. We use a vector of ascendant growth rates with the last element 
% of this vector as the experimental growth rate we want to reach. Then, we iterate over this vector and we fix the biomass reaction to the ith-element of the growth rates
% vector. Finally, we Optimize the LP problem by minimizing the number of changes on the bounds of the enzymatic reactions required to achieve the experimental growth rate.
%
% OPTIMIZATION 1:
% Minimum number of reactions to be relaxed in order to achieve a higher growth rate.
%	Minimize ∑ alphas
%	   s.t.  S*v = 0
%		 v_biomass = [growth rates vector]
%                vl - alpha*1000 <= v_enzymatic <= vu + alpha*1000
%		 alpha = {0,1}, binary variable
%
% OPTIMIZATION 2:
% Minimum Change on the bounds (quantity added to the bounds) of the reactions in order to achieve a higher growth rate.
%	Minimize ∑ alphas
%	   s.t.  S*v = 0
%		 v_biomass = [growth rates vector]
%                vl - alpha <= v_enzymatic <= vu + alpha
%		 alpha = ℝ
% Input:
%
%   tested_condition = The "label" of the Growth condition we want to test. We constrain this condition using Gene Expression and Environmental constraints.
%   baseline_condition =  The "label" of the Baseline condition, i.e. Training point of our model.
%   optimization =  1, for runnin Optimization 1. 2, for running Optimization 2.
%
% Output:
%
%  OPTIMIZATION 1 :
%  solution = LP solution returned by Gurobi. Solution contains the flux solution, including the values for the alpha variables.
%  number_rxns = Number of relaxed reactions per iteration over the growth rates vector
%
%  OPTIMIZATION 2 :
%  solution = LP solution returned by Gurobi
%  number_rxns = Number of relaxed reactions per iteration over the growth rates vector
%  percentage_change = Ratio of the Sum of Fluxes at fixed Growth rate divided by the Sum of Fluxes at initial growth rate 
%  vector_distribution_pFBA = pFBA flux distribution per iteration over the growth rates vector
%  vector_distribution_FBA = FBA flux distribution per iteration over the growth rates vector
%
% Rogelio Rodriguez 8/08/2017

%% In the first part of this script we Simulate the environmental condition of the Baseline, Run FVA on the baseline model,
%% Simulate environmental condition of the Tested model, Run FVA on the Tested model, Identify and Solve reactions  
%% directionality discrepenacies between the Baseline and the Tested model. Finally, we add gene expression to the 
%% bounds of the enzymatic reactions of the tested model and optimize the model by maximizing biomass reactions of the tested model.


% Read Bacillus subtilis metabolic model
model = readCbModel('gb-2009-10-6-r69-s4', 1e4);
[nMet, nRxns] = size(model.S);

% Read media composition file, the file contains the specific set of exchange reactions
% for each growth medium, allowing the simulation of the different conditions
[~, ~, data] = xlsread('media_composition');
conditions_ids = data(2, 1:3:32)';
media_rxns = data(5:end, 3:3:end);

% Read Observed Growth rates file, this file contains the Experimental Growth Rates
% for the different conditions.
[~, ~, data] = xlsread('Model_data_and_metadata', 'Growth rates');
growth_conditions = data(4:end, 1);
Observed_growths = data(4:end, 3);
Observed_growths =  cell2mat(Observed_growths);

baseline_index = find(strcmp(growth_conditions, baseline_condition) == 1);
observed_growth_baseline = Observed_growths(baseline_index);

tested_index = find(strcmp(growth_conditions, tested_condition) == 1);
observed_growth_tested = Observed_growths(tested_index);

% Read Exchange rates file, this file contains the specific uptake rates for 8 conditions (M9 + carbon source).
% We use the uptake rate to constrain the lower and upper bound of the Carbon-Source Exchange reaction.
[~, ~, data] = xlsread('Model_data_and_metadata', 'Exchange rates');
exchange_rates_conditions = data(3, 3:2:end);
exchange_rates_carbonNames = data(4, 3:2:end);
exchange_rates_metID = data(6:end, 1);
exchange_rates_metName = data(6:end, 2);
exchange_rates = data(6:end, 3:2:end);


%% Simulate Baseline environmental condition

% Extract the Baseline Exchange Reactions
medium_column = find(strcmp(conditions_ids, baseline_condition) == 1);

is_ex_baseline = cellfun(@(rxn) ischar(rxn), media_rxns(:, medium_column), 'UniformOutput', false);
is_ex_baseline = cell2mat(is_ex_baseline);
exchange_rxns_baseline = media_rxns(is_ex_baseline, medium_column);

% We simulate Growth on the Baseline condition as follows:
% 1. Fix the lower bound of the Carbon-Metabolite Exchange Rxns and Non-Carbon Metabolite Exchange Rxns to -12, and -50, respectively. 
% 2. Fix the Experimental Uptake rate on the lower bound the carbon-source Exchange Rxn (Only when Baseline is Minimal medium + carbon source).
% 3. Fix the lower bound of the Biomass Rxn to the Experimental Growth rate.
% 4. Run FBA 
% Following the above procedure we get biologically meaningful fluxes for CO2, O2, NH3, PO4.

exchange_rxns = model.rxns(cellfun(@(id) strcmp(id(1:3), 'EX_'), model.rxns));
is_ex_rxn = ismember(model.rxns, exchange_rxns);
model.lb(is_ex_rxn) = 0;
model.ub(is_ex_rxn) = 50;

is_ex_rxn_baseline =  ismember(model.rxns, exchange_rxns_baseline);
model.lb(is_ex_rxn_baseline) = -50;

% Fix the lower bound of the Carbon Exchange rxns to -12 (sugars, amino acids, other carbon metabolites)
is_carbon = cellfun(@(id) numel(regexp(id, '^C{1}[\d]')), model.metFormulas, 'UniformOutput', false);
is_carbon = cell2mat(is_carbon);

is_ext_carbon_met = cellfun(@(id, f) id(10) == 'e' && (f == 1), model.mets, num2cell(is_carbon));
carbon_ex_rxn_ids = cellfun(@(met) ['EX_' met(1:8) '(e)'], model.mets(is_ext_carbon_met), 'UniformOutput', false);
is_carbon_ex_rxn = ismember(model.rxns, carbon_ex_rxn_ids);
model.lb(is_carbon_ex_rxn) = max(-12, model.lb(is_carbon_ex_rxn));
model.ub(is_carbon_ex_rxn) = min( 12, model.ub(is_carbon_ex_rxn));

% Constrain to zero the outflux of Carbon-Metabolites, except for co2, acetate, pyruvate, lactate.
is_ext_carbon_met(ismember(model.mets, {'cpd00011_e[Extracellular]', 'cpd00029_e[Extracellular]', 'cpd00020_e[Extracellular]', 'cpd00159_e[Extracellular]'})) = false; %co2, acetate, pyruvate, lactate
is_carbon_ex_rxn = ismember(model.rxns, cellfun(@(met) ['EX_' met(1:8) '(e)'], model.mets(is_ext_carbon_met), 'UniformOutput', false));
%model.ub(is_carbon_ex_rxn) = 0;

% Fix the Uptake Rate of the Carbon Source (Only when Baseline is M9 + carbon-source conditions).

if find(ismember(exchange_rates_conditions, baseline_condition))
	if strcmp(baseline_condition,exchange_rates_conditions(end-1))
		model.lb(find(strcmp(model.rxns, 'EX_cpd00130(e)'))) = exchange_rates{1,7};
		model.ub(find(strcmp(model.rxns, 'EX_cpd00130(e)'))) = exchange_rates{1,7};
		model.lb(find(strcmp(model.rxns, 'EX_cpd00027(e)'))) = exchange_rates{6,7};
		model.ub(find(strcmp(model.rxns, 'EX_cpd00027(e)'))) = exchange_rates{6,7};

	elseif strcmp(baseline_condition, exchange_rates_conditions(end))
		model.lb(find(strcmp(model.rxns, 'EX_cpd00023(e)'))) = exchange_rates{8,end};
		model.ub(find(strcmp(model.rxns, 'EX_cpd00023(e)'))) = exchange_rates{8,end};
		model.lb(find(strcmp(model.rxns, 'EX_cpd00036(e)'))) = exchange_rates{end, end};
		model.ub(find(strcmp(model.rxns, 'EX_cpd00036(e)'))) = exchange_rates{end,end};
	else
		ind = find(strcmp(exchange_rates_conditions, baseline_condition));
		met_id = exchange_rates_metID(find(strcmp(exchange_rates_metName, exchange_rates_carbonNames(ind))));
		uptake = exchange_rates{find(strcmp(exchange_rates_metName, exchange_rates_carbonNames(ind))), ind};
		model.ub(find(strcmp(model.rxns, strcat('EX_',met_id,'(e)')))) = uptake;
		model.lb(find(strcmp(model.rxns, strcat('EX_',met_id,'(e)')))) = uptake;
		disp(model.lb(find(strcmp(model.rxns, strcat('EX_',met_id,'(e)')))))
 	end
end

% Allow oxygen consumption
model.ub(find(strcmp(model.rxns,'EX_cpd00007(e)') == 1)) = 0;
model.lb(find(strcmp(model.rxns, 'EX_cpd00007(e)') == 1)) = -1000;


% Run FVA on the Baseline condition, i.e. Minimize and Maximize the flux through the reactions 
% while fixing the Biomass reaction to the experimental growth rate. 

ind_biomass = find(strcmp(model.rxnNames, 'Biomass'));
model.lb(ind_biomass) = 0;
model.ub(ind_biomass) = observed_growth_baseline;

[vmin,vmax] = fluxVariability(model, 90);

% Assert Vmin < Vmax, vmin <= 0, and vmax >= 0
for i = 1:length(model.rxns)
	if vmin(i) > vmax(i)
		new_vmax = vmin(i);
		vmin(i) = vmax(i);
		vmax(i) = new_vmax;
	end

	upper = max(0, vmax(i));
	lower = min(0, vmin(i));
	vmin(i) = lower;
	vmax(i) = upper;
end


%% Simulate the Environmental condition for the Tested model 
Tested_model = readCbModel('gb-2009-10-6-r69-s4', 1e4);

% Extract the Exchange Reactions for the Tested condition
medium_column = find(strcmp(conditions_ids, tested_condition) == 1);

is_ex_tested = cellfun(@(rxn) ischar(rxn), media_rxns(:, medium_column), 'UniformOutput', false);
is_ex_tested = cell2mat(is_ex_tested);
exchange_rxns_tested = media_rxns(is_ex_tested, medium_column);

% We simulate Growth on the Tested condition as follows:
% 1. Fix the lower bound of the Carbon-Metabolite Exchange Rxns and Non-Carbon Metabolite Exchange Rxns to -30, and -100, respectively. 

% We use large values as lower bounds for the Exchange Rxns of the Tested model because we don't want to constrain
% We prefer that the Gene Expression Ratio controls has a major role on the Fluxes prediction rat.

is_ex_rxn = ismember(Tested_model.rxns, exchange_rxns);
Tested_model.lb(is_ex_rxn) = 0;
Tested_model.ub(is_ex_rxn) = 1000;

is_ex_rxn_tested =  ismember(Tested_model.rxns, exchange_rxns_tested);
Tested_model.lb(is_ex_rxn_tested) = -100;

% Fix the lower bound of the Carbon Exchange Rxns to -30 (sugars, amino acids, other carbon metabolites)
Tested_model.lb(is_carbon_ex_rxn) = max(-30, Tested_model.lb(is_carbon_ex_rxn));
Tested_model.ub(is_carbon_ex_rxn) = min( 30, Tested_model.ub(is_carbon_ex_rxn));

% Constrain to zero the outflux of Carbon-Metabolites, exceptfor co2, acetate, pyruvate, lactate.
%Tested_model.ub(is_carbon_ex_rxn) = 0;

% Allow oxygen consumption
Tested_model.ub(find(strcmp(Tested_model.rxns,'EX_cpd00007(e)') == 1)) = 0;
Tested_model.lb(find(strcmp(Tested_model.rxns, 'EX_cpd00007(e)') == 1)) = -1000;



% Run FVA on the Tested model. We retrieve the vmin and vmax bounds from the tested model
% in order to compare reaction directionality between the Baseline and the Tested model. 
% We identify and solve reaction directionality discrepancies between the the enzymatic reactions
% of the Baseline and Tested model. The above Discrepancies came from the differences on metabolite 
% composition of the Baseline and the Tested growth conditions. For example some enzymatic reactions 
% involved on transporting amino acids, sugars, etc are blocked on the Baseline and are active on the Tested condition.
% We solve this


% Fix Biomass flux to the Observed growth of the Baseline condition
Tested_model.lb(ind_biomass) = 0;
Tested_model.ub(ind_biomass) = observed_growth_tested;

% Run FVA on the Tested model
[vmin_tested,vmax_tested] = fluxVariability(Tested_model, 90);

% Assert Vmin < Vmax, vmin <= 0, vmax >= 0
for i = 1:length(Tested_model.rxns)
	if vmin_tested(i) > vmax_tested(i)
		new_vmax = vmin_tested(i);
		vmin_tested(i) = vmax_tested(i);
		vmax_tested(i) = new_vmax;
	end

	upper = max(0, vmax_tested(i));
	lower = min(0, vmin_tested(i));
	vmin_tested(i) = lower;
	vmax_tested(i) = upper;
end

Tested_model.lb(ind_biomass) = 0;
Tested_model.ub(ind_biomass) = 10000;

% Read Enzymatic reactions file
enzymatic_reactions = textread('enzymatic_rxns.txt', '%s', 'delimiter', '\n');
enzymatic_index = find(ismember(model.rxns, enzymatic_reactions));


% Solve directionality discrepancies of the enzymatic reactions.
% Most of the discrepancies came from the differences between the Baseline and Tested medium composition.

blocked_rxns_index = [];
forward_rxns_index = [];
reverse_rxns_index = [];
reversible_rxns_index = [];

for i = 1:length(enzymatic_index)
	ind = enzymatic_index(i);
	% Blocked reactions on tested condition
    if (vmax_tested(ind) == 0) & (vmin_tested(ind) == 0)
            if vmin(ind) ~= 0
                blocked_rxns_index = [blocked_rxns_index; ind];
            elseif vmax(ind) ~=  0
                blocked_rxns_index = [blocked_rxns_index; ind];
            end
        % Reactions with forward direction on tested condition    
    elseif (vmin_tested(ind) == 0) & (vmax_tested(ind) > 0)
            if  ~(vmax(ind) >  0)
                forward_rxns_index = [forward_rxns_index; ind];
            end
        % Reactions with reverse direction on tested condition
    elseif (vmin_tested(ind) < 0) & (vmax_tested(ind) == 0)
            if ~(vmin(ind) < 0)
                reverse_rxns_index = [reverse_rxns_index; ind];
            end
        % Reversible reactions on tested condition
    elseif (vmin_tested(ind)) < 0 & (vmax_tested(ind) > 0)
            if ~(vmin(ind) < 0) | ~(vmax(ind) > 0)
                reversible_rxns_index = [reversible_rxns_index; ind];
            end
    end
end

discrepancies_rxns_index = [];

if numel(blocked_rxns_index)
	discrepancies_rxns_index = [discrepancies_rxns_index;blocked_rxns_index];
end

if numel(forward_rxns_index)
	discrepancies_rxns_index = [discrepancies_rxns_index;forward_rxns_index];
end

if numel(reverse_rxns_index)
	discrepancies_rxns_index = [discrepancies_rxns_index;reverse_rxns_index];
end

if numel(reversible_rxns_index)
	discrepancies_rxns_index = [discrepancies_rxns_index;reversible_rxns_index];
end
	

% Save the index of the enzymatic reactions with no directionality problem. We later constrain 
% these reactions by Gene Expression.
enzymatic_reactions_index = enzymatic_index(~ismember(enzymatic_index, discrepancies_rxns_index));

%% Constrain the bounds of the enzymatic reactions of the Tested model

% vmin * GeX ratio <= enzymatic_reaction <= vmax * GeX ratio

% The bounds of the enzymatic reactions are constrain with vmin, vmax bounds from the Baseline condition and
% we multiply those bounds by the Gene expression ratio (calculated by gene_ratios.py script). 

% Add the Path to the directory containing gene_ratios.py script. The script reads the metabolic model, parses
% the gene reaction rules (identifying gene complexes or isoenzymes) and calculates Gene expression ratio 
% (between two conditions) for all the enzymatic reactions.

Gex_ratios = py.gene_ratios.Calculate_Gexratios(tested_condition, baseline_condition);

% Convert python list into double (MATLAB object)
Gex_ratios = double(py.array.array('d',Gex_ratios));

Gex_ratios_index = find(ismember(enzymatic_reactions, Tested_model.rxns(enzymatic_reactions_index)));

% Constrain the bounds of the enzymatic reactions

Tested_model.lb(enzymatic_reactions_index) = vmin(enzymatic_reactions_index).*(Gex_ratios(Gex_ratios_index))';
Tested_model.ub(enzymatic_reactions_index) = vmax(enzymatic_reactions_index).*(Gex_ratios(Gex_ratios_index))';

% Solve the Tested model and save the solution
solution_tested = optimizeCbModel(Tested_model);

% We Create the Growth Rates vector "gr". The LP problem iterates over this vector and minimizes the number of changes
% on the bounds of the enzymatic reactions needed to achieve the ith-growth of the 'n' size Growth Rates vector.

gr = round(solution_tested.f,2):0.01:observed_growth_tested;


%% In the second part of the script we create the structure of the LP problem, the structure contains the stoichiometric matrix, flux variables,
%% and bounds for each reaction (i.e. vmin and vmax multiplied by the Gene Expression Ratio). The structure of the LP changes whether we implement 
%% Optimization 1 or Optimization 2.


% LP problem Structure
%  A      LHS matrix 
%  b      RHS vector
%  c      Objective coeff vector
%  lb     Lower bound vector [0, Inf]
%  ub     Upper bound vector [Inf]
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).

% Optimization 1 Matrix transformation:

% |S       0|   | v |  | = 0       |
% |I   1000I| * |   |  |>= vmin*GeX|
% |I  -1000I|   | a |  |<= vmax*GeX|

% Here alpha (a) is a binary variable {0,1}

% Optimization 2 Matrix transformation:

% |S       0|   | v |  | = 0       |
% |I       I| * |   |  |>= vmin*GeX|
% |I      -I|   | a |  |<= vmax*GeX|

% Here alpha (a) is a Real number


% LP.A, Left hand dide Matrix

% Right side of the LHS Matrix
I = eye(nRxns);
IS = [model.S;I;I];

if optimization == 1
	I(find(I == 1))= 1000;
	I_thousand = I;
else
	I_positive = I;
end

I = eye(nRxns);

if optimization == 1
	I(find(I == 1))= -1000;
	I_negativethousand = I;
else
	I(find(I == 1))= -1;
	I_negative = I;
end


%Left part of the LHS Matrix
zero = zeros(size(model.S));
if optimization == 1
	mat_zero = [zero;I_thousand;I_negativethousand];
else
	mat_zero = [zero;I_positive;I_negative];
end

% Concatenate the right and left Matrices in a single Matrix (LP.A)
new_matrix =  [IS mat_zero];


% LP.c
% Our objective coeff vector, is a vector of zeros of size (nRxns + 1), except for the last entry of the vector with value 1. 
% We add an extra column to the matrix, this column works as an "artificial reaction" representing the Sum of alphas:
%          a1 a2 a3.. an 
% |                     0|
% |                     0|
% |                     .|
% |                     .|
% |0 0... +1 +1 +1...  -1| = 0 

% We Minimize the flux through the "artificial reaction" , such that [∑ alphas - "artificial reaction" flux = 0], the absolute value of the flux of the artificial reaction is
% equal to the ∑ alphas.

new_matrix(end+1,:) = zeros(size(new_matrix(1,:)));
alpha_index = nRxns+1:size(new_matrix,2);
new_matrix(end,alpha_index) = +1;


new_matrix(:,end + 1) = zeros(size(new_matrix(:,1)));
new_matrix(end, end) = -1;

% The objective function is the last column of the Matrix
c = [model.c;model.c;1];
c(427) = 0;
c(2112) = 0;


% Define LP.lb and LP.ub, i.e. vectors with the lower and upper bounds for the reactions and for the alpha variables:

% -100|-30 <= Growth_Medium Exchange reaction <= 30|1000|0
%    0 <= Non-Growth_medium Exchange reaction <= 1000|0
% -1e4 <= Reversible reaction <= 1e4
%    0 <= Irreversible reaction <= 1e4
% -1e6 <= Enzymatic reaction <= 1e6
% If optimization 1, alpha is binary, {0,1}
% If optimization 2, alpha is a real number, 0 <= alpha <= 1000
% | = OR

model = readCbModel('gb-2009-10-6-r69-s4', 1e4);

is_ex_rxn = ismember(model.rxns, exchange_rxns);
model.lb(is_ex_rxn) = Tested_model.lb(is_ex_rxn);
model.ub(is_ex_rxn) = Tested_model.ub(is_ex_rxn);

for i = 1:length(enzymatic_reactions_index) 
	index = enzymatic_reactions_index(i);
	if model.rev(index) ~= 0
		model.ub(index) = 1e6;
		model.lb(index) = -1e6;
	else
		model.ub(index) = 1e6;
		model.lb(index) = 0;
	end
end

% LP. lb
alpha_zeros = zeros(size(model.lb));
LB = [model.lb; alpha_zeros];
LB(end+1) = 0;


% LP.ub
alpha_ones= ones(size(model.ub));
if optimization ~= 1
	alpha_ones(1:nRxns,1) = 1000;
end
UP = [model.ub; alpha_ones];
UP(end+1) = 10000;


% LP.csense, Vector with Constraint Senses of the linear equations (Rows of the Matrix)
csense_equality(1:nMets,1) = 'E';
csense_greater(1:nRxns,1) = 'G';
csense_less(1:nRxns,1) = 'L';
csense =[csense_equality; csense_greater; csense_less];
csense(end+1) = 'E';

% Lp.b, the Right hand side vector containing the [vmin * GeX] vector, and [vmax * GeX] vector.
b_stoimatrix(1:nMets,1) = 0;
b_upper = Tested_model.ub;
b_lower = Tested_model.lb;
B =[b_stoimatrix; b_lower; b_upper];
B(end+1) = 0;


%%% In the third part of the script we optimize the LP problem. We minimize the ∑ alphas, while fixing the Biomass flux to the ith-element of the Growth Rates Vector.

%%%%%%%%%%%%%%%%%%%%%%%% Optimization 1 %%%%%%%%%%%%%%%%%%%%%%
if optimization == 1
	%LP.vtype, type of Variables, we define Rxns as continuous Variables and Alphas as Binary Variables
	[nMets,nRxns] = size(model.S);
	vtype_rxns = repmat('C', nRxns, 1);
	vtype_alphas = repmat('B', nRxns, 1);
	vtype =[vtype_rxns; vtype_alphas];
	vtype(end+1) = 'C';

	number_rxns=[];

	for i = 1:length(gr)
    		% We fix the growth rate
    		LB(ind_biomass) = gr(i);
    		UP(ind_biomass) = gr(i);
    		% Minimize ∑ alphas
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
		% Save the number of relaxed reactions on a vector
    		number_rxns = [number_rxns n_rxns];
	end

	LP =[];

    	s = strcat('number_rxns_opt1_',baseline_condition, 'to', tested_condition, '.txt');
    	fileID = fopen(s,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%% Optimization 2 %%%%%%%%%%%%%%%%%%%%%%
else
	Gex_constrained_model = Tested_model;
	Gex_constrained_model.lb(find(Gex_constrained_model.c)) = solution_tested.f-0.0001;
	Gex_constrained_model.ub(find(Gex_constrained_model.c)) = solution_tested.f-0.0001;
	% Run pFBA on the GeX constrained model to obtain the Sum of the total metabolic fluxes.
	[minFlux] = pFBA_simple(Gex_constrained_model);
	sum_fluxes_original = minFlux.f;

	vector_distribution_pFBA = [];
	vector_distribution_FBA = [];
	percentage_change=[];
	number_rxns(1:length(gr),1) = 0;
	rxns = cell(length(gr),1);

	vector_distribution_pFBA = [minFlux.x(1:nRxns)];
	vector_distribution_FBA = [solution_tested.x];

	for i = 1:length(gr)
	    % We fix the growth rate
	    LB(ind_biomass) = gr(i);
	    UP(ind_biomass) = gr(i);
	    % Minimize ∑ alphas
	    LP.osense = +1;
	    LP.A = new_matrix;
	    LP.c = c;
	    LP.lb = LB;
	    LP.ub = UP;
	    LP.csense = csense;
	    LP.b = B;
	    solution = solveCobraLP(LP);

	    % Identify the number of relaxed reactions, i.e. alphas with values different from zero.
	    x = length(find(solution.full(1686:3370) ~=0));
	    number_rxns(i) = x;

	    % Save the name of the relaxed reactions and Identify their indices
	    rxns{i} = model.rxns(find(solution.full(1686:3370) ~=0));
	    
	    ind = (find(solution.full(1686:3370) ~=0));
	    alphas = solution.full(find(solution.full(1686:3370) ~=0) + 1685);
	    model_test = Tested_model;

	    % Run pFBA on the relaxed model
	    model_test.lb(ind) = model_test.lb(ind) - alphas;
	    model_test.ub(ind) = model_test.ub(ind) + alphas;

	    fba = optimizeCbModel(model_test);
	    model_test.lb(ind_biomass) = fba.f-0.0001;
	    model_test.ub(ind_biomass) = fba.f-0.0001;

	    sol = pFBA_simple(model_test);
	    % Calculate the ratio of change of the Sum of Fluxes (calculated by pFBA) between the relaxed model and the original GeX constrained model
	    percentage = (sol.f/sum_fluxes_original);
	    percentage_change =[percentage_change percentage];

	    % Save the distribution of fluxes (calculated by pFBA and FBA) on matrices.
	    vector_distribution_pFBA = [vector_distribution_pFBA sol.x(1:nRxns)];
	    vector_distribution_FBA = [vector_distribution_FBA solution.full(1:nRxns)];

	end

	LP =[];

    	s1 = strcat('number_rxns_opt2_',baseline_condition, 'to', tested_condition, '.txt');
    	fileID = fopen(s1,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);

    	s2 = strcat('flux_change_opt2_',baseline_condition, 'to', tested_condition, '.txt');
    	fileID = fopen(s2,'w');
    	for row = 1:length(percentage_change)
        	fprintf(fileID,  '%.2f\n', percentage_change(row));
    	end
    	fclose(fileID);
	
	s3 = strcat('name_rxns_opt2_',baseline_condition, 'to', tested_condition, '.txt');
	fileID = fopen(s3,'w');
	for row = 1:length(rxns{length(gr)})
		fprintf(fileID,  '%s\n', rxns{length(gr)}{row,:});
	end
	fclose(fileID);
   

	%dlmwrite('fluxes_pFBA.txt', vector_distribution_pFBA);
	%dlmwrite('fluxes_FBA.txt', vector_distribution_FBA);

end
end
