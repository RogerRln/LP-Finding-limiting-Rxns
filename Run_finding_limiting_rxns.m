% Run the Finding_limitingrxns MATLAB function file.

% You need to have the below files on your working directory (Directory where Finding_limitingrxns.m file is placed) in order to run the function
% 1. Bacillus subtilis metabolic model (gb-2009-10-6-r69-s4.xml)
% 2. Media composition file (media_composition.xlsx)
% 3. Model data and Meta data file (Model_data_and_metadata)
% 4. File with the Enzymatic Reactions Indices (enzymatic_rxns.txt)
% 5. Python file 'gene_ratios.py', python function for calculating the Gene Expression Ratios.
% 6. Gene Expression Matrix (expression_normalized.csv, input for the python function).
% 7. pFBA_simple.m (matlab function for Running pFBA and retrieve the parsimonious flux distribution)


% Set up gurobi as LP solver (Move to the Directory where Gurobi is installed and Run the following command)
gurobi_setup

% Intialize Cobra Toolbox (Move to the Directory where the Cobra Toolbox functions are and Run the following command)
initCobraToolbox

% Run Finding_limitingrxns.m function.
% Arguments for the function: Baseline label, Tested condition label, and Optimization Number.
% (Read Finding_limitingrxns.m file for more information about Optimization Options).

% Media conditions labels and Names              
% Fru = M9 + Fructose
% Glc = M9 + Glucose
% Gly = M9 + Glycerol
% Glucon= M9+ Gluconate
% Mal = M9 + Malate
% Pyr = M9 + pyruvate
% Mal/Glc = M9 + Malate + Glucose
% LB = Luria-Bertani medium
% SMM = Spizizen's minimal salts medium
% CH = Casein hydrolysate medium
% Glut/Succ= M9 + Glutamate + Succinate 
% Optimization = {1,2}   

[solution, number_rxns, percentage_change, vector_distribution_pFBA, vector_distribution_FBA] = Finding_limitingrxns (tested_condition, baseline_condition, optimization);
