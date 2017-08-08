# LP-Finding-limiting-Rxns
Function implemented in MATLAB for identifying the needed changes on the bounds of the reactions in order to achieve the desired growth rate.

We designed this LP problem with the aim of identifying the minimum number of changes on the bounds of the enzymatic reactions that allow a Biomasss flux prediction be equal to the experimental growth rate of any of the 11 environmental conditions for which we have experimental data. We implement this LP problem on a Gene expression-bounded metabolic model. We use a vector of ascendant growth rates with the last element of this vector as the experimental growth rate we want to reach. Then, we iterate over this vector and we fix the biomass reaction to the ith-element of the growth rates vector. Finally, we Optimize the LP problem by Minimizing the number of changes on the bounds of the enzymatic reactions required to achieve the experimental growth rate.

Please install:
1. MATLAB.
2. Cobra Toolbox
3. LP solver (Gurobi, glpk, cplex) in order to run the above function.

Read:
1. Please read the Run_finding_limiting_rxns.m (go to codes directory) for more information about how to run the Finding_limitingrxns.m function. This file also contains a list of requirements for properly running the function.
