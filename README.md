# LP-Finding-limiting-Rxns
LP Problem implemented in MATLAB that help us to identify the Enzymatic Reactions that we must relaxe their bounds in order to achieve the experimental Growth Rate on certain growth conditions.

Here, I provided the data (see data directory) for the simulation of growth on CH condition using SMM condition as Baseline.

Requirements for running the above simulation:
1. MATLAB
2. Cobra-Toolbox library
3. solveCobraLP.m (code that I modified from Cobra-toolbox library in order to allow binary variables in the LP formulation)
4. OptimizedGrowth.m (function for finding the relaxed reactions)
