# LP-Finding-limiting-Rxns
LP Problem implemented in MATLAB that help us to identify the Enzymatic Reactions that we must relaxe their bounds in order to achieve the experimental Growth Rate on certain growth conditions.

Here, I provided the data (inside data directory) for the simulation of growth on CH condition using SMM condition as Baseline. The example covers the issue of trying to predict growth on a rich medium using as a poor medium as a Baseline (training condition for our model).

Requirements for running the above simulation:
1. MATLAB
2. Cobra Toolbox library
3. solveCobraLP.m (function modified from Cobra toolbox library that solves the LP we gave as argument to this function)
4. OptimizedGrowth.m (I wrote this function for finding the relaxed reactions)
5. pFBA_simple.m (function modified from Cobra toolbox that finds the minimum flux through the network and returns the minimized flux and an irreversible model)
