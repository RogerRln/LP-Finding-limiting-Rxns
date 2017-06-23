# LP-Finding-limiting-Rxns
LP Problem implemented in MATLAB that identify the Enzymatic Reactions whose bounds need to be expanded in order to achieve the experimental Growth Rate on different media conditions.

Here, I provided the data (data directory) for the Growth simulation on CH (rich medium) using SMM (poor medium) as Baseline. This example covers the issue of trying to predict growth on a rich medium using a poor medium as Baseline (training condition for our model).

Requirements for running the above simulation:
1. MATLAB
2. Cobra Toolbox library
3. solveCobraLP.m (function modified from Cobra toolbox library to solve an LP problem)
4. Finding_limitingrxns.m (I wrote this function for finding the relaxed reactions)
5. pFBA_simple.m (function modified from Cobra toolbox for finding the minimum flux through the network and returning the minimized flux and an irreversible model)
