function minFlux = pFBAsimple(model)
% This function finds the minimum flux through the network and returns the
% minimized flux and an irreversible model

% convert model to irrev
%       This is buggy, doesn't deal well with constraints on irreversible
%       reactions.
modelIrrev = convertToIrreversible(model);
% add pseudo-metabolite to measure flux through network
if nargin==1,GeneOption=0;
end
if GeneOption==0, % signal that you want to minimize the sum of all gene and non-gene associated fluxes
    modelIrrev.S(end+1,:) = ones(size(modelIrrev.S(1,:)));
elseif GeneOption==1, % signal that you want to minimize the sum of only gene-associated fluxes
    %find all reactions which are gene associated
    Ind=find(sum(modelIrrev.rxnGeneMat,2)>0);
    modelIrrev.S(end+1,:) = zeros(size(modelIrrev.S(1,:)));
    modelIrrev.S(end,Ind) = 1;
end
modelIrrev.b(end+1) = 0;
modelIrrev.mets{end+1} = 'fluxMeasure';

% add a pseudo reaction that measures the flux through the network
modelIrrev = addReaction(modelIrrev,'netFlux',{'fluxMeasure'},[-1],false,0,inf,0,'','');

% set the flux measuring demand as the objective
modelIrrev.c = zeros(length(modelIrrev.rxns),1);
modelIrrev = changeObjective(modelIrrev, 'netFlux');

% minimize the flux measuring demand (netFlux)
minFlux = optimizeCbModel(modelIrrev,'min');
minFlux.xIr = minFlux.x;
minFlux.x = zeros(size(model.rxns));

% recompile into a useful vector
count = 0
for i = 1:length(modelIrrev.rxns)
    revs(i) = strcmp(modelIrrev.rxns{i}(end-1:end),'_r');
    if revs(i) == 1
        minFlux.xIr(i) = -minFlux.xIr(i);
    end
    if strcmp(modelIrrev.rxns{i}(end-1:end),'_f') == 1;
        minFlux.x(i-count) = (minFlux.xIr(i) - minFlux.xIr(i+1));
        count = count+1;
    elseif strcmp(modelIrrev.rxns{i}(end-1:end),'_b') == 1;
        p=0;
    else
        minFlux.x(i-count) = (minFlux.xIr(i));
    end
end   
end

