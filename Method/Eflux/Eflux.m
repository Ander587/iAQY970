% E-Flux method (Modify according to iCC541)
model=nnewmodel
% free-living state
if ~exist('fvaBounds.mat', 'file')
   modelFVA = model;
   modelFVA.lb(modelFVA.lb > 0) = 0;
   modelFVA.ub(modelFVA.ub < 0) = 0;
   [minFlux, maxFlux] = fluxVariability(modelFVA);
   [minFlux, maxFlux] = fluxVariability(modelFVA, 0, 'max', [], 2, 0 );
   % You need to ensure that the lower bound <= 0 and the upper bound >= 0.
   assert(all(minFlux <= 0) & all(maxFlux >= 0))
   save('fvaBounds.mat', 'minFlux', 'maxFlux')
else
     d = load('fvaBounds.mat');
    [minFlux, maxFlux] = deal(d.minFlux, d.maxFlux);
end
%% Prepare for mapExpressionToReactions
fbamodel = model;
% I need to change LB and UB according to Eflux requirements: 
% LB = -10000 || 0 and  UB = 10000 || 0 will be used
for i = 1:length(fbamodel.rxns)
    rxnInd = findRxnIDs(fbamodel,fbamodel.rxns{i,1});
    
    if fbamodel.lb(rxnInd) < 0 
        fbamodel.lb(rxnInd) = -10000;
    end
end

fbamodel.ub(:) = 10000;    

if (~all(fbamodel.lb == 0 | fbamodel.lb == -10000) ...
    && ~all(fbamodel.ub == 0 | fbamodel.ub == 10000))
    disp('FBA model bounds must be 0 or +/- Infinity');
    return;
end
nrxn = length(fbamodel.rxns);
fbamodel.present = true(nrxn, 1); 

expressionData_yaq.gene = fbamodel.genes;
[expressionRxns_yaq, parsedGPR_yaq] = mapExpressionToReactions(fbamodel, expressionData_yaq);

expressionData_yaqm.gene = fbamodel.genes;
[expressionRxns_yaqm, parsedGPR_yaqm] = mapExpressionToReactions(fbamodel, expressionData_yaqm);

%% E-flux using FVA bounds
% for each reaction, either both have expression values or none of them have expression values
[lbC, lbW] = deal(minFlux);
[ubC, ubW] = deal(maxFlux);

rxnWtExpression = expressionRxns_yaq >0 & expressionRxns_yaqm >0;

% for each reaction, get the condition with lower expression 
[~, ind] = min([expressionRxns_yaq, expressionRxns_yaqm], [], 2);
% for each reaction, set the bounds for the condition with lower expression
% as bound * expression ratio
% reactions with lower expression in logarithmic phase
rxnLowExpryaq = rxnWtExpression & ind == 1;
lbW(rxnLowExpryaq) = minFlux(rxnLowExpryaq) .* expressionRxns_yaq(rxnLowExpryaq) ./ expressionRxns_yaqm(rxnLowExpryaq);
ubW(rxnLowExpryaq) = maxFlux(rxnLowExpryaq) .* expressionRxns_yaq(rxnLowExpryaq) ./ expressionRxns_yaqm(rxnLowExpryaq);
% reactions with lower expression in stable phase
rxnLowExpryaqm = rxnWtExpression & ind == 2;
lbC(rxnLowExpryaqm) = minFlux(rxnLowExpryaqm) .* expressionRxns_yaqm(rxnLowExpryaqm) ./ expressionRxns_yaq(rxnLowExpryaqm);
ubC(rxnLowExpryaqm) = maxFlux(rxnLowExpryaqm) .* expressionRxns_yaqm(rxnLowExpryaqm) ./ expressionRxns_yaq(rxnLowExpryaqm);
[modelyaq, modelyaqm] = deal(model);
[modelyaq.lb, modelyaq.ub, modelyaqm.lb, modelyaqm.ub] = deal(lbW, ubW, lbC, ubC);

solution1=optimizeCbModel(modelyaq,'max','one');
solution2=optimizeCbModel(modelyaqm,'max','one');
[f1, f2] = deal(struct());
for j = 1:numel(model.rxns)
    f1.(model.rxns{j}) = solution1.x(j);
    f2.(model.rxns{j}) = solution2.x(j);
end
% symbiotic state
if ~exist('fvaBounds.mat', 'file')
   modelFVA = model;
   modelFVA.lb(modelFVA.lb > 0) = 0;
   modelFVA.ub(modelFVA.ub < 0) = 0;
   [minFlux, maxFlux] = fluxVariability(modelFVA);
   [minFlux, maxFlux] = fluxVariability(modelFVA, 0, 'max', [], 2, 0 );
   % You need to ensure that the lower bound <= 0 and the upper bound >= 0.
   assert(all(minFlux <= 0) & all(maxFlux >= 0))
   save('fvaBounds.mat', 'minFlux', 'maxFlux')
else
     d = load('fvaBounds.mat');
    [minFlux, maxFlux] = deal(d.minFlux, d.maxFlux);
end
%% Prepare for mapExpressionToReactions
fbamodel = model;
% I need to change LB and UB according to Eflux requirements: 
% LB = -10000 || 0 and  UB = 10000 || 0 will be used
for i = 1:length(fbamodel.rxns)
    rxnInd = findRxnIDs(fbamodel,fbamodel.rxns{i,1});
    
    if fbamodel.lb(rxnInd) < 0 
        fbamodel.lb(rxnInd) = -10000;
    end
end

fbamodel.ub(:) = 10000;    

if (~all(fbamodel.lb == 0 | fbamodel.lb == -10000) ...
    && ~all(fbamodel.ub == 0 | fbamodel.ub == 10000))
    disp('FBA model bounds must be 0 or +/- Infinity');
    return;
end
nrxn = length(fbamodel.rxns);
fbamodel.present = true(nrxn, 1); 

expressionData_yaqc.gene = fbamodel.genes;
[expressionRxns_yaqc, parsedGPR_yaqc] = mapExpressionToReactions(fbamodel, expressionData_yaq_c);

expressionData_yaqw.gene = fbamodel.genes;
[expressionRxns_yaqw, parsedGPR_yaqw] = mapExpressionToReactions(fbamodel, expressionData_yaq_w);

%% E-flux using FVA bounds
% for each reaction, either both have expression values or none of them have expression values
[lbC, lbW] = deal(minFlux);
[ubC, ubW] = deal(maxFlux);

rxnWtExpression = expressionRxns_yaqc >0 & expressionRxns_yaqw >0;

% for each reaction, get the condition with lower expression 
[~, ind] = min([expressionRxns_yaqc, expressionRxns_yaqw], [], 2);
% for each reaction, set the bounds for the condition with lower expression
% as bound * expression ratio
% reactions with lower expression in logarithmic phase
rxnLowExpryaqc = rxnWtExpression & ind == 1;
lbW(rxnLowExpryaqc) = minFlux(rxnLowExpryaqc) .* expressionRxns_yaqc(rxnLowExpryaqc) ./ expressionRxns_yaqw(rxnLowExpryaqc);
ubW(rxnLowExpryaqc) = maxFlux(rxnLowExpryaqc) .* expressionRxns_yaqc(rxnLowExpryaqc) ./ expressionRxns_yaqw(rxnLowExpryaqc);
% reactions with lower expression in stable phase
rxnLowExpryaqw = rxnWtExpression & ind == 2;
lbC(rxnLowExpryaqw) = minFlux(rxnLowExpryaqw) .* expressionRxns_yaqw(rxnLowExpryaqw) ./ expressionRxns_yaqc(rxnLowExpryaqw);
ubC(rxnLowExpryaqw) = maxFlux(rxnLowExpryaqw) .* expressionRxns_yaqw(rxnLowExpryaqw) ./ expressionRxns_yaqc(rxnLowExpryaqw);
[modelyaqc, modelyaqw] = deal(model);
[modelyaqc.lb, modelyaqc.ub, modelyaqw.lb, modelyaqw.ub] = deal(lbW, ubW, lbC, ubC);

solution1=optimizeCbModel(modelyaqc,'max','one');
solution2=optimizeCbModel(modelyaqw,'max','one');
[f1, f2] = deal(struct());
for j = 1:numel(model.rxns)
    f1.(model.rxns{j}) = solution1.x(j);
    f2.(model.rxns{j}) = solution2.x(j);
end