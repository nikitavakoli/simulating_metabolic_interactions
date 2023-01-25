function [pairwiseInteractions, pairwiseSolutions] = simulatePairwiseInteractions_Sampling(pairwiseModelFolder, varargin)
% This function predicts the outcome of pairwise simulations in every
% combination from a given list of pairwise models. The pairwise models
% need to be created first with the function joinModelsPairwiseFromList.
% This script requires the COBRA Toolbox function solveCobraLP. Due to the
% coupling constraints on the pairwise models, the simulations cannot
% currently be run with optimizeCbModel.
% Below is a description of all possible consequences for the two joined
% organisms. Please note that the outcomes depend on the two genome-scale
% reconstructions joined and are highly dependent on the applied
% constraints.
%
% * Competition: both organisms grow slower in co-growth than separately
%   (same outcome for both).
%
% * Parasitism: one organism grows faster in co-growth than separately, while
%   the other grows slower in co-growth than separately. Outcome for
%   faster-growing organism: Parasitism_Taker, for slower-growing organism:
%   Parasitism_Giver
%
% * Amensalism: one organism's growth is unaffected by co-growth,  while
%   the other grows slower in co-growth than separately. Outcome for
%   unaffected organism: Amensalism_Unaffected, for slower-growing organism:
%   Amensalism_Affected
%
% * Neutralism: both organisms' growths are unaffected by co-growth (same
%   outcome for both)
%
% * Commensalism: one organism's growth is unaffected by co-growth,  while
%   the other grows fatser in co-growth than separately. Outcome for
%   unaffected organism: Commensalism_Giver, for slower-growing organism:
%   Commensalism_Taker
%
% * Mutualism: both organisms growth faster in co-growth than separately
%   (same outcome for both)
%
% USAGE:
%     [pairwiseInteractions, pairwiseSolutions] = simulatePairwiseInteractions(pairwiseModelFolder, varargin)
%
% INPUTS:
%     pairwiseModelFolder:     Folder where pairwise models and pairwise model info file are located
%
% OPTIONAL INPUTS:
%     inputDiet:               Cell array of strings with three columns containing
%                              exchange reaction model abbreviations, lower bounds,
%                              and upper bounds.
%                              If no diet is input then the input pairwise models
%                              will be used with unchanged constraints.
%     saveSolutionsFlag:       If true, flux solutions are stored (may result in
%                              large output files)
%     numWorkers:              Number of workers in parallel pool if desired
%     sigD:                    Difference in growth rate that counts as significant
%                              change (by default: 10%)
%
% OUTPUTS:
%     pairwiseInteractions:    Table with computed pairwise and single growth
%                              rates for all entered microbe-microbe models
%
% OPTIONAL OUTPUT:
%     pairwiseSolutions:       Table with all computed flux solutions saved
%                              (may result in large output files)
%
% .. Authors:
%       - Almut Heinken, 02/2018
%       - Almut Heinken, 02/2020: Inputs changed for more efficient computation 

parser = inputParser();  % Parse input parameters
parser.addRequired('pairwiseModelFolder', @ischar);
parser.addParameter('inputDiet', {}, @iscell)
parser.addParameter('saveSolutionsFlag', false, @(x) isnumeric(x) || islogical(x))
parser.addParameter('numWorkers', 0, @(x) isnumeric(x))
parser.addParameter('sigD', 0.1, @(x) isnumeric(x))

parser.parse(pairwiseModelFolder, varargin{:});

pairwiseModelFolder = parser.Results.pairwiseModelFolder;
inputDiet = parser.Results.inputDiet;
saveSolutionsFlag = parser.Results.saveSolutionsFlag;
numWorkers = parser.Results.numWorkers;
sigD = parser.Results.sigD;

%% initialize COBRA Toolbox and parallel pool
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox
end
solver = CBT_LP_SOLVER;

if numWorkers > 0
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end
environment = getEnvironment();

if saveSolutionsFlag == true
    pairwiseSolutions{1, 1} = 'pairedModelID';
    pairwiseSolutions{1, 2} = 'PairwiseSolution';
    pairwiseSolutions{1, 3} = 'SingleModel1Solution';
    pairwiseSolutions{1, 4} = 'SingleModel2Solution';
end

pairwiseSolutions = {};
pairwiseInteractions{1, 1} = 'pairedModelID';
pairwiseInteractions{1, 2} = 'ModelID1';
pairwiseInteractions{1, 3} = 'ModelID2';
pairwiseInteractions{1, 4} = 'pairedGrowth_Model1_sampled';
pairwiseInteractions{1, 5} = 'pairedGrowth_Model2_sampled';
pairwiseInteractions{1, 6} = 'singleGrowth_Model1_sampled';
pairwiseInteractions{1, 7} = 'singleGrowth_Model2_sampled';
pairwiseInteractions{1, 8} = 'Outcome_Model1_sampled';
pairwiseInteractions{1, 9} = 'Outcome_Model2_sampled';
pairwiseInteractions{1, 10} = 'Total_Outcome_sampled';
pairwiseInteractions{1, 11} = 'Most Common Outcome_Model1_sampled';
pairwiseInteractions{1, 12} = 'Most Common Outcome_Model2_sampled';
pairwiseInteractions{1, 13} = 'Most Common Total_Outcome_sampled';
pairwiseInteractions{1, 14} = 'pairedGrowth_Model1_sampled';
pairwiseInteractions{1, 15} = 'pairedGrowth_Model2_sampled';
pairwiseInteractions{1, 16} = 'singleGrowth_Model1_sampled';
pairwiseInteractions{1, 17} = 'singleGrowth_Model2_sampled';
pairwiseInteractions{1, 18} = 'Median Outcome_Model1_sampled';
pairwiseInteractions{1, 19} = 'Median Common Outcome_Model2_sampled';
pairwiseInteractions{1, 20} = 'Median Common Total_Outcome_sampled';

growthRates={};
solutionsTmp={};
nSamples = 1000;
if isfile([pairwiseModelFolder filesep 'pairedModelInfo.mat'])
load([pairwiseModelFolder filesep 'pairedModelInfo.mat'],'pairedModelInfo')
else
    error('Pairwise models have not been created. Please run joinModelsPairwiseFromList first.')
end

info=pairedModelInfo;
parfor i = 1:size(info, 1)
    restoreEnvironment(environment);
    changeCobraSolver(solver, 'LP', 0, -1);
    changeCobraSolverParams('LP', 'logFile', 0);
    
    % load the model
    filename=[pairwiseModelFolder filesep info{i,1}];
    pairedModel=readCbModel(filename);
    pairedModelOrg=pairedModel;
    % if a diet was input
    if ~isempty(inputDiet)
        pairedModel = useDiet(pairedModel, inputDiet);
    end
    % for each paired model, set both biomass objective functions as
    % objectives
    biomass1 = strcat(info{i, 2}, '_', info{i, 3});
    biomass2 = strcat(info{i, 4}, '_', info{i, 5});
    model1biomass = find(ismember(pairedModel.rxns, biomass1));
    pairedModel.c(model1biomass, 1) = 1;
    model2biomass = find(ismember(pairedModel.rxns, biomass2));
    pairedModel.c(model2biomass, 1) = 1;
    % enforce minimal growth
    pairedModel = changeRxnBounds(pairedModel, biomass1, 0.001, 'l');
    pairedModel = changeRxnBounds(pairedModel, biomass2, 0.001, 'l');
    % calculate joint biomass
    solutionPaired = solveCobraLP(buildLPproblemFromModel(pairedModel,false));
    fprintf('Sampling Paired Model')
    initSampler
    P = struct;
    P.lb = pairedModel.lb;
    P.ub = pairedModel.ub;
    P.beq = pairedModel.b;
    P.Aeq = pairedModel.S;
    Model_samplingResult = sample(P, nSamples);
    Paired_samples = Model_samplingResult.samples;
    
    % separate growth
    % silence model 2 and optimize model 1
    % load the model again to avoid errors
    pairedModel=pairedModelOrg;
    % if a diet was input
    if ~isempty(inputDiet)
        pairedModel = useDiet(pairedModel, inputDiet);
    end
    pairedModel = changeObjective(pairedModel, biomass1);
    % disable flux through the second model
    pairedModel = changeRxnBounds(pairedModel, pairedModel.rxns(strmatch(strcat(info{i, 4}, '_'), pairedModel.rxns)), 0, 'b');
    % calculate single biomass
    fprintf('Sampling Model 1 alone')
    solutionSingle1 = solveCobraLP(buildLPproblemFromModel(pairedModel,false));
    P = struct;
    P.lb = pairedModel.lb;
    P.ub = pairedModel.ub;
    P.beq = pairedModel.b;
    P.Aeq = pairedModel.S;
    Model_samplingResult = sample(P, nSamples);
    Single1_samples = Model_samplingResult.samples;
    
    % silence model 1 and optimize model 2
    % load the model again to avoid errors
    pairedModel=pairedModelOrg;
    % if a diet was input
    if ~isempty(inputDiet)
        pairedModel = useDiet(pairedModel, inputDiet);
    end
    pairedModel = changeObjective(pairedModel, biomass2);
    % disable flux through the first model
    pairedModel = changeRxnBounds(pairedModel, pairedModel.rxns(strmatch(strcat(info{i, 2}, '_'), pairedModel.rxns)), 0, 'b');
    % calculate single biomass
    solutionSingle2 = solveCobraLP(buildLPproblemFromModel(pairedModel,false));
    fprintf('Sampling Model 2 alone')
    P = struct;
    P.lb = pairedModel.lb;
    P.ub = pairedModel.ub;
    P.beq = pairedModel.b;
    P.Aeq = pairedModel.S;
    Model_samplingResult = sample(P, nSamples);
    Single2_samples = Model_samplingResult.samples;
    % save results temporarily
        %growthRates{i}{1,1} = solutionPaired.full(model1biomass);
        growthRates{i}{1,1} = Paired_samples(model1biomass,:);
        %growthRates{i}{1,2} = solutionPaired.full(model2biomass);
        growthRates{i}{1,2} = Paired_samples(model2biomass,:);

       % growthRates{i}{1,3} = solutionSingle1.full(model1biomass);
        growthRates{i}{1,3} = Single1_samples(model1biomass,:);

        %growthRates{i}{1,4} = solutionSingle2.full(model2biomass);
        growthRates{i}{1,4} =Single2_samples(model2biomass,:);
        
 
    solutionsTmp{i}{1,1} = solutionPaired;
    solutionsTmp{i}{1,2} = solutionSingle1;
    solutionsTmp{i}{1,3} = solutionSingle2;

    
    
end

% backup results for large number of simulations
if size(pairedModelInfo, 1) > 1000
save('growthRates','growthRates','-v7.3');
save('solutionsTmp','solutionsTmp','-v7.3');
end

for i = 1:size(pairedModelInfo, 1)
    solPaired = solutionsTmp{i}{1,1};
    solSingle1 = solutionsTmp{i}{1,2};
    solSingle2 = solutionsTmp{i}{1,3};
    
    grPaired1 = growthRates{i}{1,1};
    grPaired2 = growthRates{i}{1,2};
    grSingle1 = sort(growthRates{i}{1,3});
    grSingle2 = sort(growthRates{i}{1,4});
    
    mPaired1 = median(growthRates{i}{1,1});
    mPaired2 = median(growthRates{i}{1,2});
    mSingle1 = median(growthRates{i}{1,3});
    mSingle2 = median(growthRates{i}{1,4});
    if solPaired.stat ~= 0 && solSingle1.stat ~= 0 && solSingle2.stat ~= 0
        [iAFirstModel_Med, iASecondModel_Med, iATotal_Med] = analyzePairwiseInteractions(mPaired1, mPaired2, mSingle1, mSingle2, sigD);
        
    end
    
    
    First = {};
    Second = {};
    Total = {};
    mins = min([size(grSingle1,2), size(grSingle2,2), size(grPaired1,2)]);
    %for z = 1:size(grSingle1,2)
    for z = 1:mins  
        grSingle1z = grSingle1(z);
        %for y = 1:size(grSingle2,2)
        for y = 1:mins
            grSingle2y = grSingle2(y);
            FirstModel_all = {};
            SecondModel_all = {};
            TotalAll = {};
            %for x = 1:size(grPaired1,2)
            for x = 1:mins
                grPaired1x = grPaired1(x);
                grPaired2x = grPaired2(x);
                
                if solPaired.stat ~= 0 && solSingle1.stat ~= 0 && solSingle2.stat ~= 0
                    [iAFirstModel, iASecondModel, iATotal] = analyzePairwiseInteractions(grPaired1x, grPaired2x, grSingle1z, grSingle2y, sigD);
                    FirstModel_all(x) = cellstr(iAFirstModel);
                    SecondModel_all(x) = cellstr(iASecondModel);
                    TotalAll(x) = cellstr(iATotal);
                    
                end
            end
            [s,~,j]=unique(FirstModel_all);
            First(z,y) = cellstr(s{mode(j)});
            [s,~,j]=unique(SecondModel_all);
            Second(z,y) = cellstr(s{mode(j)});
            [s,~,j]=unique(TotalAll);
            Total(z,y) = cellstr(s{mode(j)});
          
        end
    end
    
    [s,~,j]=unique(First);
    mostCommon1 =cellstr(s{mode(j)});
    [s,~,j]=unique(Second);
    mostCommon2 = cellstr(s{mode(j)});
    [s,~,j]=unique(Total);
    mostCommonPaired = cellstr(s{mode(j)});

%     % interpret computed growth rates-only if solutions are feasible
%     if solPaired.stat ~= 0 && solSingle1.stat ~= 0 && solSingle2.stat ~= 0
%         [iAFirstModel, iASecondModel, iATotal] = analyzePairwiseInteractions(grPaired1, grPaired2, grSingle1, grSingle2, sigD);
%     end
%     
    %% list all results
    % fill out the results table with the model information
    pairwiseInteractions{i + 1, 1} = pairedModelInfo{i, 1};
    pairwiseInteractions{i + 1, 2} = pairedModelInfo{i, 2};
    pairwiseInteractions{i + 1, 3} = pairedModelInfo{i, 4};
    % combined growth
    pairwiseInteractions{i + 1, 4} = grPaired1;
    pairwiseInteractions{i + 1, 5} = grPaired2;
    pairwiseInteractions{i + 1, 6} = grSingle1;
    pairwiseInteractions{i + 1, 7} = grSingle2;
    % make sure that infeasible solutions are flagged
    if solPaired.stat == 0 || solSingle1.stat == 0 || solSingle2.stat == 0
        pairwiseInteractions{i + 1, 8} = 'Infeasible';
        pairwiseInteractions{i + 1, 9} = 'Infeasible';
        pairwiseInteractions{i + 1, 10} = 'Infeasible';
        pairwiseInteractions{i + 1, 11} = 'Infeasible';
        pairwiseInteractions{i + 1, 12} = 'Infeasible';
        pairwiseInteractions{i + 1, 13} = 'Infeasible';
    else
        % save the results
        pairwiseInteractions{i + 1, 8} = First;
        pairwiseInteractions{i + 1, 9} = Second;
        pairwiseInteractions{i + 1, 10} = Total;
        pairwiseInteractions{i + 1, 11} = mostCommon1;
        pairwiseInteractions{i + 1, 12} = mostCommon2;
        pairwiseInteractions{i + 1, 13} = mostCommonPaired;
        pairwiseInteractions{i + 1, 14} = mPaired1;
        pairwiseInteractions{i + 1, 15} = mPaired2;
        pairwiseInteractions{i + 1, 16} = mSingle1;
        pairwiseInteractions{i + 1, 17} = mSingle2;             
        pairwiseInteractions{i + 1, 18} = iAFirstModel_Med;
        pairwiseInteractions{i + 1, 19} = iASecondModel_Med;
        pairwiseInteractions{i + 1, 20} = iATotal_Med;        
    end
    %% store the solutions if storeSolutionFlag==true
    if nargin > 3 && saveSolutionsFlag == true
        pairwiseSolutions{i + 1, 1} = pairedModelInfo{i, 1};
        pairwiseSolutions{i + 1, 2} = solPaired.full;
        pairwiseSolutions{i + 1, 3} = solSingle1.full;
        pairwiseSolutions{i + 1, 4} = solSingle2.full;
    end
end

end

function [iAFirstOrg, iASecondOrg, iATotal] = analyzePairwiseInteractions(pairedGrowthOrg1, pairedGrowthOrg2, singleGrowthOrg1, singleGrowthOrg2, sigD)
% This function evaluates the outcome of pairwise growth of two organisms
% compared with the two organisms separately. There are six possible
% outcomes of the interaction, and nine possible outcomes from the
% perspective of each organism.

if abs(1 - (pairedGrowthOrg1 / singleGrowthOrg1)) < sigD
    % first microbe unaffected - all possible cases resulting froms
    % second microbe's growth
    if abs(1 - (pairedGrowthOrg2 / singleGrowthOrg2)) < sigD
        % second microbe unaffected
        iAFirstOrg = 'Neutralism';
        iASecondOrg = 'Neutralism';
        iATotal = 'Neutralism';
    elseif abs((pairedGrowthOrg2 / singleGrowthOrg2)) > 1 + sigD
        % second microbe grows better
        iAFirstOrg = 'Commensalism_Giver';
        iASecondOrg = 'Commensalism_Taker';
        iATotal = 'Commensalism';
    elseif abs((singleGrowthOrg2 / pairedGrowthOrg2)) > 1 + sigD
        % second microbe grows slower
        iAFirstOrg = 'Amensalism_Unaffected';
        iASecondOrg = 'Amensalism_Affected';
        iATotal = 'Amensalism';
    else
        % if no case fits - needs inspection!
        iAFirstOrg = 'No_Result';
        iASecondOrg = 'No_Result';
        iATotal = 'No_Result';
    end
elseif abs((pairedGrowthOrg1 / singleGrowthOrg1)) > 1 + sigD
    % first microbe grows better - all possible cases resulting froms
    % second microbe's growth
    if abs(1 - (pairedGrowthOrg2 / singleGrowthOrg2)) < sigD
        % second microbe unaffected
        iAFirstOrg = 'Commensalism_Taker';
        iASecondOrg = 'Commensalism_Giver';
        iATotal = 'Commensalism';
    elseif abs((pairedGrowthOrg2 / singleGrowthOrg2)) > 1 + sigD
        % second microbe grows better
        iAFirstOrg = 'Mutualism';
        iASecondOrg = 'Mutualism';
        iATotal = 'Mutualism';
    elseif abs((singleGrowthOrg2 / pairedGrowthOrg2)) > 1 + sigD
        % second microbe grows slower
        iAFirstOrg = 'Parasitism_Taker';
        iASecondOrg = 'Parasitism_Giver';
        iATotal = 'Parasitism';
    else
        % if no case fits - needs inspection!
        iAFirstOrg = 'No_Result';
        iASecondOrg = 'No_Result';
        iATotal = 'No_Result';
    end
elseif abs((singleGrowthOrg1 / pairedGrowthOrg1)) > 1 + sigD
    % first microbe grows slower - all possible cases resulting froms
    % second microbe's growth
    if abs(1 - (pairedGrowthOrg2 / singleGrowthOrg2)) < sigD
        % second microbe unaffected
        iAFirstOrg = 'Amensalism_Affected';
        iASecondOrg = 'Amensalism_Unaffected';
        iATotal = 'Amensalism';
    elseif abs((pairedGrowthOrg2 / singleGrowthOrg2)) > 1 + sigD
        % second microbe grows better
        iAFirstOrg = 'Parasitism_Giver';
        iASecondOrg = 'Parasitism_Taker';
        iATotal = 'Parasitism';
    elseif abs((singleGrowthOrg2 / pairedGrowthOrg2)) > 1 + sigD
        % second microbe grows slower
        iAFirstOrg = 'Competition';
        iASecondOrg = 'Competition';
        iATotal = 'Competition';
    else
        % if no case fits - needs inspection!
        iAFirstOrg = 'No_Result';
        iASecondOrg = 'No_Result';
        iATotal = 'No_Result';
    end
else
    % if no case fits - needs inspection!
    iAFirstOrg = 'No_Result';
    iASecondOrg = 'No_Result';
    iATotal = 'No_Result';
end
end
