%% Computation and analysis of TME metabolic interactions
%% 
% Import a file with information on the AGORA organisms including reconstruction 
% names and taxonomy.
clear all
close all 
clc
initCobraToolbox(false)
infoFile = {'ModelNames'; 'CRC_model'; 'fibro'; 'HealthyColon_model'; 'm1_model'; 'm2_model'; 'naiveTcell'; 'Th1Cell'; 'Th2Cell'; 'Th17Cell'}

%infoFile = {'ModelNames';'fibro'; 'm1_model'; 'm2_model'; 'naiveTcell'; 'Th1Cell'; 'Th2Cell'; 'Th17Cell'}

modelList=infoFile(2:end,1);
c = 400;
u = 0;
mergeGenes = false;
numWorkers = 4;
mkdir('Results');
resPath = [pwd filesep 'Results'];
modPath = pwd;

%joinModelsPairwiseFromList(modelList,modPath,'pairwiseModelFolder', resPath,'c',c,'u',u,'mergeGenesFlag',mergeGenes,'numWorkers',numWorkers);
sigD = 0.1;
[pairwiseInteractions]=simulatePairwiseInteractions(resPath,'sigD',sigD,'saveSolutionsFlag', false,'numWorkers', numWorkers);
save('pairwiseInteractions.mat', 'pairwiseInteractions')
[pairwiseInteractions_Sampled]=simulatePairwiseInteractions_Sampling(resPath,'sigD',sigD,'saveSolutionsFlag', false,'numWorkers', numWorkers);
save('pairwiseInteractions_Sampled.mat', 'pairwiseInteractions_Sampled')
% Simulate the pairwise interactions on the four dietary conditions.

%% Analysis of computed pairwise interactions
% The computed microbe-microbe interactions will be plotted by type and analyzed 
% in the context of the taxonomy of the joined strains. There are six possible 
% types of interactions total that can result in increased growth (+), no change 
% in growth (=) or decreased growth (-) compared with the single condition for 
% each joined microbe.
%% 
% * Competition (-/-)
% * Parasitism (+/-)
% * Amensalism (=/-)
% * Neutralism (=/=)
% * Commensalism (+/=)
% * Mutualism (+/+)
%% 
% This results in nine different outcomes total from the perspective of each 
% joined microbe.
% 
% Plot the percentage of interactions computed.

figure('rend','painters','pos',[10 10 900 600])
typesIA=unique(pairwiseInteractions(2:end,10));
for i = 1:length(conditions)
    pairwiseInteractions=Interactions.(conditions{i});
    listIA=pairwiseInteractions(2:end,10);
    for j=1:length(typesIA)
        dat(j)=sum(strcmp(listIA(:),typesIA{j}));
    end
    subplot(2,2,i)
    pie(dat)
    set(gca,'FontSize',10)
    h=title(conditions{i});
    set(h,'interpreter','none')
    title(conditions{i})
end
legend1=legend(typesIA);
set(legend1,'Position',[0.42 0.45 0.2 0.2],'FontSize',12)
sgtitle('Percentage of computed pairwise interactions')


%% Costless Growth
% load Models
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/CRC_model.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/fibro.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/HealthyColon_model.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/m1_model.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/m2_model.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/naiveTcell.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/Th1Cell.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/Th2Cell.mat')
load('/Users/patrickgelbach/Dropbox/MicroenvironmentModels/Th17Cell.mat')

% establish models in one struct 
models.CRC_model = CRC_model;
models.fibro = fibro;
models.HealthyColon_model = HealthyColon_model;
models.m1_model = m1_model;
models.m2_model = m2_model;
models.naiveTcell = naiveTcell;
models.Th1Cell = Th1Cell;
models.Th2Cell = Th2Cell;
models.Th17Cell = Th17Cell;

save('TMEmodels.mat', 'models')
