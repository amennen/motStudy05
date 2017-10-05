% so you have the run, bad voxel number
folder = ['../data/' num2str(subject) '/Localizer']
fname = findNewestFile(folder,fullfile(folder,['loctrainedModel'  '*.mat']));
p = load(fname);
weights = p.trainedModel.ridge.betas;
fname = findNewestFile(folder,fullfile(folder,['locpatterns'  '*.mat']));
p = load(fname);
sigVox = p.patterns.sigVox;
for i = 1:length(z.patterns.allLow)
    
    
end