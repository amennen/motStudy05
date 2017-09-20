% check if in sigvox
addpath('/Data1/code/motStudy05/code/')
r = load('reg/retrieval_mask.mat');
roiInds = find(r.mask_brain);

fname = findNewestFile('Localizer',fullfile(['Localizer/locpatternsdata'  '*.mat']));
n = load(fname);
indwholething = roiInds(n.patterns.sigVox);
ind = n.patterns.sigVox;
for i = 1:3
    folder = ['motRun' num2str(i)];
    cd(folder);
    fn = dir(['motpatternsdata_' '*']);
    n = load(fullfile(fn(end).name));
    
    test_sd = std(n.patterns.raw(:,:),[],1);
    bad = find(test_sd==0);% & n.patterns.raw(iTrial,:)<20);
    i
    bad
    cd ../
end
ismember(bad,ind)