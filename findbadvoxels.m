% check if in sigvox
r = load('reg/retrieval_mask.mat');
roiInds = find(r.mask_brain);
fn = dir(['Localizer/locpatternsdata_' '*']);
n = load(fullfile(fn(end).name));
ind = roiInds(n.patterns.sigVox);

for i = 1:3
    folder = ['motRun' num2str(i)];
    cd(folder);
    fn = dir(['motpatternsdata_' '*']);
    n = load(fullfile(fn(end).name));
    
    test_sd = std(n.patterns.raw(:,:),[],1);
    bad = find(test_sd==0);% & n.patterns.raw(iTrial,:)<20);
    cd ../
end