r = load('reg/retrieval_mask.mat');
roiInds = find(r.mask_brain);
cd('Localizer');
fn = dir(['locpatternsdata_' '*']);
n = load(fullfile(fn(end).name));
ind = roiInds(l2.patterns.sigVox);
[gX gY gZ] = ind2sub([64 64 36],roiInds);
gZ
cd ../