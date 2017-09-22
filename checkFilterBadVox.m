% check that the data is filerred properly
addpath('/Data1/code/motStudy05/code/')

subject = 28;
run = 2;

cd(num2str(subject));
r = load('reg/retrieval_mask.mat');
roiInds = find(r.mask_brain);
[gX gY gZ] = ind2sub(size(r.mask_brain),roiInds);


folder = ['motRun' num2str(run)];
fname = findNewestFile(folder,fullfile(folder,['motpatterns'  '*.mat']));
z = load(fname);
nTRs = size(z.patterns.raw,1);
% first check that the ones that it's setting to zero are faulty
for i = 1:length(z.patterns.allLow)
    figure;
    plot(z.patterns.raw(:,z.patterns.allLow(i)))
    hold on;
    plot(z.patterns.raw_sm_filt_z(:,z.patterns.allLow(i)), 'r')
    title(['Voxel number: ' num2str(z.patterns.allLow(i)) ' or ' '(' num2str(gX(z.patterns.allLow(i))) ',' num2str(gY(z.patterns.allLow(i))) ',' num2str(gZ(z.patterns.allLow(i))) ')'])
end

% plot brain
figure;
for i = 1:36
    subplot(6,6,i)
    imagesc(r.mask_brain(:,:,i));
end

% check low var
lv = [];
for t = 1:nTRs
    lv = [lv z.patterns.lowVar{t}];
end
lv = unique(lv)


% check that it didn't miss any
allStd = std(z.patterns.raw,[],1);
figure;
[sorted, ind] = sort(allStd);
plot(ind,sorted, '.')
title('Standard Deviation by Voxel')
xlabel('Voxel #')

vox = 84;
figure;
plot(z.patterns.raw(:,vox))
hold on;
plot(z.patterns.raw_sm_filt(:,vox), 'r')
title(['Voxel number: ' num2str(vox) ' or ' '(' num2str(gX(vox)) ',' num2str(gY(vox)) ',' num2str(gZ(vox)) ')'])
% now make sure these voxels for these trials are not in brain

cd ../