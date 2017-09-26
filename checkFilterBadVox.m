% check that the data is filerred properly
addpath('/Data1/code/motStudy05/code/')

subject = 30;
run = 3;

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
    plot(z.patterns.raw_sm_filt(:,z.patterns.allLow(i)), 'm')
    legend('raw', 'z', 'filtered')
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
for t = 1:length(z.patterns.lowVar)
    if ~isempty(z.patterns.lowVar{t})
        lv = [lv z.patterns.lowVar{t}];
    end
end
lv = unique(lv)
highLowVar = find(~ismember(lv,z.patterns.allLow))
for i = 1:length(highLowVar)
    figure;
    plot(z.patterns.raw(:,highLowVar(i)))
    hold on;
    plot(z.patterns.raw_sm_filt_z(:,highLowVar(i)), 'r')
    plot(z.patterns.raw_sm_filt(:,highLowVar(i)), 'm')
    legend('raw', 'z', 'filtered')
    title(['Voxel number: ' num2str(highLowVar(i)) ' or ' '(' num2str(gX(highLowVar(i))) ',' num2str(gY(highLowVar(i))) ',' num2str(gZ(highLowVar(i))) ')'])
end

% check that it didn't miss any
allStd = std(z.patterns.raw,[],1);
figure;
[sorted, ind] = sort(allStd);
plot(ind,sorted, '.')
title('Standard Deviation by Voxel')
xlabel('Voxel #')

vox = 1;
figure;
plot(z.patterns.raw(:,vox))
hold on;
plot(z.patterns.raw_sm_filt(:,vox), 'r')
title(['Voxel number: ' num2str(vox) ' or ' '(' num2str(gX(vox)) ',' num2str(gY(vox)) ',' num2str(gZ(vox)) ')'])
% now make sure these voxels for these trials are not in brain
tr = 1:10;
vox = 3;
r = z.patterns.raw(tr,vox);
t = z.patterns.raw_sm_filt(tr,vox);
a = z.patterns.raw_sm_filt_z(tr,vox);

m = mean(z.patterns.raw_sm_filt(tr,vox),1);
s = std(z.patterns.raw_sm_filt(tr,vox),1,1);

cd ../