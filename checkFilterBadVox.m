% check that the data is filerred properly
z = load('../data/28/motRun1/motpatternsdata_21_20170920T170707.mat');
figure;
plot(z.patterns.raw_sm_filt_z(:,z.patterns.allLow))

for i = 1:length(z.patterns.allLow)
    figure;
    plot(z.patterns.raw(:,z.patterns.allLow(i)))
end

plot(z.patterns.raw_sm_filt_z(:,z.patterns.allLow))

% now make sure these voxels for these trials are not in brain