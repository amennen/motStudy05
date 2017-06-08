%zscore runs of data over all runs--so zscore over all timepoints
%assumes matrix is nTRs x nVoxels

function zscoredData = zscore_runs(inputData)


runMean = mean(inputData,1);
runStd = std(inputData,1);

zscoredData = (inputData - repmat(runMean,size(inputData,1),1))./repmat(runStd,size(inputData,1),1);

end