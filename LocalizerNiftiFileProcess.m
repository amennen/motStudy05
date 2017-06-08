function LocalizerNiftiFileProcess(subjectNum,crossval,featureSelect,prev,scanNow,scanNum,SESSION,runNum)
% now going to be loading in localizer data
if IsLinux
    biac_dir = '/Data1/packages/BIAC_Matlab_R2014a/';
    bxhpath='/opt/BXH/1.11.1/bin/';
    fslpath='/opt/fsl/5.0.9/bin/';
    dcm2path = '/opt/MRICROGL/2-2016/';
end

%add necessary package
if ~exist('readmr','file')
    addpath(genpath(biac_dir));
    addpath([biac_dir '/mr/']);
    addpath([biac_dir '/general/'])
end
projectName = 'motStudy04';
multipath = '/Data1/code/multibandutils/';
addpath(genpath(multipath));
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; 
process_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/' 'reg' '/'];
patterns_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/' 'Localizer' '/'];
if ~exist(patterns_dir, 'dir')
    mkdir(patterns_dir)
end
cd(patterns_dir)
behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
roi_dir = ['/Data1/code/' projectName '/data/'];
code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
addpath(genpath(code_dir));

if ~prev %if getting data today
    subjectName = [datestr(now,5) datestr(now,7) datestr(now,11) num2str(runNum) '_' projectName];
    dicom_dir = ['/Data1/subjects/' datestr(now,10) datestr(now,5) datestr(now,7) '.' subjectName '.' subjectName '/'];
else
    %for subject 16 doing again
    allDates = {'4-5-17'};
    %allDates = {'7-1-2016' '3-26-2016', '3-29-2016', '4-1-2016', '4-27-2016', '4-29-2016', '5-05-2016'};
    subjectName = [datestr(allDates{1},5) datestr(allDates{1},7) datestr(allDates{1},11) num2str(runNum) '_' projectName];
    dicom_dir = ['/Data1/subjects/' datestr(allDates{1},10) datestr(allDates{1},5) datestr(allDates{1},7) '.' subjectName '.' subjectName '/'];
end

roi_name = 'retrieval';
while ~exist(fullfile(process_dir,[roi_name '_mask.mat']), 'file')%wait for mask to appear
end
temp = load(fullfile(process_dir,[roi_name '_mask.mat'])); %should be stretched_brain
roi = logical(temp.mask_brain);
assert(exist('roi','var')==1);
roiDims = size(roi);
roiInds = find(roi);

%% Boilerplate

seed = sum(100*clock); %get random seed
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));%set seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));%set seed

%initialize system time calls
GetSecs;
%% preprocessing parameters
FWHM = 5;
cutoff = 160;
TR = 1;

%% block sequence
nTRsTotal = 1376;
remove = 20/TR; % because it's 20 seconds we want to take off 
nTRs = nTRsTotal - remove; 
patterns.fileAvail = zeros(1,nTRs);
patterns.newFile = cell(1,nTRs);
patterns.timeRead = cell(1,nTRs);
patterns.fileload = NaN(1,nTRs);
patterns.raw = NaN(nTRs,numel(roiInds));
patterns.raw_sm = NaN(nTRs,numel(roiInds)); %here we use the corrected number of TRs**
patterns.raw_sm_filt = NaN(nTRs,numel(roiInds)); 
patterns.raw_sm_z = NaN(nTRs,numel(roiInds)); 

patterns.block = ones(1,nTRsTotal);
%% Output File Setup
dataFile = fullfile(patterns_dir, 'fileprocessing.txt');
printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'* Mot Real Time v.1.0\n');
printlog(dataFile,['* Date/Time: ' datestr(now,0) '\n']);
printlog(dataFile,['* Seed: ' num2str(seed) '\n']);
printlog(dataFile,['* Subject Number: ' num2str(subjectNum) '\n']);
printlog(dataFile,['* Subject Name: ' subjectName '\n']);
printlog(dataFile,['* Run Number: ' num2str(runNum) '\n']);
printlog(dataFile,['* Real-Time Data: ' num2str(scanNow) '\n']);
printlog(dataFile,'*********************************************\n\n');


%% Start Experiment

% prepare for trial sequence
printlog(dataFile,'run\tblock\tTR\tloaded\n');


%% acquire files

fileCounter = 0;
for iTrial = 1:nTRs % the first 10 TRs have been taken out to detrend
    % increase count of TRs
    %increase the count of TR pulses
    fileCounter = fileCounter+1;
    thisTR = iTrial + remove; %account for taking out TRs
    while ~patterns.fileAvail(iTrial)
        [patterns.fileAvail(iTrial) patterns.newFile{iTrial}] = GetSpecificFMRIFile(dicom_dir,scanNum,thisTR);
    end
    
    %if desired file is recognized, pause for 200ms to complete transfer
    if scanNow==1 && patterns.fileAvail(iTrial)
       pause(.2);
    end
    
    % if file available, load it
    if (patterns.fileAvail(iTrial))
        
        t0 = GetSecs;
        niftiname = sprintf('nifti%3.3i', thisTR);
        unix(sprintf('%sdcm2niix %s -f %s -o %s -s y %s%s',dcm2path,dicom_dir,niftiname,patterns_dir,dicom_dir,patterns.newFile{iTrial}))
        t1 = GetSecs;
        unix(sprintf('%smcflirt -in %s.nii -reffile %sexfunc_re.nii',fslpath,niftiname,process_dir))
        t2 = GetSecs;
        moco = t2-t1;
%         
        niftiname = sprintf('nifti%3.3i_mcf.nii.gz', thisTR); 
        niftidata = readnifti(niftiname);
        newVol = niftidata(roi);
        patterns.raw(iTrial,:) = newVol;  % keep patterns for later training
        
        if (any(isnan(patterns.raw(iTrial,:)))) && (iTrial>1)
            patterns.fileload(iTrial) = 0; %mark that load failed
            patterns.raw(iTrial,:) = patterns.raw(iTrial-1,:); %replicate last complete pattern
        else
            patterns.fileload(iTrial) = 1;
        end
        
        patterns.raw(iTrial,:) = newVol;  % keep patterns for later training
        
        if (any(isnan(patterns.raw(iTrial,:)))) && (iTrial>1)
            patterns.fileload(iTrial) = 0; %mark that load failed
            patterns.raw(iTrial,:) = patterns.raw(iTrial-1,:); %replicate last complete pattern
        else
            patterns.fileload(iTrial) = 1;
        end
        
    end
    %smooth files
    patterns.raw_sm(iTrial,:) = SmoothRealTime(patterns.raw(iTrial,:),roiDims,roiInds,FWHM);
    
    % print trial results
    printlog(dataFile,'%d\t%d\t%d\t%d\n',runNum,patterns.block(iTrial),thisTR,patterns.fileAvail(iTrial));
    %fprintf('%d\t%d\t%d\t%d\n',runNum,patterns.block(iTrial),thisTR,patterns.fileAvail(iTrial));
end

%% pre-process

printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'beginning model preprocessing...\n');

preprocStart = tic; %start timing 

% detrend
patterns.raw_sm_filt = HighPassBetweenRuns(patterns.raw_sm,TR,cutoff);

% z-score
patterns.runMean = mean(patterns.raw_sm_filt,1);  %mean across all volumes per voxel
patterns.runStd = std(patterns.raw_sm_filt,[],1); %std dev across all volumes per voxel
patterns.raw_sm_filt_z = (patterns.raw_sm_filt - repmat(patterns.runMean,size(patterns.raw_sm_filt,1),1))./repmat(patterns.runStd,size(patterns.raw_sm_filt,1),1);
%localizerStd = patterns.runStd;
preprocTime = toc(preprocStart);  %end timing

% print training timing and results

save([patterns_dir 'locpreprocpatternsdata_' num2str(runNum) '_' datestr(now,30)],'patterns');
%save([save_dir 'localizerStd'], 'localizerStd');
printlog(dataFile,'data preprocessing time: \t%.3f\n',preprocTime);

%% training: cross-validation
%first cross-validate
if crossval
%print xval results
printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'beginning model cross-validation...\n');

%parameters
penalty = 100;
keepTR = 30; %%changed it on 1/20 for subjects 6 on!
shiftTR = 4/TR;
startXVAL = tic;

%first get session information
[newpattern t] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir,keepTR);
patterns.regressor.allCond = newpattern.regressor.allCond;
patterns.regressor.twoCond = newpattern.regressor.twoCond;
patterns.selector.xval = newpattern.selector.xval;
patterns.selector.allxval = newpattern.selector.allxval;
nIter = size(patterns.selector.allxval,1);
%shift regressor
nCond = size(patterns.regressor.twoCond,1);
for j = 1:nIter
   selector = patterns.selector.allxval(j,:); 
   easyIdx = find(patterns.regressor.allCond(2,:));
   hardIdx = find(patterns.regressor.allCond(1,:));
   trainIdx = find(selector == 1);
   %trainIdx = intersect(hardIdx,trainIdx);
   testIdx = find(selector == 2);
    
   % now shift indices forward
   %trainIdx = trainIdx + shiftTR;
   %testIdx = testIdx + shiftTR;
   
   trainPats = patterns.raw_sm_filt_z(trainIdx+shiftTR,:);
   testPats = patterns.raw_sm_filt_z(testIdx+shiftTR,:);
   trainTargs = patterns.regressor.twoCond(:,trainIdx);
   testTargs = patterns.regressor.twoCond(:,testIdx);
   
   if featureSelect
       thr = 0.1;
       p = run_mathworks_anova(trainPats',trainTargs);
       sigVox = find(p<thr);
       trainPats = trainPats(:,sigVox);
       testPats = testPats(:,sigVox);
   end
   
   scratchpad = train_ridge(trainPats,trainTargs,penalty);
   [acts scratchpad] = test_ridge(testPats,testTargs,scratchpad);
   %acts is nCond x nVoxels in the mask
   
   
   %calculate AUC for JUST TARGET vs. LURE
   for i = 1:length(acts)
      condition = find(testTargs(:,i));
      if condition == 1
          labels{i} = 'target';
      elseif condition == 2
          labels{i} = 'lure';
      end
   end
   [X,Y,t,AUC(j)] = perfcurve(labels,acts(1,:), 'target');
   
   %calculate AUC SEPARATELY for easy targets vs. lure && hard targets vs.
   %lure
   testTargsFour = patterns.regressor.allCond(:,testIdx);
   hardIdx = find(testTargsFour(1,:)==1);
   easyIdx = find(testTargsFour(2,:)==1);
   lureIdx = find(testTargs(2,:)==1);
   
   actsHard = acts(1,[hardIdx lureIdx]);
   actsEasy = acts(1,[easyIdx lureIdx]);
   for i = 1:length(actsHard)
       if i <= length(hardIdx)
           labelsHard{i} = 'target';
       else
           labelsHard{i} = 'lure';
       end
   end
   [X,Y,t,AUC_hard(j)] = perfcurve(labelsHard,actsHard, 'target');
   [X,Y,t,AUC_easy(j)] = perfcurve(labelsHard,actsEasy, 'target');
   fprintf(['* Completed Iteration ' num2str(j) '; AUC = ' num2str(AUC(j)) '\n']); 
   fprintf(['* Hard vs. Lure AUC = ' num2str(AUC_hard(j)) '\n']);
   fprintf(['* Easy vs. Lure AUC = ' num2str(AUC_easy(j)) '\n']);    
end

average_AUC = mean(AUC);
std_AUC = std(AUC)/sqrt(nIter-1);
average_hardAUC = mean(AUC_hard);
average_easyAUC = mean(AUC_easy);
std_hardAUC = std(AUC_hard)/sqrt(nIter-1);
std_easyAUC = std(AUC_easy)/sqrt(nIter-1);
xvaltime = toc(startXVAL); %end timing
%print cross-validation results
printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'finished cross-validation...\n');
printlog(dataFile,['* Average AUC over Iterations: ' num2str(average_AUC) ' +- ' num2str(std_AUC) '\n']);
printlog(dataFile,['* Average Hard vs. Lure AUC over Iterations: ' num2str(average_hardAUC) ' +- ' num2str(std_hardAUC) '\n']);
printlog(dataFile,['* Average Easy vs. Lure AUC over Iterations: ' num2str(average_easyAUC) ' +- ' num2str(std_easyAUC) '\n']);
printlog(dataFile, 'Cross-validation model training time: \t%.3f\n',xvaltime);

end
%% training on all data, no cross-validation
%print training results
printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'beginning model training...\n');

%parameters
penalty = 100;
shiftTR = 4/TR;
keepTR = 30; %changed on 1/20 for subject 6 onwards to train the localizer task on all 15 TR's instead of just the first one
trainStart = tic;

%first get session information
[newpattern t] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir,keepTR);
patterns.regressor.twoCond = newpattern.regressor.twoCond;
patterns.regressor.allCond = newpattern.regressor.allCond;

trainIdx = find(newpattern.selector.xval); %find all nonzero timepoints to train on
%trainLabels = patterns.regressor.twoCond([1 2],trainIdx);
trainPats = patterns.raw_sm_filt_z(trainIdx+shiftTR,:);
trainTargs = patterns.regressor.twoCond(:,trainIdx);

% try feature selection
if featureSelect
    thr = 0.1;
    p = run_mathworks_anova(trainPats',trainTargs);
    sigVox = find(p<thr);
    trainPats = trainPats(:,sigVox);
    patterns.sigVox = sigVox;
end
trainedModel = train_ridge(trainPats,trainTargs,penalty); %the weights are trainedModel.ridge.betas

trainingOnlyTime = toc(trainStart);  %end timing

if ~crossval %only save if we're not cross-validating (assuming cross-val is just for testing)
save([patterns_dir 'locpatternsdata_' num2str(runNum) '_' datestr(now,30)],'patterns');
save([patterns_dir 'loctrainedModel_' num2str(runNum) '_' datestr(now,30)],'trainedModel','trainPats','t');
end

%print training timing and results
printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'finished model training...\n');
printlog(dataFile,'model training time: \t%.3f\n',trainingOnlyTime);
printlog(dataFile, ['model saved to: ' patterns_dir 'locpatternsdata_' num2str(runNum) '_' datestr(now,30) '\n']);
cd(code_dir)
end
