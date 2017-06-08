
function [patterns, t] = RecallFileProcess(subjectNum,runNum,scanNum,SESSION,date,featureSelect,saveNew,trainedFile) %,rtfeedback)
% function [patterns] = RealTimeMemoryFileProcess(subjectNum,subjectName,runNum,scanNum,rtData)
%
% this function describes the file processing procedure for the realtime
% fMRI attentional training experiment
%
%
% REQUIRED INPUTS:
% - subjectNum:  participant number [any integer]
%                if subjectNum = 0, no information will be saved
% - subjectName: ntblab subject naming convention [MMDDYY#_REALTIME02]
% - runNum:      run number [any integer]
% - scanNum:     whether collecting scanNum data [scannumber if yes/0 if
% not] (MOT runs are 13, 15, 17)
% - rtData:      whether data acquired in realtime or previously collected [2/1/0]
%                0: offline data
%                1: realtime data from scanner
%                2: simulated realtime data
% - rtfeedback:  whether feedback delivered in realtime or not [1/0]
%
% OUTPUTS
% - patterns: elapsed time for each iteration of SVM testing
%
% Written by: Megan deBettencourt
% Version: 1.0, July 2013

% runNum = 1;
% subjectNum = 1;
% scanNum = 13;
% scanIndex = 1; %which scan they were that day
% SESSION = 21;
% rtData = 2;
%% check inputs
%check that there is a sufficient number of inputs
% if nargin < 5;error('5 inputs are required: subjectNum, subjectName, runNum, scanNum, rtData');end
% 
% if ~isnumeric(subjectNum);error('subjectNum must be a number');end
% if ~ischar(subjectName);error('subjectName must be a string');end
% if ~isnumeric(runNum);error('runNum must be a number');end
% if ~isnumeric(scanNum);error('scanNum must be a number - equal to the next motion-corrected scan number');end
% if (rtData~=1) && (rtData~=0)&& (rtData~=2);error('rtData must be either 2 (if simulated realtime) or 1 (if realtime data acquisition) or 0 (if offline data)');end




%% initialize path prefix for different replyDrive

projectName = 'motStudy03';
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
process_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/' 'reg' '/'];
code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
locPatterns_dir = fullfile(save_dir, 'Localizer/');
behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
addpath(genpath(code_dir));
subjectName = [datestr(date,5) datestr(date,7) datestr(date,11) num2str(runNum) '_' projectName];
dicom_dir = ['/Data1/subjects/' datestr(date,10) datestr(date,5) datestr(date,7) '.' subjectName '.' subjectName '/'];

%check that dicom_dir exists
assert(logical(exist(dicom_dir,'dir')));
fprintf('fMRI files being read from: %s\n',dicom_dir);

runHeader = [save_dir '/recallrun' num2str(runNum)];
classOutputDir = [runHeader '/classOutput/'];

%get ROI
imgmat = 64; %image matrix size
anat_mask = 1;
roi_name = 'retrieval';
stretched = 0;
if anat_mask
    if stretched
        temp = load(fullfile(process_dir,[roi_name '_anat_mask' '.mat'])); %should be stretched_brain
        roi = logical(temp.stretched_brain);
    else
        temp = load(fullfile(process_dir,[roi_name '_anat_mask_orig' '.mat']));
        roi = logical(temp.mask_brain);
    end
    assert(exist('roi','var')==1);
    roiDims = size(roi);
    roiInds = find(roi);
end


%load trained model:original way
trainedFile = 'loctrainedModel'; %uncomment if doing original way
fname = findNewestFile(locPatterns_dir,[locPatterns_dir '/' trainedFile '*']);
load(fname);

%allfn = dir([locPatterns_dir '/loctrainedModel' '*']);
%take the last model saved
%load(fullfile(locPatterns_dir, allfn(end).name));
fprintf('\n*********************************************\n');
fprintf(['* Loaded ' fname '\n']);

%load featureSelected voxels
allLast = findNewestFile(locPatterns_dir,[locPatterns_dir 'locpatternsdata_' '*']);
loc = load(allLast);

%load this run's regressors and information (do this after load so loading
%doesn't overwrite)
[newpattern, t] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir); %check that this works, make sure it gets if conditions are REALTIME or OMIT!
patterns.regressor = newpattern.regressor;
%% Boilerplate

seed = sum(100*clock); %get random seed
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));%set seed

%initialize system time calls
GetSecs;

%% preprocessing parameters
FWHM = 5;
cutoff = 160;
rtData = 2;
TR = 2;
shiftTR = 2;
%% Block Sequence

patterns.nTRs = size(patterns.regressor.twoCond,2); %already has first 10 removed
patterns.firstTestTR = find(patterns.regressor.twoCond(1,:)+patterns.regressor.twoCond(2,:),1,'first') ; %(because took out first 10)
patterns.fileAvail = zeros(1,patterns.nTRs);
patterns.newFile = cell(1,patterns.nTRs);
patterns.timeRead = cell(1,patterns.nTRs);
patterns.fileload = NaN(1,patterns.nTRs);
patterns.raw = nan(patterns.nTRs,numel(roiInds));
patterns.raw_sm = nan(patterns.nTRs,numel(roiInds));
patterns.raw_sm_filt = nan(patterns.nTRs,numel(roiInds));
patterns.raw_sm_z = nan(patterns.nTRs,numel(roiInds));
patterns.categSep = nan(1,patterns.nTRs);
patterns.realtimeMean = nan(1,numel(roiInds));
patterns.realtimeStd = nan(1,numel(roiInds));
patterns.realtimeY = nan(1,numel(roiInds));
patterns.realtimeLastMean = nan(1,numel(roiInds));
patterns.realtimeLastStd = nan(1,numel(roiInds));
patterns.realtimeLastY = nan(1,numel(roiInds));
patterns.realtimeVar = nan(1,numel(roiInds));
patterns.predict = nan(1,patterns.nTRs);
patterns.activations = nan(numel(patterns.regressor.twoCond(:,1)),patterns.nTRs);
normalization_const = zeros(1,numel(roiInds));
if SESSION == 20 %whatever recall 1 is
    patterns.block = ones(1,patterns.nTRs);
else
    patterns.block = 2*ones(1,patterns.nTRs);
end
%Output Files Setup
%% 

% open and set-up output file
dataFile = fopen([save_dir 'fileprocessing.txt'],'a');
fprintf(dataFile,'\n*********************************************\n');
fprintf(dataFile,'* Induce MOT Training: Curtain Practice Runs v.1.0\n');
fprintf(dataFile,['* Date/Time: ' datestr(now,0) '\n']);
fprintf(dataFile,['* Seed: ' num2str(seed) '\n']);
fprintf(dataFile,['* Subject Number: ' num2str(subjectNum) '\n']);
fprintf(dataFile,['* Subject Name: ' subjectName '\n']);
fprintf(dataFile,['* Run Number: ' num2str(runNum) '\n']);
fprintf(dataFile,['* Real-Time Data: ' num2str(rtData) '\n']);
fprintf(dataFile,'*********************************************\n\n');

% print header to command window
fprintf('\n*********************************************\n');
fprintf('* MOT Study v.1.0\n');
fprintf(['* Date/Time: ' datestr(now,0) '\n']);
fprintf(['* Seed: ' num2str(seed) '\n']);
fprintf(['* Subject Number: ' num2str(subjectNum) '\n']);
fprintf(['* Subject Name: ' subjectName '\n']);
fprintf(['* Run Number: ' num2str(runNum) '\n']);
fprintf(['* Real-Time Data: ' num2str(rtData) '\n']);
fprintf('*********************************************\n\n');
%% Start Experiment

% prepare for trial sequence
fprintf(dataFile,'run\tblock\tTR\tloaded\n');
fprintf('run\tblock\tTR\tloaded\n');

%% acquiring files

zscoreNew = 1;
useHistory = 1;
for iTrial = 1:patterns.nTRs % the first 10 TRs have been taken out to detrend
    
    tstart(iTrial) = tic;
    zscoreLen = double(iTrial);
    zscoreLen1 = double(iTrial - 1);
    zscoreConst = 1.0/zscoreLen;
    zscoreConst1 = 1.0/zscoreLen1;
    % increase count of TRs
    %increase the count of TR pulses
    thisTR = iTrial + 10; %account for taking out TRs
    while ~patterns.fileAvail(iTrial)
        [patterns.fileAvail(iTrial) patterns.newFile{iTrial}] = GetSpecificFMRIFile(dicom_dir,scanNum,thisTR);
        %timing.fileAppear(iTrial) = toc(tstart(iTrial));
    end
    
    %if desired file is recognized, pause for 200ms to complete transfer
    if rtData==1 || exist('reply','var')
        pause(.2);
    end
    
    % if file available, load it
    if (patterns.fileAvail(iTrial))
        [newVol patterns.timeRead{iTrial}] = ReadFile([dicom_dir patterns.newFile{iTrial}],imgmat,roi); % NTB: only reads top file
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
    fprintf(dataFile,'%d\t%d\t%d\t%d\n',runNum,patterns.block(iTrial),thisTR,patterns.fileAvail(iTrial));
    fprintf('%d\t%d\t%d\t%d\n',runNum,patterns.block(iTrial),thisTR,patterns.fileAvail(iTrial));
end
%% pre-process
fprintf(dataFile,'\n*********************************************\n');
fprintf(dataFile,'beginning recall preprocessing...\n');
fprintf('\n*********************************************\n');
fprintf('beginning recall preprocessing...\n');

preprocStart = tic; %start timing 

% detrend
patterns.raw_sm_filt = HighPassBetweenRuns(patterns.raw_sm,TR,cutoff);

% z-score
patterns.runMean = mean(patterns.raw_sm_filt,1);  %mean across all volumes per voxel
patterns.runStd = std(patterns.raw_sm_filt,[],1); %std dev across all volumes per voxel
patterns.raw_sm_filt_z = (patterns.raw_sm_filt - repmat(patterns.runMean,size(patterns.raw_sm_filt,1),1))./repmat(patterns.runStd,size(patterns.raw_sm_filt,1),1);
preprocTime = toc(preprocStart);  %end timing
% print training timing and results

%save([save_dir 'recallpreprocpatternsdata_' num2str(runNum) '_' datestr(now,30)],'patterns');
fprintf(dataFile,'data preprocessing time: \t%.3f\n',preprocTime);
fprintf('data preprocessing time: \t%.3f\n',preprocTime);

%% test each TR
fprintf(dataFile,'\n*********************************************\n');
fprintf(dataFile,'beginning recall model testing...\n');
fprintf('\n*********************************************\n');
fprintf('beginning recall model testing...\n');

% prepare for trial sequence
fprintf(dataFile,'run\tblock\tTR\tloaded\n');
fprintf('run\tblock\tTR\tloaded\n');
    for iTrial = 1+shiftTR:patterns.nTRs
       % if any(patterns.regressor.twoCond(:,iTrial-shiftTR)) %take out so
       % it will test every TR
        if featureSelect
            %not localizer so use
           % goodVox=sigVox; %use this version if you're inputting file!!
            goodVox = loc.patterns.sigVox;
            [patterns.predict(iTrial), patterns.activations(1:2,iTrial)] = predict_ridge(patterns.raw_sm_filt_z(iTrial,goodVox),trainedModel);
        else
            [patterns.predict(iTrial), patterns.activations(1:2,iTrial)] = predict_ridge(patterns.raw_sm_filt_z(iTrial,:),trainedModel);
            
        end

       patterns.categSep(iTrial) = patterns.activations(1,iTrial) - patterns.activations(2,iTrial); %this is what we want to use as our output date point
       classOutput = patterns.categSep(iTrial);
       % if we want to save
       %save(fullfile(classOutputDir, ['vol_' num2str(thisTR)]), 'classOutput');
        
    % print trial results
    fprintf(dataFile,'%d\t%d\t%d\t%d\n',runNum,patterns.block(iTrial),iTrial,patterns.fileAvail(iTrial));
    fprintf('%d\t%d\t%d\t%d\n',runNum,patterns.block(iTrial),iTrial,patterns.fileAvail(iTrial));
    end
    if saveNew
        save([save_dir 'recallpatternsdata_' num2str(SESSION) '_' datestr(now,30)],'patterns', 't');
    end
