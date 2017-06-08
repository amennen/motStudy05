function [patterns, t] = RealTimeNiftiFileProcess(subjectNum,featureSelect,prev,scanNow,scanNum,SESSION,blockNum,runNum) %,rtfeedback)
% function [patterns] = RealTimeMemoryFileProcess(subjectNum,subjectName,runNum,scanNum,scanNow)
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
% - scanNow:      whether data acquired in realtime or previously collected [2/1/0]
%                0: offline data
%                1: realtime data from scanner
%                2: simulated realtime data
% - rtfeedback:  whether feedback delivered in realtime or not [1/0]
%
% OUTPUTS
% - patterns: elapsed time for each iteration of SVM testing
%
% Written by: Megan deBettencourt//Anne Mennen made the version to work
% with nifti files
% Version: 1.0, July 2013

% runNum = 1;
% subjectNum = 1;
% scanNum = 13;
% scanIndex = 1; %which scan they were that day
% SESSION = 21;
% scanNow = 2;
%% check inputs
%check that there is a sufficient number of inputs
% if nargin < 5;error('5 inputs are required: subjectNum, subjectName, runNum, scanNum, scanNow');end
% 
% if ~isnumeric(subjectNum);error('subjectNum must be a number');end
% if ~ischar(subjectName);error('subjectName must be a string');end
% if ~isnumeric(runNum);error('runNum must be a number');end
% if ~isnumeric(scanNum);error('scanNum must be a number - equal to the next motion-corrected scan number');end
% if (scanNow~=1) && (scanNow~=0)&& (scanNow~=2);error('scanNow must be either 2 (if simulated realtime) or 1 (if realtime data acquisition) or 0 (if offline data)');end




%% initialize path prefix for different replyDrive
projectName = 'motStudy04'; %CHANGE BACK!!
biac_dir = '/Data1/packages/BIAC_Matlab_R2014a/';
bxhpath='/opt/BXH/1.11.1/bin/';
fslpath='/opt/fsl/5.0.9/bin/';
dcm2path = '/opt/MRICROGL/2-2016/';
%add necessary package
if ~exist('readmr','file')
    addpath(genpath(biac_dir));
    addpath([biac_dir '/mr/']);
    addpath([biac_dir '/general/'])
end
multipath = '/Data1/code/multibandutils/';
addpath(genpath(multipath));
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
process_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/' 'reg' '/'];
roi_dir = ['/Data1/code/' projectName '/data/'];
code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
lastRunHeader = fullfile(save_dir, ['motRun' num2str(blockNum-1) '/']);
locPatterns_dir = fullfile(save_dir, 'Localizer/');
behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
addpath(genpath(code_dir));
%runNum = 1; %assume first person that day
if ~prev %if getting data today
    subjectName = [datestr(now,5) datestr(now,7) datestr(now,11) num2str(runNum) '_' projectName];
    dicom_dir = ['/Data1/subjects/' datestr(now,10) datestr(now,5) datestr(now,7) '.' subjectName '.' subjectName '/'];
else
    allDates = {'4-5-17'};
    %allDates = {'7-1-2016' '3-26-2016', '3-29-2016', '4-1-2016', '4-27-2016', '4-29-2016', '5-05-2016'};
    subjectName = [datestr(allDates{1},5) datestr(allDates{1},7) datestr(allDates{1},11) num2str(runNum) '_' projectName];
    dicom_dir = ['/Data1/subjects/' datestr(allDates{1},10) datestr(allDates{1},5) datestr(allDates{1},7) '.' subjectName '.' subjectName '/'];
end
%check that dicom_dir exists
assert(logical(exist(dicom_dir,'dir')));
fprintf('fMRI files being read from: %s\n',dicom_dir);
%check that the fMRI dicom files do NOT exist (if real-time)
if scanNow == 1
    [testFile testFileName] = GetSpecificFMRIFile(dicom_dir,scanNum,1);
    if exist([dicom_dir testFileName],'file');
        reply = input('Files with this scan number already exist. Do you want to continue? Y/N [N]: ', 's');
        if isempty(reply)
            reply = 'N';
        end
        if ~(strcmp(reply,'Y') || strcmp(reply,'y'))
            return
        end
    end
     if exist( fullfile(behavioral_dir, ['SessionInfo' '_' num2str(SESSION) '.mat']), 'file')
        reply = input('Delete previous session information? Y/N [N]: ', 's');
        if isempty(reply)
            reply = 'N';
        end
        if (strcmp(reply,'Y') || strcmp(reply,'y'))
            unix(sprintf('rm %sSessionInfo_%s.mat', behavioral_dir, num2str(SESSION))); %remove if there's already a file for that person
        end
    end
    
end

classOutputDir = fullfile(runHeader, 'classOutput/');
if ~exist(classOutputDir, 'dir')
    mkdir(classOutputDir);
end
cd(runHeader); % go to run header so that we're saving all the nifti files here
%get ROI
roi_name = 'retrieval';
temp = load(fullfile(process_dir,[roi_name '_mask.mat'])); %should be stretched_brain
roi = logical(temp.mask_brain);
assert(exist('roi','var')==1);
roiDims = size(roi);
roiInds = find(roi);


%load trained model
allfn = dir([locPatterns_dir 'loctrainedModel_' num2str(runNum) '*']); %
%take the last model saved
load(fullfile(locPatterns_dir, allfn(end).name));
%fname = findNewestFile(ppt_dir2,fullfile(ppt_dir2, ['mot_realtime01_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
%load(fname);
fprintf('\n*********************************************\n');
fprintf(['* Loaded ' allfn(end).name '\n']);
%load localizer run's standard deviation and voxelInfo
allLast = dir([locPatterns_dir 'locpatternsdata_' '*']);
loc = load(fullfile(locPatterns_dir, allLast(end).name));
if SESSION == 21 %changed! this should be MOT_RUN1!!
   patterns.lastStd = loc.patterns.runStd;
   %patterns.lastStd = .1; %for testing
else
    allLast = dir([lastRunHeader 'motpatternsdata_' num2str(SESSION-1) '*']);
    last = load(fullfile(lastRunHeader,allLast(end).name));
    patterns.lastStd = last.patterns.runStd;
end
%put this back in!!
fprintf('\n*********************************************\n');
fprintf(['* Loaded ' 'last run''s Std in ' allLast(end).name '\n']);

%load this run's regressors and information (do this after load so loading
%doesn't overwrite)
[newpattern t] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir);
patterns.regressor = newpattern.regressor;
%% Boilerplate

seed = sum(100*clock); %get random seed
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));%set seed
%initialize system time calls
GetSecs;

%% preprocessing parameters
FWHM = 5;
cutoff = 160;
TR = 1; % changed here back!
shiftTR = 4/TR; % this will now be 4 TR's (was 2 TRs for a 2 s TR)
%% Block Sequence

patterns.nTRs = size(patterns.regressor.twoCond,2); %already has first 20 removed
patterns.firstTestTR = find(patterns.regressor.twoCond(1,:)+patterns.regressor.twoCond(2,:),1,'first') ; %(because took out first 10)
patterns.fileAvail = zeros(1,patterns.nTRs);
patterns.newFile = cell(1,patterns.nTRs);
patterns.timeRead = cell(1,patterns.nTRs);
patterns.fileload = NaN(1,patterns.nTRs);
patterns.raw = nan(patterns.nTRs,numel(roiInds));
patterns.raw_sm = nan(patterns.nTRs,numel(roiInds));
patterns.raw_sm_filt = nan(patterns.nTRs,numel(roiInds));
patterns.raw_sm_z = nan(patterns.nTRs,numel(roiInds));
patterns.categsep = nan(1,patterns.nTRs);
patterns.realtimeMean = nan(1,numel(roiInds));
patterns.realtimeStd = nan(1,numel(roiInds));
patterns.realtimeY = nan(1,numel(roiInds));
patterns.realtimeLastMean = nan(1,numel(roiInds));
patterns.realtimeLastStd = nan(1,numel(roiInds));
patterns.realtimeLastY = nan(1,numel(roiInds));
patterns.realtimeVar = nan(1,numel(roiInds));
patterns.predict = nan(1,patterns.nTRs);
patterns.activations = nan(numel(patterns.regressor.twoCond(:,1)),patterns.nTRs);
patterns.block = (SESSION - 19)*ones(1,patterns.nTRs);
patterns.outputClass = zeros(1,patterns.nTRs);
%% Output Files Setup

% open and set-up output file
dataFile = fullfile(runHeader, 'fileprocessing.txt');
printlog(dataFile,'\n*********************************************\n');
printlog(dataFile,'* Induce MOT Training: Curtain Practice Runs v.1.0\n');
printlog(dataFile,['* Date/Time: ' datestr(now,0) '\n']);
printlog(dataFile,['* Seed: ' num2str(seed) '\n']);
printlog(dataFile,['* Subject Number: ' num2str(subjectNum) '\n']);
printlog(dataFile,['* Subject Name: ' subjectName '\n']);
printlog(dataFile,['* Run Number: ' num2str(blockNum) '\n']);
printlog(dataFile,['* Real-Time Data: ' num2str(scanNow) '\n']);
printlog(dataFile,'*********************************************\n\n');

%% Start Experiment

% prepare for trial sequence
printlog(dataFile,'Block\tTR\tFileAvail\tOutputClass\tCategSep\n'); 
%% acquiring files
remove = 20/TR;
zscoreNew = 1;
useHistory = 1;
firstBlockTRs = 128/TR; %total number of TRs to take for standard deviation of last run
for iTrial = 1:patterns.nTRs % the first 10 TRs have been taken out to detrend
    
    tstart(iTrial) = tic;
    zscoreLen = double(iTrial);
    zscoreLen1 = double(iTrial - 1);
    zscoreConst = 1.0/zscoreLen;
    zscoreConst1 = 1.0/zscoreLen1;
    % increase count of TRs
    %increase the count of TR pulses
    thisTR = iTrial + remove; %account for taking out TRs
    while ~patterns.fileAvail(iTrial)
        [patterns.fileAvail(iTrial) patterns.newFile{iTrial}] = GetSpecificFMRIFile(dicom_dir,scanNum,thisTR);
        timing.fileAppear(iTrial) = toc(tstart(iTrial));
    end
    
    %if desired file is recognized, pause for 200ms to complete transfer
    if scanNow==1 && patterns.fileAvail(iTrial)
        pause(.2);
    end
    
    % if file available, load it
    if (patterns.fileAvail(iTrial))
        %[newVol patterns.timeRead{iTrial}] = ReadFile([dicom_dir patterns.newFile{iTrial}],imgmat,roi); % NTB: only reads top file
        t0 = GetSecs;
        niftiname = sprintf('nifti%3.3i', thisTR);
        
        unix(sprintf('%sdcm2niix %s -f %s -o %s -s y %s%s',dcm2path,dicom_dir,niftiname,runHeader,dicom_dir,patterns.newFile{iTrial}))
        t1 = GetSecs;
        unix(sprintf('%smcflirt -in %s.nii -reffile %sexfunc_re.nii',fslpath,niftiname,process_dir))
        t2 = GetSecs;
        moco = t2-t1;
        
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
        
    end
    %smooth files
    patterns.raw_sm(iTrial,:) = SmoothRealTime(patterns.raw(iTrial,:),roiDims,roiInds,FWHM);
    
    %filter, then calculate mean and standard deviation
    if iTrial == (patterns.firstTestTR-1)
        
        patterns.raw_sm_filt(1:iTrial,:) = HighPassBetweenRuns(patterns.raw_sm(1:iTrial,:),TR,cutoff);
        patterns.btwnrunsfiltered(1:iTrial) = 1;
        
        patterns.realtimeMean(1,:) = mean(patterns.raw_sm_filt(1:iTrial,:),1);
        patterns.realtimeY(1,:) = mean(patterns.raw_sm_filt(1:iTrial,:).^2,1);
        %make sure to use population standard deviation, divide by N
        patterns.realtimeStd(1,:) = std(patterns.raw_sm_filt(1:iTrial,:),1,1);
        
        patterns.realtimeVar(1,:) = patterns.realtimeStd(1,:).^2;
    end
    
    if iTrial > (patterns.firstTestTR - 1)

        %filter, then calculate mean and standard deviation
        patterns.raw_sm_filt(iTrial,:) = HighPassRealTime(patterns.raw_sm(1:iTrial,:),TR,cutoff);
        
        patterns.realtimeMean(1,:) = mean(patterns.raw_sm_filt(1:iTrial,:),1);
        patterns.realtimeY(1,:) = mean(patterns.raw_sm_filt(1:iTrial,:).^2,1);
        patterns.realtimeStd(1,:) = std(patterns.raw_sm_filt(1:iTrial,:),1,1); %flad to use N instead of N-1
        patterns.realtimeVar(1,:) = patterns.realtimeStd(1,:).^2;
        
        %now zscore
        
        %record last history
        patterns.realtimeLastMean(1,:) = patterns.realtimeMean(1,:);
        patterns.realtimeLastY(1,:) = patterns.realtimeY(1,:);
        patterns.realtimeLastVar(1,:) = patterns.realtimeVar(1,:);
        %update mean
        patterns.realtimeMean(1,:) = (patterns.realtimeMean(1,:).*zscoreLen1 + patterns.raw_sm_filt(iTrial,:)).*zscoreConst;
        %update y = E(X^2)
        patterns.realtimeY(1,:) = (patterns.realtimeY(1,:).*zscoreLen1+ patterns.raw_sm_filt(iTrial,:).^2).*zscoreConst;
        %update var
        if useHistory
            patterns.realtimeVar(1,:) = patterns.realtimeLastVar(1,:) ...
                + patterns.realtimeLastMean(1,:).^2 - patterns.realtimeMean(1,:).^2 ...
                + patterns.realtimeY(1,:) - patterns.realtimeLastY(1,:);
        else
            % update var
            patterns.realtimeVar(1,:) = patterns.realtimeVar(1,:) - patterns.realtimeMean(1,:).^2 ...
                + ((patterns.realtimeMean(1,:).*zscoreLen - patterns.raw_sm_filt(iTrial,:)).*zscoreConst1).^2 ...
                + (patterns.raw_sm_filt(iTrial,:).^2 - patterns.realtimeY(1,:)).*zscoreConst1;
        end
        if iTrial > firstBlockTRs
            patterns.raw_sm_filt_z(iTrial,:) = (patterns.raw_sm_filt(iTrial,:) - patterns.realtimeMean(1,:))./patterns.realtimeStd(1,:);
        else
            
            patterns.raw_sm_filt_z(iTrial,:) = (patterns.raw_sm_filt(iTrial,:) - patterns.realtimeMean(1,:))./patterns.lastStd(1,:);
        end
        % now test if it's when we want to
        if any(patterns.regressor.twoCond(:,iTrial)) || any(patterns.regressor.twoCond(:,iTrial-(shiftTR+2))) %go a little extra
            if featureSelect
                goodVox = loc.patterns.sigVox;
                [patterns.predict(iTrial), patterns.activations(1:2,iTrial)] = predict_ridge(patterns.raw_sm_filt_z(iTrial,goodVox),trainedModel);
            else
                [patterns.predict(iTrial), patterns.activations(1:2,iTrial)] = predict_ridge(patterns.raw_sm_filt_z(iTrial,:),trainedModel);
                
            end
            patterns.categsep(iTrial) = patterns.activations(1,iTrial) - patterns.activations(2,iTrial);
            classOutput = patterns.categsep(iTrial);
            if scanNow ==1
                save(fullfile(classOutputDir, ['vol_' num2str(thisTR)]), 'classOutput');
                patterns.outputClass(iTrial) = 1;
            end
        end
    end
    % print trial results
    printlog(dataFile,'%d\t%d\t%d\t\t%d\t\t%6.3f\n',patterns.block(iTrial),thisTR,patterns.fileAvail(iTrial),patterns.outputClass(iTrial),patterns.categsep(iTrial) );
end
patterns.runStd = std(patterns.raw_sm_filt,[],1); %std dev across all volumes per voxel
save([runHeader 'motpatternsdata_' num2str(SESSION) '_' datestr(now,30)],'patterns','timing', 't');
cd(code_dir); % go back to original area
end


