% get regressors (and selectors if want to cross validate)
% inputs are subject number, which session number, if want to cross
% validate or not
function [patterns trials stimOrder hardSpeed acc rt] = GetSessionInfo(subjNum,SESSION,varargin)
%subjNum = 1;
%SESSION = 18; %for localizer
if SESSION == 18
    crossval = 1;
else
    crossval = 0;
end

if IsLinux
    behav_dir = '/Data1/code/motStudy02/BehavioralData/';
elseif ismac
    behav_dir = '/Users/amennen/Documents/Norman/MOT/mot_study/Temp_participantData/';
    cd('/Users/amennen/Documents/Norman/MOT/mot_study/');
end
LOC = 18;
MOT = [18 20:22]; %(can change the way the files are named in the future)
RECALL = [19 23];
MOT_PREP = 5;
MAX_SPEED = 30;
hardSpeed = nan;
acc = nan;
rt = NaN;

subjectFolder = [behav_dir num2str(subjNum)];
if ~isempty(varargin)
    N_TRS_LOC = cell2mat(varargin);
else
    N_TRS_LOC = 15; %set to all if don't specify
end
NCOND = 4;

if ismember(SESSION,MOT)
    fn = dir(fullfile(subjectFolder,['EK' num2str(SESSION) '_DOT*mat']));
else
    fn = dir(fullfile(subjectFolder,['EK' num2str(SESSION) '_SUB*mat']));
end
data = load(fullfile(subjectFolder,fn(end).name));
fn_stim = dir(fullfile(subjectFolder,['mot_realtime01_' num2str(subjNum) '_' num2str(SESSION) '*mat']));
s = load(fullfile(subjectFolder,fn_stim(end).name));
nTRs = s.config.nTRs.perBlock + 5; %includes 5 seconds at the end

% get hard dot speed
fileSpeed = dir([subjectFolder '/' 'mot_realtime01_' num2str(subjNum) '_' num2str(MOT_PREP)  '*.mat']);
if ~isempty(fileSpeed)
    matlabOpenFile = [subjectFolder '/' fileSpeed(end).name];
    lastRun = load(matlabOpenFile);
    hardSpeed = MAX_SPEED - lastRun.stim.tGuess(end);
end

VARIATIONS_MAT = zeros(NCOND,nTRs); %regressor with all four conditions
SELECTOR_XVAL = zeros(1,nTRs); %which TRs are for training and testing for cross-validation

conditions = data.datastruct.trials.cond;
stimOrderAll = data.datastruct.trials.stim_id;
nTrials = length(conditions);
TH = find(conditions==1); %for recall this is fast dot motion
TE = find(conditions==2); %for recall this is slow dot motion
LH = find(conditions==3); %for recall this is omit trials
LE = find(conditions==4); %for recall there's no condition 4
trials.hard = TH;
trials.easy = TE;
trials.lure = [LH LE];
[~,stimOrder.hard] = sort(stimOrderAll(TH));
[~,stimOrder.easy] = sort(stimOrderAll(TE));
[~,stimOrder.lure] = sort(stimOrderAll(LH));
% get the conditions by trial here
if ismember(SESSION,LOC) %MOT is only 1 condition
    accuracy = data.datastruct.trials.acc;
    acc.hard = accuracy(TH);
    acc.easy = accuracy(TE);
    acc.lure = accuracy(LH);
    
    reaction = data.datastruct.trials.rt;
    rt.hard = reaction(TH);
    rt.easy = reaction(TE);
    rt.lure = reaction(LH);    
end
% now we have to match these to TRs to get the actual regressors
if ismember(SESSION,MOT)
    iTR.start = convertTR(s.timing.trig.wait,s.timing.actualOnsets.motion(1,:),s.stim.TRlength);
    trialDur = s.timing.plannedOnsets.probe(1) - s.timing.plannedOnsets.motion(1,1) +4;
else
    iTR.start = convertTR(s.timing.trig.wait,s.timing.actualOnsets.prompt,s.stim.TRlength);
    trialDur = s.timing.plannedOnsets.record(1) - s.timing.plannedOnsets.prompt(1); %try full recording TR just to see what it says!
    %trialDur = s.timing.plannedOnsets.vis(1) - s.timing.plannedOnsets.prompt(1);
end
trialDurTR = (trialDur/s.stim.TRlength) - 1; %20s/2 = 10 - 1 = 9 TRs
if SESSION == 18 && N_TRS_LOC > 0 %shift over a little bit more
    trialDurTR = N_TRS_LOC - 1;
end

iTR.TH = iTR.start(TH);
iTR.TE = iTR.start(TE);
iTR.LH = iTR.start(LH);
iTR.LE = iTR.start(LE);

%make matrix of target hard, target easy, lure hard, lure easy
for i=1:length(iTR.TH);
    VARIATIONS_MAT(1,iTR.TH(i):iTR.TH(i) + trialDurTR) = 1;
    VARIATIONS_MAT(2,iTR.TE(i):iTR.TE(i) + trialDurTR) = 1;
    SELECTOR_XVAL(iTR.TH(i):iTR.TH(i) + trialDurTR)= i;
    SELECTOR_XVAL(iTR.TE(i):iTR.TE(i) + trialDurTR)= i;
end
for i=1:length(iTR.LH)
    VARIATIONS_MAT(3,iTR.LH(i):iTR.LH(i)+trialDurTR) = 1;
    SELECTOR_XVAL(iTR.LH(i):iTR.LH(i) + trialDurTR)= i;
end
for i=1:length(iTR.LE)
    VARIATIONS_MAT(4,iTR.LE(i):iTR.LE(i)+trialDurTR) = 1;
    SELECTOR_XVAL(iTR.LE(i):iTR.LE(i) + trialDurTR)= i;
end


targets = sum(VARIATIONS_MAT(1:2,:));
lures = sum(VARIATIONS_MAT(3:4,:));
patterns.regressor.allCond = VARIATIONS_MAT(:,11:end);
REGRESSORS = [targets;lures];
patterns.regressor.twoCond = REGRESSORS(:,11:end); %get rid of first 10 TRs

patterns.selector.xval = SELECTOR_XVAL(11:end);
% make the separate selectors
if crossval
    nIterations = length(iTR.TH);
    for j = 1:nIterations
        allnonzero = find(patterns.selector.xval);
        thisIndex = find(patterns.selector.xval==j);
        temp = patterns.selector.xval;
        temp(allnonzero) = 1;
        temp(thisIndex) = 2;
        patterns.selector.allxval(j,:) = temp;
    end
end
end
