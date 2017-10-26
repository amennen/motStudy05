%%so now this would be all the commands you would want to do ONLY for
%%fmri session
%first these are all the session numbers

SUBJECT = 24; %experimental subject number
prev = 1; %5if today's date (0) or previous date (1)
scanNow = 0; %if using triggers (1)
runNum = 1; %what number subject they are today

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for rerunning code:
% prev = 1;
% scanNow = 0;
% change date inside RealTimeNiftiFileProcess!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SPTB_PATH = ['/Data1/code/SPTBanne'];
addpath(genpath(SPTB_PATH));

NUM_TASK_RUNS = 3;
% orientation session
SETUP = 1; % stimulus assignment 1
FAMILIARIZE = SETUP + 1; % rsvp study learn associates 2
TOCRITERION1 = FAMILIARIZE + 1; % rsvp train to critereon 3
MOT_PRACTICE = TOCRITERION1 + 1;%4
MOT_PREP = MOT_PRACTICE + 1;%5

% day 1
FAMILIARIZE2 = MOT_PREP + 2; % rsvp study learn associates %7
TOCRITERION2 = FAMILIARIZE2 + 1; % rsvp train to critereon
TOCRITERION2_REP = TOCRITERION2 + 1;
RSVP = TOCRITERION2_REP + 1; % rsvp train to critereon

% day 2
STIM_REFRESH = RSVP + 2; %12
SCAN_PREP = STIM_REFRESH + 1; %13
MOT_PRACTICE2 = SCAN_PREP + 1; %14
RECALL_PRACTICE = MOT_PRACTICE2 + 1;
%SCAN_PREP = RECALL_PRACTICE + 1;
RSVP2 = RECALL_PRACTICE + 1; % rsvp train to critereon
FAMILIARIZE3 = RSVP2 + 1; % rsvp study learn associates
TOCRITERION3 = FAMILIARIZE3 + 1; % rsvp train to critereon
MOT_LOCALIZER = TOCRITERION3 + 1; % category classification
RECALL1 = MOT_LOCALIZER + 1;
counter = RECALL1 + 1; MOT = [];
for i=1:NUM_TASK_RUNS
    MOT{i} = counter;
    counter = counter + 1;
end
RECALL2 = MOT{end} + 1; % post-scan rsvp memory test
DESCRIPTION = RECALL2 + 1; %26
ASSOCIATES = DESCRIPTION + 1; %27
%last input is scan number

% 1-4: SCOUT
% 5: MPRAGE
% 6: Example functional
% 7-8: FIELDMAP 
% 9: LOCALIZER
% 10: RECALL 1
% 11: MOT 1
% 12: MOT 2
% 13: MOT 3
% 14: RECALL 2

%% LOCALIZER FILE PROCESS
% number of TR's total: (should be 688 originally)
scanNum = 9;
crossval = 0;
featureSelect = 1;
OLDLocalizerNiftiFileProcess(SUBJECT,crossval,featureSelect,prev,scanNow,scanNum,MOT_LOCALIZER,runNum)

%% MOT RUN 1 FILE PROCESS
scanNum = 11;%normally 15;
blockNum = 1;
featureSelect = 1;
OLDRealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{1},blockNum,runNum);

%% MOT RUN 2 FILE PROCESS
scanNum = 12;
featureSelect = 1;
blockNum = 2;
OLDRealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{2},blockNum,runNum);
%% MOT RUN 3 FILE PROCESS
scanNum = 13;
featureSelect = 1;
blockNum = 3;
OLDRealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{3},blockNum,runNum);
