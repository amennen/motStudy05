%%so now this would be all the commands you would want to do ONLY for
%%fmri session
%first these are all the session numbers

SUBJECT = 32; %experimental subject number
prev = 0; %if today's date (0) or previous date (1)
scanNow = 1; %if using triggers (1)
runNum = 1; %what number subject they are today
svec = [31];
s2 = findMatch(SUBJECT,svec); %enter in s2!!
SPTB_PATH = ['/Data1/code/SPTBanne'];
addpath(genpath(SPTB_PATH));
% if prev

%     allScanNums = [7:2:19];
% else
%     allScanNums = [7 11:2:21];
% end
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
SCAN_PREP = RSVP + 2;
MOT_PRACTICE2 = SCAN_PREP + 1; %12
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
ASSOCIATES = RECALL2 + 1;
%last input is scan number
%scanning numbers should be 21 total3
% 1-4: SCOUT
% 5-8
% 5: MPRAGE
% (6)7: EXFUNCTIONAL
% 8-9: FIELDMAP
% (10)11: LOCALIZER
% (12)13: RECALL 1
% (14)15: MOT 1
% (16)17: MOT 2
% (18)19: MOT 3
% (20)21: RECALL 2

%% RUN MP_RAGE FIRST
%% RUN VARIOUS BEHAVIORAL TASKS
%first MOT_PRACTICE and RECALL PRACTICE
mot_realtime01b(SUBJECT,MOT_PRACTICE2, [],0,scanNow,s2); %will move automatically into RECALL_PRACTICE
%then start RSVP task5
%% SCAN_PREP: instructions and also 8 seconds
scanNum = 7;
mot_realtime01b(SUBJECT,SCAN_PREP,[],scanNum,scanNow,s2)

%% SCAN_PREP FILE PROCESS
scanNum = 7; %change 
processNew = 1;
ProcessMask(SUBJECT,processNew,prev,scanNum,runNum) %have it so it waits until it finds the file

%% NOW RUN FIELD MAPS WHILE NEXT BEHAVIORAL TASKS (RSVP2,FAMILIARIZE3,TOCRITERION3)

mot_realtime01b(SUBJECT,RSVP2,[],0,scanNow,s2) %will continue until TOCRITERION3
%look for mask and test it

%% LOCALIZER DISPLAY
scanNum = 11;
mot_realtime01b(SUBJECT,MOT_LOCALIZER,[],scanNum,scanNow,s2);

%% LOCALIZER FILE PROCESS
scanNum = 11;
crossval = 0;
%crossval = 1;
featureSelect = 1;
LocalizerFileProcess(SUBJECT,crossval,featureSelect,prev,scanNow,scanNum,MOT_LOCALIZER,runNum)

%% RECALL 1
scanNum = 13;
mot_realtime01b(SUBJECT,RECALL1,[],scanNum,scanNow,s2);

%% MOT RUN 1 DISPLAY
scanNum = 15; %new would be 15
mot_realtime01b(SUBJECT,MOT{1},[],scanNum,scanNow,s2);
%% MOT RUN 1 FILE PROCESS
scanNum = 15;
blockNum = 1;
featureSelect = 1;
RealTimeMemoryFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{1},blockNum,runNum);

%% MOT RUN 2 DISPLAY
scanNum = 17;
mot_realtime01b(SUBJECT,MOT{2},[],scanNum,scanNow,s2);
%% MOT RUN 2 FILE PROCESS
scanNum = 17;
featureSelect = 1;
blockNum = 2;
RealTimeMemoryFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{2},blockNum,runNum);

%% MOT RUN 3 DISPLAY
scanNum = 19;
mot_realtime01b(SUBJECT,MOT{3},[],scanNum,scanNow,s2);
%% MOT RUN 3 FILE PROCESS
scanNum = 19;
featureSelect = 1;
blockNum = 3;
RealTimeMemoryFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{3},blockNum,runNum);

%% RECALL 2
scanNum = 21;
mot_realtime01b(SUBJECT,RECALL2,[],scanNum,scanNow,s2);