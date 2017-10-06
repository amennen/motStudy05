%%so now this would be all the commands you would want to do ONLY for
%%fmri session

%first these are all the session numbers

SUBJECT = 33; %experimental subject number
prev = 0; %5if today's date (0) or previous date (1)
scanNow = 1; %if using triggers (1)
runNum = 3; %what number subject they are today
base_path = [fileparts(which('mot_realtime05.m')) filesep];
subjDir = [base_path 'BehavioralData' filesep num2str(SUBJECT) filesep];
load([subjDir 'matchSubj.mat'])

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
%% RUN MP_RAGE FIRST
%% RUN VARIOUS BEHAVIORAL TASKS
%first MOT_PRACTICE and RECALL PRACTICE
mot_realtime05b(SUBJECT,MOT_PRACTICE2, [],0,scanNow,s2); %will move automatically into RECALL_PRACTICE
%then start RSVP task5
%% SCAN_PREP: instructions and also 8 seconds
scanNum = 5; % here run process Nifti but look at outputs--make sure magnitude bet is okay and functional bet!
mot_realtime05b(SUBJECT,SCAN_PREP,[],scanNum,scanNow,s2)
% run the AP scan, PA scan, and functional scan here!! then just quit
% qhenever they finish
%% SCAN_PREP FILE PROCESS
% scanNum = 5; %change 
% processNew = 1;
% ProcessNiftiMask(SUBJECT,processNew,prev,scanNum,runNum) %have it so it waits until it finds the file

%% NOW RUN FIELD MAPS WHILE NEXT BEHAVIORAL TASKS (RSVP2,FAMILIARIZE3,TOCRITERION3)
%changed from RSVP2 change back!!
% run after functional scans! as soon as it finishes!
% (could maybe have two windows to process anatomical and other functional
% ones but eh)
mot_realtime05b(SUBJECT,RSVP2,[],0,scanNow,s2) %will continue until TOCRITERION3
%look for mask and test it

%% LOCALIZER DISPLAY
scanNum = 9;
mot_realtime05b(SUBJECT,MOT_LOCALIZER,[],scanNum,scanNow,s2);

%% LOCALIZER FILE PROCESS
% number of TR's total: (should be 688 originally)
scanNum = 9;
crossval = 0;
featureSelect = 1;
LocalizerNiftiFileProcess(SUBJECT,crossval,featureSelect,prev,scanNow,scanNum,MOT_LOCALIZER,runNum)

%% RECALL 1
% number of TR's total: 474
scanNum = 10;
mot_realtime05b(SUBJECT,RECALL1,[],scanNum,scanNow,s2);

%% MOT RUN 1 DISPLAY
% number of TR's total 452 (should be 226 originally)
scanNum = 11; %11new would be 15
mot_realtime05b(SUBJECT,MOT{1},[],scanNum,scanNow,s2);
%% MOT RUN 1 FILE PROCESS
scanNum = 11;%normally 15;
blockNum = 1;
featureSelect = 1;
RealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{1},blockNum,runNum);

%% MOT RUN 2 DISPLAY
scanNum = 12;
mot_realtime05b(SUBJECT,MOT{2},[],scanNum,scanNow,s2);
%% MOT RUN 2 FILE PROCESS
scanNum = 12;
featureSelect = 1;
blockNum = 2;
RealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{2},blockNum,runNum);

%% MOT RUN 3 DISPLAY
scanNum = 13;
mot_realtime05b(SUBJECT,MOT{3},[],scanNum,scanNow,s2);
%% MOT RUN 3 FILE PROCESS
scanNum = 13;
featureSelect = 1;
blockNum = 3;
RealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{3},blockNum,runNum);
%% RECALL 2
scanNum = 14;
mot_realtime05b(SUBJECT,RECALL2,[],scanNum,scanNow,s2);
%% ANALYZE RECALL DATA
% do for recall 1 and recall 2
makeFile = 1;
scanNum1 = 7;
scanNum2 = 11;
featureSelect = 1;
if prev
   date = '2-1-17';
end
RecallNiftiFileProcess(SUBJECT,runNum,scanNum1,RECALL1,date,featureSelect,makeFile,1);
RecallNiftiFileProcess(SUBJECT,runNum,scanNum2,RECALL2,date,featureSelect,makeFile,2);
%% ANALYZE RECALL DATA FOR ANOTHER MASK

scanNum1 = 10;
scanNum2 = 14;
datevec = {'6-15-17', '6-21-17', '6-23-17', '6-25-17', '6-25-17', '6-28-17'};
svec = [1 3:6 8];
runvec = ones(1,length(svec));
irun2 = find(svec==6);
runvec(irun2) = 2;
nsub = length(svec);
for s = 6
    SUBJECT = svec(s);
    if SUBJECT == 1
        scanNum1 = 11;
        scanNum2 = 15;
    else
        scanNum1 = 10;
        scanNum2 = 14;
    end
    date = datevec{s};
    runNum = runvec(s);
    AnatRecallFileProcess(SUBJECT,runNum,scanNum1,RECALL1,date,1);
    AnatRecallFileProcess(SUBJECT,runNum,scanNum2,RECALL2,date,2);
end
