%%so now this would be all the commands you would want to do ONLY for
%%fmri session
%first these are all the session numbers

SUBJECT = 8; %experimental subject number
prev = 0; %if today's date (0) or previous date (1)
scanNow = 1; %if using triggers (1)
runNum = 1; %what number subject they are today

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

% 1: SCOUT
% 2: MPRAGE
% 3: AP Scan
% 4: PA Scan
% 5: Example functional
% 6: LOCALIZER
% 7: RECALL 1
% 8: MOT 1
% 9: MOT 2
% 10: MOT 3
% 11: RECALL 2

%% RUN MP_RAGE FIRST
%% RUN VARIOUS BEHAVIORAL TASKS
%first MOT_PRACTICE and RECALL PRACTICE
mot_realtime04MB(SUBJECT,MOT_PRACTICE2, [],0,scanNow); %will move automatically into RECALL_PRACTICE
%then start RSVP task5
%% SCAN_PREP: instructions and also 8 seconds
scanNum = 5; % here run process Nifti but look at outputs--make sure magnitude bet is okay and functional bet!
mot_realtime04MB(SUBJECT,SCAN_PREP,[],scanNum,scanNow)
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
mot_realtime04MB(SUBJECT,RSVP2,[],0,scanNow) %will continue until TOCRITERION3
%look for mask and test it

%% LOCALIZER DISPLAY
scanNum = 6;
mot_realtime04MB(SUBJECT,MOT_LOCALIZER,[],scanNum,scanNow);

%% LOCALIZER FILE PROCESS
% number of TR's total: 1376 (should be 688 originally)
scanNum = 6;
crossval = 0;
featureSelect = 1;
LocalizerNiftiFileProcess(SUBJECT,crossval,featureSelect,prev,scanNow,scanNum,MOT_LOCALIZER,runNum)

%% RECALL 1
% number of TR's total: 474
scanNum = 7;
mot_realtime04MB(SUBJECT,RECALL1,[],scanNum,scanNow);

%% MOT RUN 1 DISPLAY
% number of TR's total 452 (should be 226 originally)
scanNum = 8; %new would be 15
mot_realtime04MB(SUBJECT,MOT{1},[],scanNum,scanNow);
%% MOT RUN 1 FILE PROCESS
scanNum = 8;%normally 15;
blockNum = 1;
featureSelect = 1;
RealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{1},blockNum,runNum);

%% MOT RUN 2 DISPLAY
scanNum = 9;
mot_realtime04MB(SUBJECT,MOT{2},[],scanNum,scanNow);
%% MOT RUN 2 FILE PROCESS
scanNum = 9;
featureSelect = 1;
blockNum = 2;
RealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{2},blockNum,runNum);

%% MOT RUN 3 DISPLAY
scanNum = 10;
mot_realtime04MB(SUBJECT,MOT{3},[],scanNum,scanNow);
%% MOT RUN 3 FILE PROCESS
scanNum = 10;
featureSelect = 1;
blockNum = 3;
RealTimeNiftiFileProcess(SUBJECT,featureSelect,prev,scanNow,scanNum,MOT{3},blockNum,runNum);
%% RECALL 2
scanNum = 11;
mot_realtime04MB(SUBJECT,RECALL2,[],scanNum,scanNow);
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

scanNum1 = 13;
scanNum2 = 21;
datevec = { '1-13-17', '1-14-17', '1-14-17', '1-20-17', '1-21-17', '1-22-17', '1-26-17', '1-28-17', '1-30-17', '2-1-17'};
svec = [3 4 5 6 7 8 9 11 12 13];
runvec = ones(1,length(svec));
irun2 = find(svec==5);
runvec(irun2) = 2;
nsub = length(svec);
for s = 1:nsub
    SUBJECT = svec(s);
    date = datevec{s};
    runNum = runvec(s);
    AnatRecallFileProcess(SUBJECT,runNum,scanNum1,RECALL1,date);
    AnatRecallFileProcess(SUBJECT,runNum,scanNum2,RECALL2,date);
end
