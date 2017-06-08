base_path = [fileparts(which('mot_realtime04MB.m')) filesep];
cd(base_path);

SUBJECT = 8;


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

%% first practice set
mot_realtime04MB(SUBJECT, SETUP, [], 0, 0);

% this will continue to train test and practice MOT, then move on to
% MOT_Practice, MOT_PREP
%mot_realtime02(SUBJECT,MOT_PRACTICE,[],0,0);
mot_realtime04MB(SUBJECT, FAMILIARIZE2, [], 0, 0); %continue because want to not go through the break

%% DAY TWO

%% Refresh on day 2
mot_realtime04MB(SUBJECT, STIM_REFRESH, [], 0, 0);

%% after scanner, test associates and descriptions

mot_realtime04MB(SUBJECT,DESCRIPTION, [], 0, 0); 
mot_realtime04MB(SUBJECT,ASSOCIATES, [], 0, 0); 

%% now convert recog to cell
subjects = [8];
for s = 1:length(subjects)
    behavioral_dir = ['BehavioralData/' num2str(subjects(s)) '/']
    r = dir(fullfile(behavioral_dir, ['_RECOG' '*.mat']));
    r = load(fullfile(behavioral_dir,r(end).name));
    trials = table2cell(r.datastruct.trials);
    stimID = cell2mat(trials(:,8));
    cond = cell2mat(trials(:,9));
    acc = cell2mat(trials(:,11));
    rt = cell2mat(trials(:,13));
    cresp = cell2mat(trials(:,22));
    save([behavioral_dir '/' 'recogcell.mat'],'stimID', 'cond', 'acc', 'rt', 'cresp')
end