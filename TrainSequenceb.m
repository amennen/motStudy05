base_path = [fileparts(which('mot_realtime05.m')) filesep];
cd(base_path);


SUBJECT = 102;
subjDir = [base_path 'BehavioralData' filesep num2str(SUBJECT) filesep];

%all given subjects: 8,12,13,14,15,18,21,22
SVEC = [1 3:6 8:10];


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
mot_realtime05b(SUBJECT, SETUP, [], 0, 0);
%for testing
%mot_realtime01b(SUBJECT,TOCRITERION1,[],0,0);
% this will continue to train test and practice MOT, then move on to
% MOT_Practice, MOT_PREP
s2 = findMatch(SUBJECT,SVEC);
% then restart mot_realtime01 and go to
%familiarize 2 with new data
if s2 < 0 % if no match then do regular RT
    % don't need to go to setup--should just be able to continune because
    % would have made all stimuli
    originalFomat = 1;
    save([subjDir 'matchSubj.mat'],'s2','SVEC', 'originalFormat');
    mot_realtime05(SUBJECT,FAMILIARIZE2,[],0,0); % have mot go to familiarize 2 after making pairs!
else %keep going 
    oldfile = [subjDir 'mot_realtime05_subj_' num2str(SUBJECT) '_stimAssignment.mat'];
    newfile = [subjDir 'OLDstimAssignment.mat'];
    system(['cp ' subjDir 'mot_realtime05_subj_' num2str(SUBJECT) '_stimAssignment.mat ' subjDir 'OLDstimAssignment.mat']);
    originalFormat = 0;
    save([subjDir 'matchSubj.mat'],'s2','SVEC', 'originalFormat');
    mot_realtime05b(SUBJECT, FAMILIARIZE2, [], 0, 0,s2); %continue because want to not go through the break
end


%% REFRESH
%% Refresh on day 2
load([subjDir 'matchSubj.mat'])
mot_realtime05b(SUBJECT, STIM_REFRESH, 1, 0, 0,s2);

%% after scanner, test associates and descriptions
load([subjDir 'matchSubj.mat'])
mot_realtime05b(SUBJECT,DESCRIPTION, 1, 0, 0,s2); 
mot_realtime05b(SUBJECT,ASSOCIATES, 1, 0, 0,s2); 

