base_path = [fileparts(which('mot_realtime01.m')) filesep];
cd(base_path);


SUBJECT = 32;
%all given subjects: 8,12,13,14,15,18,21,22
SVEC = [31];


NUM_TASK_RUNS = 3;
% orientation session
SETUP = 1; % stimulus assignment 
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
MOT_PRACTICE2 = SCAN_PREP + 1; %13
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

%% first practice set
mot_realtime01b(SUBJECT, SETUP, [], 0, 0);

%for testing
%mot_realtime01b(SUBJECT,TOCRITERION1,[],0,0);
% this will continue to train test and practice MOT, then move on to
% MOT_Practice, MOT_PREP
s2 = findMatch(SUBJECT,SVEC);
% then restart mot_realtime01 and go to
%familiarize 2 with new data
if s2 < 0 % if no match
    subjDir = [base_path 'BehavioralData' filesep num2str(SUBJECT) filesep];
    oldfile = [subjDir 'mot_realtime01_subj_' num2str(SUBJECT) '_stimAssignment.mat'];
    newfile = [subjDir 'OLDstimAssignment.mat'];
    system(['mv ' subjDir 'mot_realtime01_subj_' num2str(SUBJECT) '_stimAssignment.mat ' subjDir 'OLDstimAssignment.mat']);
    mot_realtime01(SUBJECT,SETUP,[],0,0); % have mot go to familiarize 2 after making pairs!
else %keep going 
    mot_realtime01b(SUBJECT, FAMILIARIZE2, [], 0, 0,s2); %continue because want to not go through the break

end

% TESTING TAKE OUT
%mot_realtime01b(SUBJECT,RECALL2, 1, 0, 0,s2); %continue because want to not go through the break

%mot_realtime01b(SUBJECT, FAMILIARIZE2, [], 0, 0,s2); %continue because want to not go through the break

%% now train on actual stimulus pairs
mot_realtime01b(SUBJECT, FAMILIARIZE2, [], 0, 0);

%% after scanner, test associates
s2 = findMatch(SUBJECT,SVEC)

mot_realtime01b(SUBJECT,ASSOCIATES, [], 0, 0,s2);