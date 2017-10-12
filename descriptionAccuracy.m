%% look at descriptive ratings
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
RECALL2 = MOT{end} + 1; % post-scan rsvp m emory test
DESCRIPTION = RECALL2 + 1; %26
ASSOCIATES = DESCRIPTION + 1; %27
base_path = [fileparts(which('mot_realtime05.m')) filesep];

SESSION = DESCRIPTION;

%%

%subjectVec = [1 3 4 5 6 8 10 12 13 14 19 21 23 26 29 32];
subjectVec = 35;
for s = 1:length(subjectVec)
    correct = [];
    subjectNum = subjectVec(s);
    fprintf('subject: %i\n', subjectNum);
    behavioral_dir = [base_path 'BehavioralData/' num2str(subjectNum) '/'];
    r = dir(fullfile(behavioral_dir, ['mot_realtime05_' num2str(subjectNum) '_' num2str(SESSION) '*.mat']));
    r = load(fullfile(behavioral_dir,r(end).name));
    ID = r.stim.id(4:end);
    desc = r.stim.description;
    
    fn = dir(fullfile(behavioral_dir, ['mot_realtime05_subj_' num2str(subjectNum) '_stimAssignment.mat']));
    r = load(fullfile(behavioral_dir,fn(end).name));
    for t = 1:20
        fprintf('trial: %i\n', t);
        thisID = ID(t);
        thisPic = r.pics{thisID};
        thisDesc = desc{t+3};
        
        figure(1);
        image(imread([base_path 'stimuli/STIM/middlepairs/' thisPic]))
        fprintf('%s\n', thisDesc);
        correct(t) = input('correct?\n');
    end
    save([behavioral_dir 'descriptionAcc'], 'correct');
end
%%
% count number of incorrect details--just are false in image (not vague but
% wrong)
%subjectVec = [1 3 4 5 6 8 10 12 13 14 19 21 23 26 29 32];
%subjectVec = [23 26 29 32];
for s = 1:length(subjectVec)
    nWrong = [];
    subjectNum = subjectVec(s);
    fprintf('subject: %i\n', subjectNum);
    behavioral_dir = [base_path 'BehavioralData/' num2str(subjectNum) '/'];
    r = dir(fullfile(behavioral_dir, ['mot_realtime05_' num2str(subjectNum) '_' num2str(SESSION) '*.mat']));
    r = load(fullfile(behavioral_dir,r(end).name));
    ID = r.stim.id(4:end);
    desc = r.stim.description;
    
    fn = dir(fullfile(behavioral_dir, ['mot_realtime05_subj_' num2str(subjectNum) '_stimAssignment.mat']));
    r = load(fullfile(behavioral_dir,fn(end).name));
    for t = 1:20
        fprintf('trial: %i\n', t);
        thisID = ID(t);
        thisPic = r.pics{thisID};
        thisDesc = desc{t+3};
        
        figure(1);
        image(imread([base_path 'stimuli/STIM/middlepairs/' thisPic]))
        fprintf('%s\n', thisDesc);
        nWrong(t) = input('number wrong details?\n');
    end
    save([behavioral_dir 'descriptionWrong'], 'nWrong');
end
