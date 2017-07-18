% check yoking--created to make sure that things are okay
s1 = 11; %this is the YC subject
s2 = 6;

base_path = [fileparts(which('mot_realtime05.m')) filesep];
s1_dir = fullfile(base_path, 'BehavioralData', num2str(s1));
s2_dir = fullfile(base_path, 'BehavioralData', num2str(s2));

%% check stimulus assignments

fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' 'subj_' num2str(s1) '_stimAssignment'  '*.mat']));
s1_stim = load(fname);

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' 'subj_' num2str(s2) '_stimAssignment'  '*.mat']));
s2_stim = load(fname);

IDX_cues = find(cellfun(@isequal,s2_stim.preparedCues, s1_stim.preparedCues)) %this should be everything except 21-28!
IDX_pics = find(cellfun(@isequal,s2_stim.pics, s1_stim.pics)) %this should be everything except 21-28!

s1_stim.pics(21:28) %these should all be in the 00's
s1_stim.preparedCues(21:28) %these should be only the training words

%% check familiarization matches familiarize2, can also do familiarize 3 which is 16 (2,7,16)

SESSION = STIM_REFRESH; 


fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);

IDs = find(cellfun(@isequal,s1_f.stim.stim, s2_f.stim.stim)) %this should be ALL
IDp = find(cellfun(@isequal,s1_f.stim.picStim, s2_f.stim.picStim)) %this should be ALL


%% check to criterion matches (tc2, tc2rep, tc3-8,9,17)

SESSION = TOCRITERION2;
fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);
nstim = length(unique(s1_f.stim.id));
fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);

% so for the first 20 trials, should be the same picture and associate

IDXw = find(cellfun(@isequal,s1_f.stim.stim(1:nstim), s2_f.stim.stim(1:nstim))) %this should be ALL
IDXp = find(cellfun(@isequal,s1_f.stim.associate(1:nstim), s2_f.stim.associate(1:nstim))) %this should be ALL

c = isequal(s1_f.stim.choicePos(1:nstim,:),s2_f.stim.choicePos(1:nstim,:))

%% check to see if RSVP is working! can do RSVP and RSVP2 (10, 15 sessions)

SESSION = RSVP2;

fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);
fname = findNewestFile(s1_dir,fullfile(s1_dir, ['EK' num2str(SESSION)  '*.mat']));
s1_ek = load(fname);
s1_trials = table2cell(s1_ek.datastruct.trials);
s1_id = cell2mat(s1_trials(:,8));
s1_dur = cell2mat(s1_trials(:,20));

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['EK' num2str(SESSION)  '*.mat']));
s2_ek = load(fname);
s2_trials = table2cell(s2_ek.datastruct.trials);
s2_id = cell2mat(s2_trials(:,8));
s2_dur = cell2mat(s2_trials(:,20));

IDi = isequal(s1_id,s2_id)
IDd = max(abs(s1_dur - s2_dur))
abs(s1_dur - s2_dur)

%% CHECK IF MOT_PRACTICE 2 IS SHOWING THE RIGHT PRACTICE WORDS & matched lure words, AND HAS THE DOT SPEED FROM THE SUBJECT, motLOC = 18, mot_practice2 = 13
SESSION =MOT{1} %left off: figure out why for mot practice, the conditions aren't exactly the same--did i not match?
fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);
%L1 = s1_f.stim.lureWords;
trialSpeeds = s1_f.stim.speed;

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);
%L2 = s2_f.stim.lureWords;
%lw = find(cellfun(@isequal, L1, L2))
checkCond = sum(sum(s1_f.stim.cond == s2_f.stim.cond))/numel(s1_f.stim.cond)
allw = find(cellfun(@isequal, s1_f.stim.stim, s2_f.stim.stim))

if SESSION>=MOT_LOCALIZER %if localizer or more check that all stimuli are the same
    allw = find(cellfun(@isequal, s1_f.stim.stim, s2_f.stim.stim))
    if SESSION >=20 %if MOT or later
        %also check to see if the dot speeds are exactly the same
        checkSpeed= sum(sum(s1_f.stim.motionSpeed == s2_f.stim.motionSpeed))/numel(s1_f.stim.motionSpeed)
    end
    % nwo check all speeds
    for i = 1:10
        speeds1 = s1_f.rtData.allSpeeds{i};
        speeds2 = s2_f.rtData.allSpeeds{i};
        se(i) = sum(speeds1==speeds2)/length(speeds1);
        t1 = s1_f.rtData.allTimeChanges{i}-s1_f.timing.actualOnsets.motion(1,i);
        t2 = s2_f.rtData.allTimeChanges{i}-s2_f.timing.actualOnsets.motion(1,i);
        t1-t2
        % make this time then closer and retest!
    end
        s1_f.stim.motionSpeed == s2_f.stim.motionSpeed
end
%% CHECK IF RECALL_PRACTICE IS SHOWING WORDS, KEY PRESSES ARE RECORDED and then recall1 and 2: 14, 19 and 23
SESSION = RECALL2;
fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);

fname = findNewestFile(s1_dir,fullfile(s1_dir, ['EK' num2str(SESSION)  '*.mat']));
s1_ek = load(fname);
s1_ek.datastruct.trials; %check that key presses are correct
allw = find(cellfun(@isequal, s1_f.stim.stim, s2_f.stim.stim))

%% associates task is the same --pictures
SESSION = ASSOCIATES;
fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);

allpic = find(cellfun(@isequal, s1_f.stim.stim, s2_f.stim.stim))
checkCond = isempty(find(~isequal(s1_f.stim.cond,s2_f.stim.cond)))
checkid = isempty(find(~isequal(s1_f.stim.id,s2_f.stim.id)))
checkCond2 = isempty(find(~isequal(s1_f.stim.AAP,s2_f.stim.AAP)))
checkCond3 = isempty(find(~isequal(s1_f.stim.AAPID,s2_f.stim.AAPID)))
%% DESCRIPTION TASK
SESSION = DESCRIPTION;
fname = findNewestFile(s1_dir,fullfile(s1_dir, ['mot_realtime05_' num2str(s1) '_' num2str(SESSION)  '*.mat']));
s1_f = load(fname);

fname = findNewestFile(s2_dir,fullfile(s2_dir, ['mot_realtime05_' num2str(s2) '_' num2str(SESSION)  '*.mat']));
s2_f = load(fname);
checkid = isempty(find(~isequal(s1_f.stim.id,s2_f.stim.id)))
checkid = isempty(find(~isequal(s1_f.stim.stim,s2_f.stim.stim)))
