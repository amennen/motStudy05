
% syntax: mot_realtime05(SUBJECT,SESSION,SET_SPEED,scanNum,scanNow)
%
% This is an implementation of a MOT memory experiment designed to
% reactivate memory, adapted from a study by Ken Norman and J. Poppenk. I
% (JP) wrote this before we developed our modern experiment design
% framework and it should not be used as an example of "how to do things"
% -- yet there are times when it is more efficient to adapt awkward,
% functional older code than it is to write nice, tidy new code. It has
% been at least updated to make of the current superpsychtoolbox
% implementation.


%%

function mot_realtime05(SUBJECT,SESSION,SET_SPEED,scanNum,scanNow)

%SUBJECT: Subject number
%SESSION: task that they're going to do (listed below)
%SET_SPEED: if not empty, debug mode is on
%scanNum: (for MOT use only) if you're going to use fmri scans for real
%time
%scanNow: if you're using the scanner right now (0 if not, 1 if yes) aka
%looking for triggers


% note: TR values begin from t=1 (rather than t=0)
ListenChar(2); %suppress keyboard input to window
KbName('UnifyKeyNames');
% initialization declarations
COLORS.MAINFONTCOLOR = [200 200 200];
COLORS.BGCOLOR = [50 50 50];
WRAPCHARS = 70;
%HideCursor;
%% INITIALIZE EXPERIMENT
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
%ANATOMICAL_PREP = ASSOCIATES + 1;
% name strings
SESSIONSTRINGS{SETUP} = 'GENERATE PAIRINGS'; % set up rsvp study learn associates
SESSIONSTRINGS{FAMILIARIZE} = 'FAMILIARIZE1'; % rsvp study learn associates
SESSIONSTRINGS{TOCRITERION1} = 'TOCRITERION1'; % rsvp study learn associates
SESSIONSTRINGS{MOT_PRACTICE} = 'MOT_PRACTICE'; % get used to task during anatomical
SESSIONSTRINGS{MOT_PREP} = 'MOT_PREP'; % realtime rsvp

% day 1
SESSIONSTRINGS{FAMILIARIZE2} = 'FAMILIARIZE2'; % rsvp study learn associates
SESSIONSTRINGS{TOCRITERION2} = 'TOCRITERION2'; % learn paired associates to critereon
SESSIONSTRINGS{TOCRITERION2_REP} = 'TOCRITERION2_REP';
SESSIONSTRINGS{RSVP} = 'FRUIT_HARVEST'; % fruit harvest task practice and familiarization of words
SESSIONSTRINGS{MOT_PRACTICE2} = 'PRACTICE WITH MOT LOCALIZER'; % get used to task during anatomical

% day 2
SESSIONSTRINGS{SCAN_PREP} = 'SCANNER PREPARATION'; % scan prep
SESSIONSTRINGS{MOT_PRACTICE} = 'MOT PRACTICE'; % get used to task during anatomical
SESSIONSTRINGS{MOT_PREP} = 'MOT PREP'; % realtime rsvp
SESSIONSTRINGS{RSVP2} = 'FRUIT HARVEST 2'; % fruit harvest task practice and familiarization of words
SESSIONSTRINGS{FAMILIARIZE3} = 'FAMILIARIZE - MOT LOCALIZER'; % rsvp study learn associates
SESSIONSTRINGS{TOCRITERION3} = 'TRAIN TO CRITERION - MOT LOCALIZER'; % get used to task during anatomical
SESSIONSTRINGS{MOT_LOCALIZER} = 'MOT LOCALIZER'; % get used to task during anatomical
SESSIONSTRINGS{RECALL_PRACTICE} = 'DELIBERATE RECALL PRACTICE'; % get used to task during anatomical
for i = 1:NUM_TASK_RUNS
    SESSIONSTRINGS{MOT{i}} = ['MOT RUN ' num2str(i)];
end
SESSIONSTRINGS{RECALL_PRACTICE} = 'DELIBERATE RECALL PRACTICE'; % get used to task during anatomical
SESSIONSTRINGS{RECALL1} = 'RECALL1'; % baseline recollection, used to train recollection classifier
SESSIONSTRINGS{RECALL2} = 'RECALL2'; % post-test recollection, used to measure effectiveness of manipulation
SESSIONSTRINGS{DESCRIPTION} = 'SCENE DESCRIPTION TASK'; 
SESSIONSTRINGS{ASSOCIATES} = 'ASSOCIATES TASK'; % face-scene classification
% SETUP: prepare experiment
if ~exist('SUBJECT','var'), SUBJECT = -1; end
if ~exist('SESSION','var'),
    SESSION = -1;
end
if ~exist('SET_SPEED','var') || isempty(SET_SPEED) || SET_SPEED <= 0
    SET_SPEED = -1;
    SPEED = 1;
    debug_mode = false;
else
    SPEED = SET_SPEED;
    debug_mode = true;
end


global CATSTRINGS WORD CAR FACE SCENE NORESP
if SESSION > 0 && SESSION < length(SESSIONSTRINGS)
    exp_string_short = ['textlog_sess' int2str(SESSION)];
    exp_string_long = ['EK' int2str(SESSION)];
else exp_string_short = []; exp_string_long = [];
end
[mainWindow WINDOWSIZE COLORS DEVICE TRIGGER WORKING_DIR LOG_NAME MATLAB_SAVE_FILE ivx fmri SLACK] = ...
    initiate_rt(SUBJECT,SESSION,SESSIONSTRINGS,exp_string_short,exp_string_long,SET_SPEED,COLORS);
if isempty(mainWindow), return; end % quit after printing session map
if ~isempty(strfind(TRIGGER,'=')) || ~isempty(strfind(TRIGGER,'5'))
    CURRENTLY_ONLINE = true; %if scanning
else CURRENTLY_ONLINE = false;
end
if ~scanNow
    CURRENTLY_ONLINE = false; %when we're not looking for triggers
end
%CURRENTLY_ONLINE = false; %change this later!!!
CENTER = WINDOWSIZE.pixels/2;
ALLDEVICES = -1;
if ~fmri && DEVICE == -1 %haven't found device yet
    KEYDEVICES = findInputDevice([],'keyboard');
else
    KEYDEVICES = DEVICE;
end

TIMEOUT = 0.050; %always wait 50 ms for trigger
% find base directory
SUBJ_NAME = num2str(SUBJECT);
documents_path = WORKING_DIR;
% if fmri
%     documents_path = ['/Data1/code/motStudy01/'];
% end
data_dir = fullfile(documents_path, 'BehavioralData');
dicom_dir = fullfile('/Data1/code/motStudy05/', 'data', SUBJ_NAME); %where all the dicom information is FOR THAT SUBJECT
if SESSION >= MOT{1}
    runNum = SESSION - MOT{1} + 1;
    classOutputDir = fullfile(dicom_dir,['motRun' num2str(runNum)], 'classOutput/');
end
if ~exist(data_dir,'dir'), mkdir(data_dir); end
ppt_dir = [data_dir filesep SUBJ_NAME filesep];
if ~exist(ppt_dir,'dir'), mkdir(ppt_dir); end
base_path = [fileparts(which('mot_realtime05.m')) filesep];
MATLAB_SAVE_FILE = [ppt_dir MATLAB_SAVE_FILE];
LOG_NAME = [ppt_dir LOG_NAME];


% SETUP: declarations
%this is where stimuli are stored
% basic session info
stim.session = SESSION;
stim.subject = SUBJECT;
stim.sessionName = SESSIONSTRINGS{SESSION};
stim.num_realtime = 10;
stim.num_omit = 10;
stim.num_learn = 8;
stim.num_localizer = 16;
stim.num_total = stim.num_learn + stim.num_realtime + stim.num_omit + stim.num_localizer; 
stim.runSpeed = SPEED;
stim.TRlength = 2*SPEED; %this is the speed
stim.fontSize = 24;
NUM_MOTLURES = 0; %changed here from 5 6/23
NUM_RECALL_LURES = stim.num_realtime + stim.num_omit;
lureCounter = 1000;
% input mapping
if ~fmri
    THUMB = 'z';
    INDEXFINGER='e';
    MIDDLEFINGER='r';
    RINGFINGER='t';
    PINKYFINGER='y';
else
    KEYDEVICES = DEVICE;
    THUMB='1';
    INDEXFINGER='2';
    MIDDLEFINGER='3';
    RINGFINGER='4';
    PINKYFINGER='5';
    %     if strcmp(TRIGGER,'5'); should already be set to =
    %         TRIGGER = '=';
    %     end
end
%[keys, valid_map, cresp, cresp_map] = keyCheck(keys,cresp,strict);

% if test
%     KEYDEVICES = ALLDEVICES;
%     DEVICE = ALLDEVICES;
% end
keys = [THUMB INDEXFINGER MIDDLEFINGER RINGFINGER PINKYFINGER];
keyCell = {THUMB, INDEXFINGER, MIDDLEFINGER, RINGFINGER, PINKYFINGER};
if fmri
    allkeys = [keys TRIGGER];
else
    allkeys = keys;
end
% define key names so don't have to do it later
for i = 1:length(allkeys)
    keys.code(i,:) = getKeys(allkeys(i));
    keys.map(i,:) = zeros(1,256);
    keys.map(i,keys.code(i,:)) = 1;
end

% scales and valid key entries
mc_keys = keyCell(2:5);
mc_keycode = keys.code(2:5,:);
mc_map = sum(keys.map(1:5,:));
subj_keys = keyCell;
subj_keycode = keys.code(1:5,:);
subj_map = sum(keys.map(1:5,:));
target_keys = keyCell([1 5]);
target_keycode = keys.code([1 5],:);
target_map = sum(keys.map([1 5],:));
fruit_key = keyCell(2);
fruit_keycode = keys.code(2,:);
fruit_map = keys.map(2,:);
kbTrig_keycode= keys.code(2,:);
if fmri
    TRIGGER_keycode = keys.code(6,:);
end
recog_map = sum(keys.map([2 3],:));
mc_scale= makeMap({'farL','midL','midR','farR'}, 1:4, mc_keys);
subj_scale = makeMap({'no image', 'generic', '1 detail', '2+ details', 'full'}, 1:5, subj_keys);
target_scale = makeMap({'non-target', 'target'}, [0 1], target_keys);
recog_scale = makeMap({'left','right'}, [1 2], mc_keys(1:2));
rsvp_scale = makeMap({'target'}, 1, fruit_key);

% stimulus categories
WORD = 1;
CAR = 2;
FACE = 3;
SCENE = 4;

% cue data struct
STIMULI = 1;
CLASS = 2;
EXPOSURE_DUR = 3;
ID = 4;
ACT_READOUT = 5;
EXPOS_DELTA = 6;
LAT_MVT = 7;

% mechanical variables
minimumDisplay = 0.25;
ALLMATERIALS = -1;
MAINPROPORTIONS = [0 0 0 1];
VERTICAL = 2;
HORIZONTAL = 1;
CORRECT = 1;
INCORRECT = 0;
NORESP = -1;
WRITE = 1;
PICDIMS = [256 256];
pic_size = [200 200];
%pic_size = [256 256]; %changed from [200 200]
RESCALE_FACTOR = pic_size/PICDIMS;
% multi choice control
FIXED = 1;
VARIABLE = 0;
CHOICES = 4;

% conditions
REALTIME = 1;
OMIT = 2;
LUREWORD = 3;

% OMIT = 3;
% LUREWORD = 4;

FRUIT = 4;
PRACTICE = 5;

LEARN = 6;
LOC = 7;

% stimulus filepaths
MATLAB_STIM_FILE = [ppt_dir 'mot_realtime05_subj_' num2str(SUBJECT) '_stimAssignment.mat'];
%CUELISTFILE_TARGETS = [base_path 'stimuli/text/wordpool_ONLYTARGETS.txt']; %changed for subject 15 just to be sure that the words for practice aren't going to be assigned to more people!
TRAININGCUELIST = [base_path 'stimuli/text/trainingCues.txt'];
MOTCUELIST = [base_path 'stimuli/text/motCues.txt'];
LOCALIZERCUELIST = [base_path 'stimuli/text/localizerCues.txt'];
LOCALIZERLURELIST = [base_path 'stimuli/text/localizerLures.txt'];
CUETARGETFILE = [base_path 'stimuli/text/ed_plants.txt'];

% it reads the text file-but just get the names from going through the
% entire file! (see Megan's script)
% change for all the correct paths here
MOTFOLDER = [base_path 'stimuli/STIM/middlepairs/'];
MOTLURESFOLDER = [base_path 'stimuli/STIM/middlepairs_lures/'];
TRAININGFOLDER = [base_path 'stimuli/STIM/training/'];
TRAININGLURESFOLDER = [base_path 'stimuli/STIM/training_lures/'];
LOCALIZERFOLDER = [base_path 'stimuli/STIM/localizer/'];

% present mapping without keylabels if ppt. is in the scanner
if CURRENTLY_ONLINE && SESSION > TOCRITERION3
    KEY_MAPPING = [base_path 'stimuli/bwvividness.jpg'];
else KEY_MAPPING = [base_path 'stimuli/bwvividness.jpg'];
end

% strings
CATSTRINGS = {'word','car','face','scene'};
CONDSTRINGS = {'realtime','omit','lure','fruit','practice','learn','localizer'};
MOTSTRINGS = {'hard-targ','easy-targ','hard-lure','easy-lure',[],[],'prep'};
MOT_RT_STRINGS = {'rt-targ'};
NOTIFY = 'Great work! You finished the task.\n\nPlease notify your experimenter.';
CONGRATS = 'Great work! You finished the task.\n\nPlease wait for further instructions.';

% modify instructions if ppt. is in the scanner
STILLEXPLAIN = ['Please remember that moving your head even a little during scanning blurs our picture of your brain.'];
STILLREMINDER = ['The scan is now starting.\n\nMoving your head even a little blurs the image, so '...
    'please try to keep your head totally still until the scanning noise stops.\n\n Do it for science!'];
PROGRESS_TEXT = 'INDEX';
final_instruct_continue = ['\n\n-- Press ' PROGRESS_TEXT ' to begin once you understand these instructions --'];
%end

% timing constants
STABILIZATIONTIME = 2 * stim.TRlength;
STILLDURATION = 6 * SPEED;
CONGRATSDURATION = 3*SPEED;
INSTANT = 0.001;

% initialize random seed generator
%s = RandStream('mt19937ar','Seed','shuffle');
seed = randseedwclock();

%RandStream.setGlobalStream(s);
% find a stimulus file if this is not first session
if SESSION ~= SETUP
    try
        load(MATLAB_STIM_FILE);
    catch
        error('No stimulus assignment file detected');
    end
end
Screen('TextSize',mainWindow,stim.fontSize);

% main experiment session switch
switch SESSION
    %% 0. SETUP
    case SETUP
        
        % goal: create all stimuli for subject
        % ordered training pics and training word pairs
        % ordered localizer pics and word pairs
        % ordered mot pics and word pairs
        % extra words for localizer and recall practice?
        % reading in name and just saving as cells
        picList = {};
        pics = {};
        picsCond = {};
        pictureCount = 1;
        for cond = [LEARN LOC REALTIME]
            tempCell = {};
            if cond == LEARN
                fn = TRAININGFOLDER;
                nstim = stim.num_learn;
                dirList = dir(fn);
                dirList = dirList(3:end);
                orderStim = randperm(nstim);
                for s = 1:nstim
                    tempCell{s} = dirList(orderStim(s)).name; % for that category
                    pics{pictureCount} = dirList(orderStim(s)).name;
                    pictureCount = pictureCount + 1;
                end
                picsCond{cond} = tempCell;
                trainPics = tempCell;
            elseif cond == LOC
                fn = LOCALIZERFOLDER;
                nstim = stim.num_localizer;
                dirList = dir(fn);
                dirList = dirList(3:end);
                orderStim = randperm(nstim);
                for s = 1:nstim
                    tempCell{s} = dirList(orderStim(s)).name; % for that category
                    pics{pictureCount} = dirList(orderStim(s)).name;
                    pictureCount = pictureCount + 1;
                end
                picsCond{cond} = tempCell;
            elseif cond == REALTIME
                % get stimuli for both
                fn = MOTFOLDER;
                nstim = stim.num_realtime;
                dirList = dir(fn);
                dirList = dirList(3:end);
                orderStim = randperm(nstim);
                
                allOptions = 1:nstim;
                indoorRT = randperm(nstim,nstim/2);
                indoorOM = allOptions(~ismember(allOptions,indoorRT));
                outdoorRT = randperm(nstim,nstim/2);
                outdoorOM = allOptions(~ismember(allOptions,outdoorRT));
                outdoorRT = outdoorRT + nstim;
                outdoorOM = outdoorOM + nstim;
                RTIND = Shuffle([indoorRT outdoorRT]);
                OMIND = Shuffle([indoorOM outdoorOM]);
                for cond2 = [REALTIME OMIT]
                    if cond2 == REALTIME
                        takefrom = RTIND;
                    else
                        takefrom = OMIND;
                    end
                    for s = 1:nstim
                        tempCell{s} = dirList(takefrom(s)).name;
                        pics{pictureCount} = dirList(takefrom(s)).name;
                        pictureCount = pictureCount + 1;
                    end
                    picsCond{cond2} = tempCell;
                end
            end
        end

        % now choose all the words
        preapredCues = {};
        wordCount = 1;
        for cond=[LEARN LOC REALTIME]
            if cond==LEARN
                trainWords = readStimulusFile(TRAININGCUELIST,stim.num_learn);
                cues{STIMULI}{LEARN} = trainWords;
                for s = 1:stim.num_learn
                    preparedCues{wordCount} = trainWords{s};
                    wordCount = wordCount + 1;
                end
            elseif cond==LOC
                locCueWords = readStimulusFile(LOCALIZERCUELIST,stim.num_localizer);
                cues{STIMULI}{LOC} = locCueWords;
                for s = 1:stim.num_localizer
                    preparedCues{wordCount} = locCueWords{s};
                    wordCount = wordCount + 1;
                end
            elseif cond==REALTIME
                motCueWords = readStimulusFile(MOTCUELIST,stim.num_realtime + stim.num_omit);
                RTWORDS = motCueWords(1:stim.num_realtime);
                OMITWORDS = motCueWords(stim.num_omit+1:end);
                for cond2 = [REALTIME OMIT]
                    if cond2 == REALTIME
                        takefrom = RTWORDS;
                        cues{STIMULI}{REALTIME} = RTWORDS;
                    else
                        takefrom = OMITWORDS;
                        cues{STIMULI}{OMIT} = OMITWORDS;
                    end
                    for s = 1:stim.num_realtime
                        preparedCues{wordCount} = takefrom{s};
                        wordCount = wordCount + 1;
                    end
                end
            end
        end
        
        pairIndex = 1:stim.num_total;
        
        NUMLURES = 23; %hardcoded
        lureWords = readStimulusFile(LOCALIZERLURELIST,NUMLURES);
        stimmap = makeMap(preparedCues);
        
        % clean up
        system(['cp ' base_path 'mot_realtime05.m ' ppt_dir 'mot_realtime05_executed.m']);
        save(MATLAB_STIM_FILE, 'cues','preparedCues','pics','pairIndex','lureWords','stimmap', 'trainWords', 'trainPics');
        
        
        mot_realtime05(SUBJECT,SESSION+1,SET_SPEED,scanNum,scanNow);
        
        
        %% 1. FAMILIARIZATION
    case {FAMILIARIZE, FAMILIARIZE2, FAMILIARIZE3,STIM_REFRESH}
        % print headers to file and command window - don't need to wait for
        % trigger pulses here
        header = ['\n\n\n*******************' exp_string_long '*******************\n'];
        header2 = ['SESSION ' int2str(SESSION) ' initiated ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n'];
        printlog(LOG_NAME,['\nSESSION ' int2str(SESSION) ' initiated ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        
        fprintf(header);
        fprintf(header2);
        % declarations
        stim.cueDuration = 2*SPEED;    % cue word alone for 0ms
        stim.picDuration = 4*SPEED;    % cue with associate for 4000ms
        stim.isiDuration = 2*SPEED;
        stim.textRow = WINDOWSIZE.pixels(2)* (1/5); %changed from 5, then 3
        stim.picRow = WINDOWSIZE.pixels(2) *5/9;
        NUMRUNS = 1;
        PROGRESS = INDEXFINGER;
        PROGRESS_TEXT = 'INDEX';
        
        PF = MOTFOLDER;
        % initialization
        trial = 0;
        printlog(LOG_NAME,'session\ttrial\tpair\tonset\tdur\tcue         \tassociate   \n');

        if SESSION == FAMILIARIZE
            [stim.cond stim.condString stimList] = counterbalance_items({cues{STIMULI}{LEARN}},{CONDSTRINGS{LEARN}}); %this just gets the cue words
            picList = lutSort(stimList, preparedCues, pics);
            IDlist = lutSort(stimList, preparedCues, pairIndex);
            PF = TRAININGFOLDER;
        elseif SESSION == FAMILIARIZE2 || SESSION == STIM_REFRESH
            [stim.cond stim.condString stimList] = counterbalance_items({cues{STIMULI}{REALTIME}, cues{STIMULI}{OMIT}},CONDSTRINGS);
            picList = lutSort(stimList, preparedCues, pics);
            IDlist = lutSort(stimList, preparedCues, pairIndex);
        elseif SESSION == FAMILIARIZE3
            [stim.cond stim.condString stimList] = counterbalance_items({cues{STIMULI}{LOC}},{CONDSTRINGS{LOC}});
            picList = lutSort(stimList, preparedCues, pics);
            IDlist = lutSort(stimList, preparedCues, pairIndex);
            PF = LOCALIZERFOLDER;
        end
        
        
        if SESSION == FAMILIARIZE2
            halfway = ['Great job, youre''re halfway there! You can take a stretching or bathroom break if you need to now. \n\n-- press ' PROGRESS_TEXT ' to continue when you''re ready. --'];
            displayText(mainWindow,halfway,INSTANT,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,KEYDEVICES);
        end

        
        instruct = ['NAMED SCENES\n\nToday, you will learn the names of ' num2str(length(stimList)) ' different scenes. ' ...
            'It is important that you learn these now, as you will need to be able to picture each scene based on its name throughout our study.\n\n' ...
            'In this section, you will get a chance to study each name-scene pair, one pair at a time. To help you remember each pair, try to imagine how ' ...
            'each scene got its name - the more vivid and unique, the better your memory will be. However, you will see each scene for only four seconds, ' ...
            'so you will need to work quickly.\n\n-- press ' PROGRESS_TEXT ' to begin --'];
        instruct2 = ['Now we will repeat the NAMED SCENE task from before, but this time we will be using different scenes'];
        if SESSION == STIM_REFRESH
             instruct = ['NAMED SCENES\n\nToday, we''ll start with a quick refresher of the ' num2str(length(stimList)) ' different scenes you learned yesterday. ' ...
            'It is important that you focus, as you will need to be able to picture each scene based on its name throughout our study.\n\n' ...
            'In this section, you will get a chance to study each name-scene pair, one pair at a time. To help you remember each pair, try to imagine how ' ...
            'each scene got its name - the more vivid and unique, the better your memory will be. However, you will see each scene for only four seconds, ' ...
            'so you will need to work quickly.\n\n-- press ' PROGRESS_TEXT ' to begin --'];
        end
        displayText(mainWindow,instruct,INSTANT,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        
        
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,KEYDEVICES);
        runStart = GetSecs;
        %first isi here
        config.TR = stim.TRlength;
        config.nTRs.ISI = stim.isiDuration/stim.TRlength;
        config.nTRs.cue = stim.cueDuration/stim.TRlength;
        config.nTRs.pic = stim.picDuration/stim.TRlength;
        config.nTrials = length(stimList);
        
        config.nTRs.perTrial =  config.nTRs.ISI + config.nTRs.cue + config.nTRs.pic;
        config.nTRs.perBlock = (config.nTRs.perTrial)*config.nTrials+ config.nTRs.ISI; %includes the last ISI and 20 s fixation in the beginning
        % calculate all future onsets
        timing.plannedOnsets.preITI(1:config.nTrials) = runStart + ((0:config.nTrials-1)*config.nTRs.perTrial)*config.TR;
        timing.plannedOnsets.cue(1:config.nTrials) = timing.plannedOnsets.preITI + config.nTRs.ISI*config.TR;
        timing.plannedOnsets.pic(1:config.nTrials) = timing.plannedOnsets.cue + config.nTRs.cue*config.TR;
        timing.plannedOnsets.lastITI = timing.plannedOnsets.pic(end) + config.nTRs.pic*config.TR;
        
        
        % repeat trials for full stimulus set
        for n=1:length(stimList)
            
            %present ISI
            timespec = timing.plannedOnsets.preITI(n) - SLACK;
            timing.actualOnsets.preITI(n) = isi_specific(mainWindow,COLORS.MAINFONTCOLOR, timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.preITI(n) - timing.plannedOnsets.preITI(n));
            
            trial = trial + 1;
            stim.stim{trial} = stimList{n};
            stim.picStim{trial} = picList{n};
            stim.id(trial) = IDlist(n);
            %stim.onsetTime(trial) = GetSecs - stim.expStartTime;
            
            %present cue
            timespec = timing.plannedOnsets.cue(n) - SLACK; %subtract so it's coming at the next possible refresh
            timing.actualOnsets.cue(n) = displayText_specific(mainWindow,stim.stim{trial},stim.textRow,COLORS.MAINFONTCOLOR,WRAPCHARS, timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.cue(n) - timing.plannedOnsets.cue(n));
            % prepare cue+associate window
            %RESCALE_FACTOR = (PICDIMS/WINDOWSIZE.pixels)/4;%be 1/4 of the screen
            picIndex = prepImage(strcat(PF, stim.picStim{trial}),mainWindow);
            topLeft(HORIZONTAL) = CENTER(HORIZONTAL) - (PICDIMS(HORIZONTAL)*RESCALE_FACTOR)/2;
            topLeft(VERTICAL) = stim.picRow - (PICDIMS(VERTICAL)*RESCALE_FACTOR)/2;
            %putting word and image
            
            DrawFormattedText(mainWindow,stim.stim{trial},'center',stim.textRow,COLORS.MAINFONTCOLOR,WRAPCHARS);
            Screen('DrawTexture', mainWindow, picIndex, [0 0 PICDIMS],[topLeft topLeft+PICDIMS.*RESCALE_FACTOR]);
            % show cue with associate
            timespec = timing.plannedOnsets.pic(n) - SLACK;
            timing.actualOnsets.pic(n) = Screen('Flip',mainWindow,timespec);% pass third input to screen flip tell it when you want it to happen
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.pic(n) - timing.plannedOnsets.pic(n));
            
            %Screen('Close',picIndex);
            offsetTime = GetSecs;
            % report trial data
            printlog(LOG_NAME,'%d\t%d\t%d\t%-6s\t%-6s\t%-12s\t%-12s\n',SESSION,trial,IDlist(n),int2str((timing.actualOnsets.cue(n)-stim.subjStartTime)*1000),int2str((offsetTime-timing.actualOnsets.cue(n))*1000),stimList{n},picList{n});
        end
        
        timespec = timing.plannedOnsets.lastITI - SLACK;
        timing.actualOnsets.lastITI = isi_specific(mainWindow,COLORS.MAINFONTCOLOR, timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.lastITI- timing.plannedOnsets.lastITI);
        WaitSecs(stim.isiDuration);
        % clean up
        save(MATLAB_SAVE_FILE,'stim','config', 'timing');
        
        displayText(mainWindow,CONGRATS,CONGRATSDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        printlog(LOG_NAME,['\n\nSESSION ' int2str(SESSION) ' ended ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        printlog(LOG_NAME,'\n\n\n******************************************************************************\n');
        % return
        sca
        if SESSION ~= STIM_REFRESH
            mot_realtime05(SUBJECT,SESSION+1,SET_SPEED,scanNum,scanNow); %don't continue if just refreshing on laptop
        end
        %% 2. LEARN TO CRITERION
    case {TOCRITERION1, TOCRITERION2, TOCRITERION2_REP, TOCRITERION3}
        
        % stimulus presentation parameters
        stim.cueDuration = 4*SPEED;
        stim.feedbackDuration = 0.5*SPEED;
        stim.reStudyDuration = 4*SPEED;
        stim.isiDuration = 4*SPEED;
        stim.choiceWidth = WINDOWSIZE.pixels(HORIZONTAL) / (CHOICES+1);
        stim.gapWidth = 10;
        stim.goodFeedback = '!!!';
        stim.badFeedback = 'X';
        stim.textRow = WINDOWSIZE.pixels(2) *(2.5/9); %changed from 5, then 3
        stim.picRow = WINDOWSIZE.pixels(2) *(5/9);
        
        % other constants
        n=0;
        PROGRESS = INDEXFINGER;
        PROGRESS_TEXT = 'INDEX';
        
        % stimulus data fields
        stim.trial = 0;
        stim.triggerCounter = 1;
        stim.missedTriggers = 0;
        stim.loopNumber = 1;
        % sequence preparation
        if SESSION == TOCRITERION1
            [cond strings stimList] = counterbalance_items({cues{STIMULI}{LEARN}},{CONDSTRINGS{LEARN}});
            condmap = makeMap({'stair'});
            pics = lutSort(stimList, preparedCues, pics);
            pairIndex = lutSort(stimList, preparedCues, pairIndex);
            PF = TRAININGFOLDER;
        elseif SESSION == TOCRITERION2 || SESSION == TOCRITERION2_REP
            [cond strings stimList] = counterbalance_items({cues{STIMULI}{REALTIME}, cues{STIMULI}{OMIT}},CONDSTRINGS);
            condmap = makeMap({'realtime','omit'});
            pics = lutSort(stimList, preparedCues, pics);
            pairIndex = lutSort(stimList, preparedCues, pairIndex);
            PF = MOTFOLDER;
        elseif SESSION == TOCRITERION3
            [cond strings stimList] = counterbalance_items({cues{STIMULI}{LOC}},{CONDSTRINGS{LOC}});
            condmap = makeMap({'localizer'});
            pics = lutSort(stimList, preparedCues, pics);
            pairIndex = lutSort(stimList, preparedCues, pairIndex);
            PF = LOCALIZERFOLDER;
        end
        %first pics is all pics, preparedCues is all cues--pics is then the
        %stimuli that were used in the run
        
        preparedCues = stimList;
        stim.gotItem(1:length(preparedCues)) = NORESP;
        
        % initialize questions
        mc_promptDur = 3.5 * SPEED;
        mc_listenDur = 0 * SPEED;
        subj_triggerNext = false;
        mc_triggerNext = false;
        subj_promptDur = 4 * SPEED;
        subj_listenDur = 0 * SPEED;
        
        
        % instructions
        instruct1 = ['SCENE TEST\n\nNow we will test your memory. First, we will show a scene name. In your mind''s eye, try to picture the ' ...
            'asssociated scene in as much detail as possible, including any specific objects or features. Try to simulate actually looking ' ...
            'at the picture. Because we are studying mental imagery, it is crucial for our experiment that you really do this, both here and ' ...
            'later on'...
            '\n\n--- Press ' PROGRESS_TEXT ' once you understand these instructions ---'];
        
        instruct2 = ['After you respond, we will show you four scenes. Using your index, middle, ring and pinky fingers ' ...
            '(corresponding to the leftmost to rightmost image), please identify the scene that goes with the current name.\n\nIf you ' ...
            'suspect one scene, select it even if you''re unsure; but if you really have no idea at all which picture is correct, do not guess randomly. ' ...
            'If you have no idea, use your THUMB to SKIP the question. Thus, you should be responding every trial, using either your INDEX ' ...
            'to PINKY fingers to choose an image, or using your THUMB if you don''t know which image to choose. '...
            '\n\nDepending how you answer, you will see a green "!!!" (correct) or a red "X" (incorrect). The task will continue until you have ' ...
            'named all scenes correctly. Good luck!\n\n--- Press ' PROGRESS_TEXT ' to begin ---'];
        
        % present instructions
        if SESSION~= TOCRITERION2_REP
            DrawFormattedText(mainWindow,' ','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            displayText(mainWindow,instruct1,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        else
            instruct2 = ['We will now repeat the SCENE TEST. You will have to name each scene correctly again. '...
                'Good luck!\n\n--- Press ' PROGRESS_TEXT ' to begin ---'];
            DrawFormattedText(mainWindow,' ','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            displayText(mainWindow,instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        end
        
        
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,KEYDEVICES);
       
        objectiveEK = initEasyKeys([exp_string_long '_OB'], SUBJ_NAME,ppt_dir, ...
            'default_respmap', subj_scale, ...
            'stimmap', stimmap, ...
            'condmap', condmap, ...
            'trigger_next', mc_triggerNext, ...
            'prompt_dur', mc_promptDur, ...
            'listen_dur', mc_listenDur, ...
            'exp_onset', stim.subjStartTime, ...
            'console', false, ...
            'device', DEVICE);
        [objectiveEK] = startSession(objectiveEK);
        
        runStart = GetSecs;
        config.TR = stim.TRlength;
        config.nTRs.ISI = stim.isiDuration/stim.TRlength;
        config.nTRs.cue = stim.cueDuration/stim.TRlength;
        %config.nTRs.vis = subj_promptDur/stim.TRlength;
        config.nTRs.mc = (mc_promptDur + stim.feedbackDuration)/stim.TRlength;
        config.nTRs.reStudy = stim.reStudyDuration/stim.TRlength;
        
        config.nTRs.trial(2) = config.nTRs.ISI + config.nTRs.cue + config.nTRs.mc;
        config.nTRs.trial(1) = config.nTRs.trial(2) + config.nTRs.reStudy;
        
        while n < length(preparedCues)
            % initialize trial
            n = n+1;
            if stim.gotItem(n) ~= CORRECT
                stim.gotItem(n) = NORESP;
                stim.trial = stim.trial + 1; %this is different from n!!! n is for stimulus number, stim.trial is trial number
                
                timing.plannedOnsets.preITI(stim.trial) = runStart;
                if stim.trial > 1
                    timing.plannedOnsets.preITI(stim.trial) = timing.plannedOnsets.preITI(stim.trial - 1) + config.nTRs.trial(lastacc)*config.TR;
                end
                timespec = timing.plannedOnsets.preITI(stim.trial) - SLACK;
                timing.actualOnsets.preITI(stim.trial) = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
                fprintf('Flip time error = %.4f\n', timing.actualOnsets.preITI(stim.trial) - timing.plannedOnsets.preITI(stim.trial));
                
                % show cue window
                stim.stim{stim.trial} = preparedCues{n};
                stim.associate{stim.trial} = pics{n};
                stim.id(stim.trial) = pairIndex(n);
                stim.cond(stim.trial) = cond(n);
                stim.condString{stim.trial} = strings{n};
                icresp = 3:5;
                cresp = keyCell(icresp);
                cresp_map = sum(keys.map(icresp,:));
                
                timing.plannedOnsets.cue(stim.trial) = timing.plannedOnsets.preITI(stim.trial) + config.nTRs.ISI*config.TR;
                timespec = timing.plannedOnsets.cue(stim.trial) - SLACK;
                timing.actualOnsets.cue(stim.trial) = displayText_specific(mainWindow,preparedCues{n},'center',COLORS.MAINFONTCOLOR,WRAPCHARS,timespec);
                fprintf('Flip time error = %.4f\n', timing.actualOnsets.cue(stim.trial) - timing.plannedOnsets.cue(stim.trial));
                
                % choose lures
                cpos{stim.trial} = randi(4);
                lureIndex = setdiff(1:CHOICES,cpos{stim.trial}); %this is the index in the positioning of the lure images
                temp_pics = pics; %these are the 5 used in the learning trial
                numLures = 0;
                fullidx = 1:length(temp_pics);
                notidx = n;
                good_idx = setdiff(fullidx,notidx); %this is the index for all images in the data set
                inside = isempty(strfind(pics{n}, 'o')); %if this is true, then the trial's image is an indoor image
                outside = ~inside;
                allOtherPics = pics;
                allOtherPics(n) = [];
                indexC = strfind(allOtherPics, 'o');
                allOtherOutside = find(not(cellfun('isempty', indexC)));
                allOtherInside = setdiff(1:length(allOtherPics),allOtherOutside);
                
                % choose 2 of each inside/outside
                if inside %choose 2 outside
                    outsidePics = allOtherOutside(randperm(length(allOtherOutside),2));
                    insidePics = allOtherInside(randperm(length(allOtherInside),1));
                else % if outside, choose 2 inside
                    outsidePics = allOtherOutside(randperm(length(allOtherOutside),1));
                    insidePics = allOtherInside(randperm(length(allOtherInside),2));
                end
                lureIndices = Shuffle([outsidePics insidePics]); %shuffles the order of 2 inside, 2 outside but they're indexes
                
               
                
                for i=1:length(lureIndices)
                    stim.choicePos(stim.trial,lureIndex(i)) = lureIndices(i);
                    picLures{i} = allOtherPics{lureIndices(i)};
                    picIndex(lureIndex(i)) = prepImage(char(strcat(PF, picLures{i})),mainWindow);
                end
                
                % close and replace the lure in the spot that belongs to the target
                picIndex(cpos{stim.trial}) = prepImage(strcat(PF, pics{n}),mainWindow); %%ooooh here you take the lure out and add target
                
                stim.choicePos(stim.trial,cpos{stim.trial}) = n;
                
                % draw exemplar options
                destDims = [200 200];
                %destDims = min(PICDIMS*RESCALE_FACTOR,PICDIMS .* (stim.choiceWidth ./ PICDIMS(HORIZONTAL)));
                RESCALE_NEW = destDims(2)/PICDIMS(2);
                topLeft(HORIZONTAL) = CENTER(HORIZONTAL) - (destDims(HORIZONTAL)*CHOICES/2) - (stim.gapWidth*CHOICES/2);
                topLeft(VERTICAL) = stim.picRow - (PICDIMS(VERTICAL)*RESCALE_NEW)/2;
                Screen('FillRect', mainWindow, COLORS.BGCOLOR);
                for i=1:CHOICES
                    Screen('DrawTexture', mainWindow, picIndex(i), [0 0 PICDIMS],[topLeft topLeft+destDims]);
                    topLeft(HORIZONTAL) = topLeft(HORIZONTAL) + destDims(HORIZONTAL) + stim.gapWidth;
                end
                DrawFormattedText(mainWindow,preparedCues{n},'center',stim.textRow,COLORS.MAINFONTCOLOR,WRAPCHARS);
                %multiple choice
                timing.plannedOnsets.mc(stim.trial) = timing.plannedOnsets.cue(stim.trial) + config.nTRs.cue*config.TR;
                timespec =  timing.plannedOnsets.mc(stim.trial) - SLACK;
                %display multiple choice options here
                timing.actualOnsets.mc(stim.trial) = Screen('Flip',mainWindow, timespec);
                fprintf('Flip time error = %.4f\n', timing.actualOnsets.mc(stim.trial) - timing.plannedOnsets.mc(stim.trial));
                objectiveEK = easyKeys(objectiveEK, ...
                    'onset', timing.actualOnsets.mc(stim.trial), ...
                    'stim', stim.stim{stim.trial}, ...
                    'cond', stim.cond(stim.trial), ...
                    'nesting', [SESSION stim.loopNumber stim.trial], ...
                    'cresp', keyCell(cpos{stim.trial}+1), ...
                    'cresp_map', keys.map(cpos{stim.trial}+1,:), 'valid_map', mc_map) ;
                
                % clean up
                for i=1:CHOICES
                    Screen('Close',picIndex(i));
                end
                
                % score this trial
                timing.plannedOnsets.fb(stim.trial) = timing.plannedOnsets.mc(stim.trial) + mc_promptDur;
                timespec = timing.plannedOnsets.fb(stim.trial) - SLACK;
                if objectiveEK.trials.acc(end) == CORRECT
                    stim.gotItem(n) = CORRECT;
                    timing.actualOnsets.fb(stim.trial) = displayText_specific(mainWindow,stim.goodFeedback,'center',COLORS.GREEN,WRAPCHARS, timespec);
                    lastacc = 2;
                else
                    stim.gotItem(n) = INCORRECT;
                    timing.actualOnsets.fb(stim.trial) = displayText_specific(mainWindow,stim.badFeedback,'center',COLORS.RED,WRAPCHARS, timespec);
                    lastacc = 1;
                end
                
                fprintf('Flip time error = %.4f\n', timing.actualOnsets.fb(stim.trial) - timing.plannedOnsets.fb(stim.trial));
                
                % present cue+associate window
                if stim.gotItem(n) ~= CORRECT
                    picIndex(1) = prepImage(strcat(PF, pics{n}),mainWindow);
                    topLeft(HORIZONTAL) = CENTER(HORIZONTAL) - (PICDIMS(HORIZONTAL)*RESCALE_FACTOR/2);
                    topLeft(VERTICAL) = stim.picRow - (PICDIMS(VERTICAL)*RESCALE_FACTOR)/2;
                    DrawFormattedText(mainWindow,preparedCues{n},'center',stim.textRow,COLORS.MAINFONTCOLOR,WRAPCHARS);
                    Screen('DrawTexture', mainWindow, picIndex(1), [0 0 PICDIMS],[topLeft topLeft+PICDIMS*RESCALE_FACTOR]);
                    
                    timing.plannedOnsets.reStudy(stim.trial) = timing.plannedOnsets.mc(stim.trial) + config.nTRs.mc*config.TR;
                    timespec = timing.plannedOnsets.reStudy(stim.trial) - SLACK;
                    timing.actualOnsets.reStudy(stim.trial) = Screen('Flip',mainWindow,timespec);
                    fprintf('Flip time error = %.4f\n', timing.actualOnsets.reStudy(stim.trial) - timing.plannedOnsets.reStudy(stim.trial));
                    
                    % Screen('Close',picIndex(1));
                end
                
                
                save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
            end
            
            
            % present feedback
            if SESSION == TOCRITERION1 && stim.trial <= 3
                
                if lastacc == 2 % if correct
                    WaitSecs(stim.feedbackDuration);
                else
                    WaitSecs(stim.reStudyDuration);
                end
                onFB = GetSecs;
                
                if isnan(objectiveEK.trials.resp(end))
                    displayText(mainWindow,['In the multiple choice section, it looks like you either did not provide a response in ' ...
                        'time, or did not use an appropriate key.\n\nRemember that when you have no idea, you should use your THUMB as a SKIP.\n\nAlso remember that your index, middle, ring or pinky ' ...
                        'finger correspond to left, left middle, right middle and right.\n\n' ...
                        '-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                end
                if ~isnan(objectiveEK.trials.resp(end))
                    displayText(mainWindow,['Good work! Your multiple choice response were both detected.\n\n' ...
                        '-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                end
                offFB = GetSecs;
                timing.plannedOnsets.preITI(stim.trial) = timing.plannedOnsets.preITI(stim.trial) + (offFB-onFB); %only need to change first one bc it's within a trial that it updates (ohh change this after!)
            end
            save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
            
            
            % when training to criteron, shuffle and repeat incorrect items on list until all items correct
            if (n == length(preparedCues)) && (sum(stim.gotItem) < length(preparedCues)) % do this only on the last time
                n = 0;
                stim.loopNumber = stim.loopNumber + 1;
                revisedOrder = randperm(length(preparedCues));
                preparedCues = preparedCues(revisedOrder);
                pics = pics(revisedOrder);
                cond = cond(revisedOrder);
                strings = strings(revisedOrder);
                pairIndex = pairIndex(revisedOrder);
                stim.gotItem = stim.gotItem(revisedOrder);
            end
            
        end
        
        %must have gotten the last choice right to exit
        timing.plannedOnsets.lastITI = timing.plannedOnsets.mc(stim.trial) + config.nTRs.mc*config.TR;
        timespec = timing.plannedOnsets.lastITI - SLACK;
        timing.actualOnsets.lastITI = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.lastITI - timing.plannedOnsets.lastITI);
        WaitSecs(stim.isiDuration);
        
        % preserve IV's, DV's for later analysis by writing stuff to file here
        save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        
        % clean up
        printlog(LOG_NAME,['\n\nSESSION ' int2str(SESSION) ' ended ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        printlog(LOG_NAME,'\n\n\n******************************************************************************\n');
        %endSession(subjectiveEK, objectiveEK, CONGRATS);
        endSession(objectiveEK, CONGRATS);
        sca
        %return
        %normally would go to session 4 but instead we want to go to
        if SESSION ~= TOCRITERION3
            mot_realtime05(SUBJECT,SESSION+1,SET_SPEED,scanNum,scanNow);
        end
        
        %% 3. PRE/POST MEMORY TEST
    case {RECALL_PRACTICE,RECALL1,RECALL2}
        % stimulus presentation parameters
        stim.promptDur = 8*SPEED; % cue word alone
        stim.prepDur = 2*SPEED;
        stim.isiDuration = 4*SPEED; %
        stim.fixBlock = 20*SPEED;
        num_digit_qs = 3;
        digits_promptDur = 1.9*SPEED;
        digits_isi = 0.1*SPEED;
        digits_triggerNext = false;
        minimal_format = true;
        subj_triggerNext = false;
        keymap_image = imread(KEY_MAPPING);
        subj_promptDur = 4 * SPEED;
        subj_listenDur = 0 * SPEED;
        % stimulus data fields
        stim.triggerCounter = 1;
        stim.missedTriggers = 0;
        PROGRESS_TEXT = 'INDEX';
        
        
        % all the instructions
        stim.instruct1 = ['MENTAL PICTURES TASK\n\nThis is a memory test, though it is not multiple choice anymore. When a word appears, picture the scene it names as vividly as if you were looking ' ...
            'at it now. Picture its objects, features, and anything else you can imagine. HOLD the image, letting it continue to mature ' ...
            'and take shape for the entire ' num2str(stim.promptDur) ' seconds it is on the screen. Don''t let it fade, and don''t let ' ...
            'yourself "space out". It is essential you actually do this!\n\n' ...
            '-- Press ' PROGRESS_TEXT ' once you understand these instructions --'];
        stim.instruct2 = ['After you have had several seconds to form an image, we will ask you how detailed it was. This rating is not a "test": ' ...
            'we want to understand your real experience, even if no image formed when you felt you SHOULD have had one. Similarly, you should not ' ...
            'adjust your rating based on how sure you are it''s the correct scene. Please make your rating based only on how much scene detail comes to you, ' ...
            'ignoring other mental imagery (e.g., you pictured a beach). You will use the provided rating scale to indicate details, using all five fingers on the key pad.'...
            '\n\n---- Press ' PROGRESS_TEXT ' once you understand these instructions,\n then press it again when you are done viewing the rating scale ---'];
        
   
        % initialize stimulus order with initial warmup item
        if SESSION == RECALL_PRACTICE
            stim.stim = cues{STIMULI}{LEARN}(1:3);
            stim.cond = [PRACTICE, PRACTICE, PRACTICE];
            condmap = makeMap({'PRACTICE'});
            stim.condString = {CONDSTRINGS{PRACTICE}, CONDSTRINGS{PRACTICE}, CONDSTRINGS{PRACTICE}};
            displayText(mainWindow,['MENTAL PICTURES TASK - PRACTICE\n\nThis is a memory test with some differences ' ...
                'from earlier: you will get only one try per item and there is no feedback. Instead, you will be alizing the scenes.' ...
                'Please carefully review today''s instructions, since many things have ' ...
                'changed and it is important you follow them exactly.\n\n' ...
                '-- Please press ' PROGRESS_TEXT ' to briefly review the instructions --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
        else
            [stim.cond stim.condString stim.stim] = counterbalance_items({cues{STIMULI}{REALTIME}, cues{STIMULI}{OMIT}},CONDSTRINGS);
            condmap = makeMap({'realtime','omit'});
            if SESSION ==RECALL1
                displayText(mainWindow,['The experiment will now ONLY involve the stimuli that you studied yesterday, both ' ...
                    'for this next task and the rest of the experiment.\n\n' ...
                    '-- Please press ' PROGRESS_TEXT ' to read the instructions for your next task --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                waitForKeyboard(kbTrig_keycode,DEVICE);
            end
            if SESSION == RECALL2
                displayText(mainWindow,['MENTAL PICTURES TASK\n\nThe format of this task is the same as earlier in today''s session.\n\n' ...
                    '-- Please press ' PROGRESS_TEXT ' to briefly review the instructions again --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                waitForKeyboard(kbTrig_keycode,DEVICE);
            end
        end
        
        %generate stimulus ID's first so can add them easily
        for i = 1:length(stim.cond)
            if SESSION == RECALL_PRACTICE
                pos = find(strcmp(cues{STIMULI}{LEARN},stim.stim{i}));
                stim.id(i) = pos;
            else
                pos = find(strcmp(preparedCues,stim.stim{i}));
                stim.id(i) = pos;
            end
        end
        
        % display instructions
        DrawFormattedText(mainWindow,' ','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        displayText(mainWindow,stim.instruct1,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        waitForKeyboard(kbTrig_keycode,DEVICE);
        displayText(mainWindow,stim.instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        waitForKeyboard(kbTrig_keycode,DEVICE);
        keymap_image = imread(KEY_MAPPING);
        keymap_prompt = Screen('MakeTexture', mainWindow, keymap_image);
        Screen('DrawTexture',mainWindow,keymap_prompt,[],[],[]); %[0 0 keymap_dims],[topLeft topLeft+keymap_dims]);
        Screen('Flip',mainWindow);
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);
        %displayText(mainWindow,stim.instruct3,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        %stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);
        
        %last instructions
        if SESSION == RECALL_PRACTICE
            displayText(mainWindow,['To help you get used to the feel of this task, we will ' ...
                'now give you three practice words.\n\n' ...
                '-- Press ' PROGRESS_TEXT ' to begin once you understand these instructions --'], ...
                minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);
        end
        
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            DrawFormattedText(mainWindow,'Waiting for scanner start, hold tight!','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            Screen('Flip', mainWindow);
        end
        
        subjectiveEK = initEasyKeys([exp_string_long '_SUB'], SUBJ_NAME,ppt_dir, ...
            'default_respmap', subj_scale, ...
            'stimmap', stimmap, ...
            'condmap', condmap, ...
            'trigger_next', subj_triggerNext, ...
            'prompt_dur', subj_promptDur, ...
            'listen_dur', subj_listenDur, ...
            'exp_onset', stim.subjStartTime, ...
            'console', false, ...
            'device', DEVICE);
        
        digits_scale = makeMap({'even','odd'},[0 1],keyCell([1 5]));
        condmap = makeMap({'even','odd'});
        digitsEK = initEasyKeys('odd_even', SUBJ_NAME, ppt_dir,...
            'default_respmap', digits_scale, ...
            'condmap', condmap, ...
            'trigger_next', digits_triggerNext, ...
            'prompt_dur', digits_promptDur, ...
            'device', DEVICE);
        
        [subjectiveEK] = startSession(subjectiveEK);
        
        % fixation period for 20 s
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            [timing.trig.wait timing.trig.waitSuccess] = WaitTRPulse(TRIGGER_keycode,DEVICE);
            runStart = timing.trig.wait;
            displayText(mainWindow,STILLREMINDER,STILLDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            DrawFormattedText(mainWindow,'+','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            Screen('Flip', mainWindow)
            config.wait = stim.fixBlock; % stim.fixation - 20 s? so 10 TRs
        else
            runStart = GetSecs;
            config.wait = 0;
        end
        
        config.TR = stim.TRlength;
        config.nTRs.ISI = stim.isiDuration/stim.TRlength;
        config.nTRs.prompt = stim.promptDur/stim.TRlength;
        config.nTRs.vis = subj_promptDur/stim.TRlength;
        config.nTRs.math = (num_digit_qs*(digits_promptDur + digits_isi))/stim.TRlength;
        config.nTrials = length(stim.stim);
        config.nTRs.perTrial = (config.nTRs.ISI + config.nTRs.prompt + config.nTRs.vis + config.nTRs.math);
        config.nTRs.perBlock = config.wait/config.TR + (config.nTRs.perTrial)*config.nTrials+ config.nTRs.ISI; %includes the last ISI
        
        % calculate all future onsets
        timing.plannedOnsets.preITI(1:config.nTrials) = runStart + config.wait + ((0:config.nTrials-1)*config.nTRs.perTrial)*config.TR;
        timing.plannedOnsets.prompt(1:config.nTrials) = timing.plannedOnsets.preITI + config.nTRs.ISI*config.TR;
        timing.plannedOnsets.vis(1:config.nTrials) = timing.plannedOnsets.prompt + config.nTRs.prompt*config.TR;
        timing.plannedOnsets.math(1:config.nTrials) = timing.plannedOnsets.vis + config.nTRs.vis*config.TR;
        timing.plannedOnsets.lastITI = timing.plannedOnsets.math(end) + config.nTRs.math*config.TR;%%make sure it pauses for this one
        
        cresp = keyCell(3:5);
        cresp_map = sum(keys.map(3:5,:));
        
        stimID = stim.id;
        stimCond = stim.cond;
        sessionInfoFile = fullfile(ppt_dir, ['SessionInfo' '_' num2str(SESSION)]);
        save(sessionInfoFile, 'stimCond','stimID', 'timing', 'config'); 
        
        for n = 1:length(stim.stim)
            % initialize trial and show cue
            stim.trial = n;
            fprintf(['Trial number: ' num2str(n) '\n']);
            
            %show pre ITI
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.preITI(n), timing.trig.preITI_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE,timing.plannedOnsets.preITI(n));
            end
            timespec = timing.plannedOnsets.preITI(n) - SLACK;
            timing.actualOnsets.preITI(n) = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.preITI(n) - timing.plannedOnsets.preITI(n));
            
            %display word
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.prompt(n), timing.trig.prompt_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE,timing.plannedOnsets.prompt(n));
            end
            timespec = timing.plannedOnsets.prompt(n)-SLACK;
            timing.actualOnsets.prompt(n) = displayText_specific(mainWindow,stim.stim{stim.trial},'center',COLORS.MAINFONTCOLOR,WRAPCHARS,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.prompt(n) - timing.plannedOnsets.prompt(n));
            
            %display visualization score
            keymap_prompt = Screen('MakeTexture', mainWindow, keymap_image);
            Screen('DrawTexture',mainWindow,keymap_prompt,[],[],[]); %[0 0 keymap_dims],[topLeft topLeft+keymap_dims]);
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.vis(n), timing.trig.vis_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE,timing.plannedOnsets.vis(n));
            end
            timespec = timing.plannedOnsets.vis(n) - SLACK;
            timing.actualOnsets.vis(n) = Screen('Flip',mainWindow,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.vis(n)-timing.plannedOnsets.vis(n));
            subjectiveEK = easyKeys(subjectiveEK, ...
                'onset', timing.actualOnsets.vis(n), ...
                'stim', stim.stim{stim.trial}, ...
                'cond', stim.cond(stim.trial), ...
                'cresp', cresp, 'cresp_map', cresp_map, 'valid_map', subj_map);
            
           
            %endrecord = GetSecs;
            %display even/odd
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.math(n), timing.trig.math_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE,timing.plannedOnsets.math(n));
            end
            timespec = timing.plannedOnsets.math(n) - SLACK;
            [stim.digitAcc(stim.trial) stim.digitRT(stim.trial) timing.actualOnsets.math(n)] = odd_even(digitsEK,num_digit_qs,digits_promptDur,digits_isi,minimal_format,mainWindow,keyCell([1 5]),COLORS,DEVICE,SUBJ_NAME,[SESSION stim.trial],SLACK,timespec, keys);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.math(n)-timing.plannedOnsets.math(n));
            
            save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        end
        
        %present last ITI
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            [timing.trig.lastITI, timing.trig.lastITI_Success] = WaitTRPulse(TRIGGER_keycode,DEVICE,timing.plannedOnsets.lastITI);
        end
        timespec = timing.plannedOnsets.lastITI - SLACK;
        timing.actualOnsets.lastITI = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.lastITI-timing.plannedOnsets.lastITI);
        
        if GetSecs - timing.actualOnsets.lastITI < 2
            WaitSecs(2 - (GetSecs - timing.actualOnsets.lastITI));
        end
        %wait in scanner at end of run
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            WaitSecs(10);
        end
        
        % clean up
        save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        
        printlog(LOG_NAME,['\n\nSESSION ' int2str(SESSION) ' ended ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        
        if SESSION == RECALL2
            endSession(subjectiveEK,'Congratulations, you have completed the scan! All that is left is a short test outside the scanner. We will come and get you out in just a moment.');
        elseif SESSION == RECALL_PRACTICE
            endSession(subjectiveEK, 'Congratulations, you have completed the practice tasks!');
        else
            endSession(subjectiveEK, CONGRATS)
        end
        sca
           %% 4. Typing description task
    case {DESCRIPTION}
         % stimulus presentation parameters
        stim.promptDur = 6*SPEED; % cue word alone
        stim.prepDur = 2*SPEED;
        stim.isiDuration = 4*SPEED; %
        stim.typeDur = 20*SPEED;
        % stimulus data fields

        stim.description = {};
        PROGRESS_TEXT = 'INDEX';
        
        returnKey=40; %define return key (used to end)
        deleteKey=42; %define delete key (used to delete)
        spaceKey=44;
        periodKey=55;
        quoteKey=52;
        commaKey=54;
        % all the instructions
        stim.instruct1 = ['SCENE DESCRIPTION TASK\n\nWhen a word appears, picture the scene it names as vividly as if you were looking ' ...
            'at it now. Picture its objects, features, and anything else you can imagine. HOLD the image, letting it continue to mature ' ...
            'and take shape for the entire ' num2str(stim.promptDur) ' seconds it is on the screen. Don''t let it fade, and don''t let ' ...
            'yourself "space out". It is essential you actually do this!\n\n' ...
            '-- Press ' PROGRESS_TEXT ' once you understand these instructions --'];
        stim.instruct2 = ['After you have had several seconds to form an image, we will ask you to describe the scene. Your description (minimum of 5-10 words) '...
            'should aim to capture the key elements of the scene so that someone unfamiliar with the scene could reliably select the image corresponding to your description. '...
            'You will have 20 seconds to type your answer. After 15 seconds, you will see a red bar appear to signal that you have 5 seconds left. Please keep describing for all '...
            '20 seconds though.\n\n-- Press ' PROGRESS_TEXT ' once you understand these instructions --'];

        % let's have the first three trials be practice ones
        % prepare counterbalanced trial sequence (at most 2 in a row)
        [stim.cond stim.condString stim.stim] = counterbalance_items({cues{STIMULI}{REALTIME}, cues{STIMULI}{OMIT}},CONDSTRINGS);
        practiceWords = cues{STIMULI}{LEARN}(1:3); % practice words from before is current state?
        stim.stim = [practiceWords stim.stim];
        stim.cond = [PRACTICE PRACTICE PRACTICE stim.cond];
        stim.condString = [CONDSTRINGS{PRACTICE} CONDSTRINGS{PRACTICE} CONDSTRINGS{PRACTICE} stim.condString];
        % initialize stimulus order with initial warmup item
        
        displayText(mainWindow,['SCENE DESCRIPTION TASK\n\nThis is a memory test with some differences ' ...
                'from earlier: you will get only one try per item, there is no feedback, and you will type descriptions of the scenes instead of choosing the correct option.' ...
                'Please carefully review today''s instructions, since many things have ' ...
                'changed and it is important you follow them exactly.\n\n' ...
                '-- Please press ' PROGRESS_TEXT ' to briefly review the instructions --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        waitForKeyboard(kbTrig_keycode,DEVICE);

        % display instructions
        DrawFormattedText(mainWindow,' ','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        displayText(mainWindow,stim.instruct1,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        waitForKeyboard(kbTrig_keycode,DEVICE);
        displayText(mainWindow,stim.instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);   
        waitForKeyboard(kbTrig_keycode,DEVICE);

        %last instructions
        displayText(mainWindow,['To help you get used to the feel of this task, we will ' ...
            'now give you three practice words. Because of the system used, shift keys won''t work so don''t use symbols. ' ...
            ' You can use letters, backspace, and basic punctuation though. If you want to type numbers, type them out as ''four'' and not ''4''. \n\n' ...
            '-- Press ' PROGRESS_TEXT ' to begin once you understand these instructions --'], ...
            minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);

        runStart = GetSecs;
        config.wait = 0;

        % want: cued word then window to type description
        % then math, ISI
        % for all stimuli
        config.TR = stim.TRlength;
        config.nTRs.ISI = stim.isiDuration/stim.TRlength;
        config.nTRs.prompt = stim.promptDur/stim.TRlength;
        config.nTRs.type = stim.typeDur/stim.TRlength;
        config.nTrials = length(stim.stim);
        config.nTRs.perTrial = (config.nTRs.ISI + config.nTRs.prompt + config.nTRs.type);
        config.nTRs.perBlock = config.wait/config.TR + (config.nTRs.perTrial)*config.nTrials+ config.nTRs.ISI; %includes the last ISI
        
        % calculate all future onsets
        timing.plannedOnsets.preITI(1:config.nTrials) = runStart + config.wait + ((0:config.nTrials-1)*config.nTRs.perTrial)*config.TR;
        timing.plannedOnsets.prompt(1:config.nTrials) = timing.plannedOnsets.preITI + config.nTRs.ISI*config.TR;
        timing.plannedOnsets.type(1:config.nTrials) = timing.plannedOnsets.prompt + config.nTRs.prompt*config.TR;
        timing.plannedOnsets.lastITI = timing.plannedOnsets.type(end) + config.nTRs.type*config.TR;
        
     
        stimCond = stim.cond;
        % initialize trial
        %generate stimulus ID's first so can add them easily
        for i = 1:length(stim.cond)
            pos = find(strcmp(preparedCues,stim.stim{i}));
            if ~isempty(pos)
                stim.id(i) = pos;
            else
                stim.id(i) = -1; %so this should never go during MOT
            end
        end
        stimID = stim.id;
        sessionInfoFile = fullfile(ppt_dir, ['SessionInfo' '_' num2str(SESSION)]);
        save(sessionInfoFile, 'stimCond','stimID', 'timing', 'config'); 
        

        for n = 1:length(stim.stim)
            % initialize trial and show cue
            stim.trial = n;
            fprintf(['Trial number: ' num2str(n) '\n']);
            
            %show pre ITI
            timespec = timing.plannedOnsets.preITI(n) - SLACK;
            timing.actualOnsets.preITI(n) = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.preITI(n) - timing.plannedOnsets.preITI(n));
            
            %display word
            timespec = timing.plannedOnsets.prompt(n)-SLACK;
            timing.actualOnsets.prompt(n) = displayText_specific(mainWindow,stim.stim{stim.trial},'center',COLORS.MAINFONTCOLOR,WRAPCHARS,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.prompt(n) - timing.plannedOnsets.prompt(n));
            
            % have them type and then get the rating
            timespec = timing.plannedOnsets.type(n) - SLACK;
            instructions = 'Please describe the image.';
            Screen('TextSize',mainWindow,20); %sets textsize for instructions
            [nxi,nyi,textbox_i] = DrawFormattedText(mainWindow,instructions, 'center', CENTER(2) - CENTER(2)/3, COLORS.MAINFONTCOLOR);
            new = textbox_i;
            warning = textbox_i;
            warning(2) = warning(2) + 20;
            warning(4) = warning(4) + 20;
            movedown = 80;
            moveaway = 80;
            %new(1) = new(1) - moveaway;
            new(2) = new(2) + movedown;
            boxsize = 200;
            new(4) = new(4)+ boxsize + movedown;
            %new(3) = new(3) + moveaway; %making the textbox a bit wider
            new(1) = CENTER(1) - CENTER(1)/1.5;
            new(3) = CENTER(1) + CENTER(1)/1.5;
            Screen('FillRect', mainWindow, COLORS.GREY/5, new);
            timing.actualOnsets.type(n) = Screen('Flip',mainWindow, timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.type(n) - timing.plannedOnsets.type(n));
           
            KbQueueCreate;
            KbQueueStart;
            AsteriskBuffer=[]; %initializes buffer
            WRAPCHARS = 100;

            AsteriskBuffer = ' ';
            while GetSecs - timing.actualOnsets.type(n) < stim.typeDur % keep checking for more typing until you can't anymore
                [ pressed, firstPress]=KbQueueCheck; %checks for keys
                if (GetSecs- timing.actualOnsets.type(n)) > 15 % draw red bar if longer than 15 seconds has passed
                    Screen('FillRect', mainWindow, COLORS.RED, warning)
                end
                Screen('TextSize',mainWindow,20);  %sets textsize for instructions
                [nxi,nyi] = DrawFormattedText(mainWindow,instructions, 'center', CENTER(2) - CENTER(2)/3, COLORS.MAINFONTCOLOR);
                Screen('FillRect', mainWindow, COLORS.GREY/5, new)
                Screen('TextSize',mainWindow,25);
                
                if pressed %keeps track of key-presses and draws text
                    if firstPress(deleteKey) && length(AsteriskBuffer>1) %if delete key then erase last key-press
                        AsteriskBuffer=AsteriskBuffer(1:end-1); %erase last key-press
                    elseif firstPress(spaceKey)
                        % we want to add a space instead of writing space lol
                        AsteriskBuffer=[AsteriskBuffer ' '];
                    elseif firstPress(periodKey)
                        AsteriskBuffer=[AsteriskBuffer '.'];
                    elseif firstPress(quoteKey)
                        AsteriskBuffer=[AsteriskBuffer  '"' ];
                    elseif firstPress(commaKey)
                        AsteriskBuffer=[AsteriskBuffer  ',' ];
                    else %otherwise add to buffer
                        firstPress(find(firstPress==0))=NaN; %little trick to get rid of 0s
                        [endtime Index]=min(firstPress); % gets the RT of the first key-press and its ID
                        AsteriskBuffer=[AsteriskBuffer KbName(Index)]; %adds key to buffer
                    end
                    if isempty(AsteriskBuffer)
                        AsteriskBuffer = ' ';
                    end
                end
                [nx,ny,textbounds] = DrawFormattedText(mainWindow, AsteriskBuffer, new(1),new(2),COLORS.MAINFONTCOLOR,WRAPCHARS); %it's going where x ends and y starts
                %have it so it goes to the next line when they type the next line
                Screen('Flip',mainWindow);
                WaitSecs(.01); % put in small interval to allow other system events
            end
            Screen('TextSize',mainWindow,20); 
            stim.description{n} = AsteriskBuffer;
            %endrecord = GetSecs;
            %display even/odd
            
            if n <= 3
                % WaitSecs(1);
                OnFB = GetSecs;
                if strcmp(AsteriskBuffer,' ')
                    displayText(mainWindow,['Oops, we didn''t record a response. Please let the experimenter know if you have any questions, or just make sure to type a response in time on the remainder of trials. ' ...
                        '\n\n-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,...
                        'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                else
                    displayText(mainWindow,['Good work! Your response was detected.\n\n-- Press ' ...
                        PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR, ...
                        WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                end
                if n==3
                    displayText(mainWindow,['We will now move onto the main experiment. Please try your best for all 20 seconds! ' ...
                        '\n\n-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,...
                        'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                end
                OffFB = GetSecs;
                timing.plannedOnsets.preITI(n+1:end) = timing.plannedOnsets.preITI(n+1:end) + OffFB-OnFB;
                timing.plannedOnsets.prompt(n+1:end) = timing.plannedOnsets.prompt(n+1:end) + OffFB-OnFB;
                timing.plannedOnsets.type(n+1:end) = timing.plannedOnsets.type(n+1:end) + OffFB-OnFB;
                timing.plannedOnsets.lastITI = timing.plannedOnsets.lastITI + OffFB-OnFB;
            end
            
            KbQueueRelease;            
            save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        end
        
        %present last ITI
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            [timing.trig.lastITI, timing.trig.lastITI_Success] = WaitTRPulse(TRIGGER_keycode,DEVICE,timing.plannedOnsets.lastITI);
        end
        timespec = timing.plannedOnsets.lastITI - SLACK;
        timing.actualOnsets.lastITI = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.lastITI-timing.plannedOnsets.lastITI);
        
        if GetSecs - timing.actualOnsets.lastITI < 2
            WaitSecs(2 - (GetSecs - timing.actualOnsets.lastITI));
        end

        % clean up
        save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        
        printlog(LOG_NAME,['\n\nSESSION ' int2str(SESSION) ' ended ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        
        DrawFormattedText(mainWindow, CONGRATS, 'center', 'center')
        Screen('Flip',mainWindow);
        WaitSecs(4);
        
        % clean up stuff loaded in memory, and restore user input
        %KbQueueRelease();
        sca;
        
        %% POST PICTURES TASK
        
    case ASSOCIATES
        
        % declarations
        stim.promptDur = .75*SPEED;
        stim.listenDur = 0.5*SPEED;
        stim.isiDuration = 2*SPEED;
        PROGRESS = INDEXFINGER;
        PROGRESS_TEXT = 'INDEX';
        condmap = makeMap({'realtime','omit','lure'});
        destDims = [200 200];
        stim.gapWidth = 40;
        stim.picRow = WINDOWSIZE.pixels(2) *(5/9);
        % change it to be index for left, middle for right choice
        nPractice = 3; % number of practice trials
        
        stim.instruct1 = ['SCENE MEMORY LIGHTNING ROUND\n\nYou''re almost done! This is the final task.\n\nWe will show you pictures of various scenes and ask you ' ...
            'whether they are old ("' recog_scale.inputs{1} '" key) or new ("' recog_scale.inputs{2} '" key). Try to respond ' ...
            'as quickly and accurately as possible.\n So press ' recog_scale.inputs{1} ' = OLD IMAGE and ' recog_scale.inputs{2} ' = NEW IMAGE.\n\n-- Press ' PROGRESS_TEXT ' once you understand these instructions --'];
        
        stim.instruct2 = ['To get you started, we''re going to give you three practice trials. After this, we will move on to the task. '...
            ' Remember to press ' recog_scale.inputs{1} ' = OLD IMAGE and ' recog_scale.inputs{2} ' = NEW IMAGE. \n\n-- Press ' PROGRESS_TEXT ' once you understand these instructions --'];
        % stimulus data fields
        stim.triggerCounter = 1;
        stim.missedTriggers = 0;
        
  
        % prepare stimuli: first actual ones that were seen
        
        COMPLETED = 0;
        while ~COMPLETED
            [stim.cond stim.condString stim.associate] = counterbalance_items({cues{STIMULI}{REALTIME}, cues{STIMULI}{OMIT}},CONDSTRINGS);
            [stim.cond2 stim.condString2 stim.associate2] = counterbalance_items({cues{STIMULI}{REALTIME}, cues{STIMULI}{OMIT}},CONDSTRINGS);
            for n = 1:length(stim.cond)
                cueSearch = strcmp(preparedCues,stim.associate{n});
                stim.pos1(n) = find(cueSearch);
                stim.stim1{n} = pics{stim.pos1(n)};
                stim.id1(n) = find(strcmp(pics,stim.stim1{n}));
                
                cueSearch = strcmp(preparedCues,stim.associate2{n});
                stim.pos2(n) = find(cueSearch);
                stim.stim2{n} = pics{stim.pos2(n)};
                stim.id2(n) = find(strcmp(pics,stim.stim2{n}));
            end
            allImages = [stim.id1 stim.id2];
            for n = 1:length(stim.cond)
                thisImage = allImages(n);
                nextImage = max(find(thisImage == allImages));
                diffinpos(n) = abs(n-nextImage);
            end
            nokay = length(find(diffinpos >=10));
            if nokay == length(stim.cond)
                COMPLETED =1;
            end
        end
        
        % now choose which ones come AA' and which are A'A
        % choose 5 omit images
        OMIT_INDS = find(stim.cond==OMIT);
        omit_AAP = OMIT_INDS(randperm(10,5));
        REALTIME_INDS = find(stim.cond==REALTIME);
        rt_AAP = REALTIME_INDS(randperm(10,5));
        stim.AAP = [omit_AAP rt_AAP]; % for the first round of all 20 these will have A first
        stim.AAPID = stim.id1(stim.AAP); % these are the stimulus id numbers to present the A first!
        % Get practice images--don't worry about balancing indoor/outdoor
        % here
        stim.cond = [stim.cond stim.cond2];
        stim.stim = [stim.stim1 stim.stim2];
        stim.id = [stim.id1 stim.id2];
        stim.associate = [stim.associate stim.associate2];
        
        stim.stim = [Shuffle(pics(1:nPractice)) stim.stim]; % practice words from before is current state?
        stim.id = [-1 -1 -1 stim.id];
        stim.cond = [PRACTICE PRACTICE PRACTICE stim.cond];
        stim.condString = [CONDSTRINGS{PRACTICE} CONDSTRINGS{PRACTICE} CONDSTRINGS{PRACTICE} stim.condString stim.condString2];
        stim.associate = [CONDSTRINGS{PRACTICE} CONDSTRINGS{PRACTICE} CONDSTRINGS{PRACTICE} stim.associate];
        % get practice indoor/outdoor

       
        % display instructions
        DrawFormattedText(mainWindow,' ','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        displayText(mainWindow,stim.instruct1,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        waitForKeyboard(kbTrig_keycode,DEVICE);
        displayText(mainWindow,stim.instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);
        
        % initialize question
        triggerNext = false; use_console = true; make_files = true;
        stimmap = makeMap(pics);
        recogEK = initEasyKeys([exp_string_long '_RECOG'], SUBJ_NAME, ppt_dir,...
            'default_respmap', recog_scale, ...
            'stimmap', stimmap, ...
            'condmap', condmap, ...
            'trigger_next', triggerNext, ...
            'prompt_dur', stim.promptDur, ...
            'listen_dur', stim.listenDur, ...
            'exp_onset', stim.subjStartTime, ...
            'device', DEVICE);
        recogEK = startSession(recogEK);
        runStart = GetSecs;
        
        config.TR = stim.TRlength;
        config.nTRs.ISI = stim.isiDuration/stim.TRlength;
        config.nTRs.cue = (stim.promptDur + stim.listenDur)/stim.TRlength;
        config.nTrials = length(stim.cond);
        config.nTRs.perTrial = config.nTRs.ISI + config.nTRs.cue;
        config.nTRs.perBlock = (config.nTRs.perTrial) * config.nTrials + config.nTRs.ISI;
        
        timing.plannedOnsets.preITI(1:config.nTrials) = runStart + ((0:config.nTrials-1)*config.nTRs.perTrial)*config.TR;
        timing.plannedOnsets.cue(1:config.nTrials) = timing.plannedOnsets.preITI + config.nTRs.ISI*config.TR;
        timing.plannedOnsets.lastITI = timing.plannedOnsets.cue(end) + config.nTRs.cue*config.TR;
        % repeat
        lureProgress = 0;
        CHOICES = length(recog_scale.inputs);
        for n = 1:length(stim.cond)
            
            timespec = timing.plannedOnsets.preITI(n) - SLACK;
            timing.actualOnsets.preITI(n) = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.preITI(n) - timing.plannedOnsets.preITI(n));
            

            
            % now present the target and lure and choose where to put them!
            if n <= nPractice% randomize if lure or not
                lure = randi(2);
                if lure == 1
                    picIndex = prepImage(strcat(TRAININGFOLDER,stim.stim{n}),mainWindow);
                    cresp = recog_scale.inputs(1);
                    cresp_map = keys.map(2,:);
                    recog.cresp_string{n} = 'old';
                else
                    picIndex = prepImage(strcat(TRAININGLURESFOLDER, stim.stim{n}),mainWindow);
                    cresp = recog_scale.inputs(2);
                    cresp_map = keys.map(3,:);
                    recog.cresp_string{n} = 'new';
                end
            else
                if ismember(stim.id(n),stim.AAPID)
                    if n <=20 + nPractice
                        picIndex = prepImage(strcat(MOTFOLDER, stim.stim{n}),mainWindow);
                        cresp = recog_scale.inputs(1);
                        cresp_map = keys.map(2,:);
                        recog.cresp_string{n} = 'old';
                    else
                        picIndex = prepImage(strcat(MOTLURESFOLDER, stim.stim{n}),mainWindow);
                        cresp = recog_scale.inputs(2);
                        cresp_map = keys.map(3,:);
                        recog.cresp_string{n} = 'new';
                    end
                else
                    if n <= 20 + nPractice
                        picIndex = prepImage(strcat(MOTLURESFOLDER, stim.stim{n}),mainWindow);
                        cresp = recog_scale.inputs(2);
                        cresp_map = keys.map(3,:);
                        recog.cresp_string{n} = 'new';
                    else
                        picIndex = prepImage(strcat(MOTFOLDER, stim.stim{n}),mainWindow);
                        cresp = recog_scale.inputs(1);
                        cresp_map = keys.map(2,:);
                    	recog.cresp_string{n} = 'old';
                    end
                end
            end

            % now present the target
            topLeft(HORIZONTAL) = CENTER(HORIZONTAL) - (PICDIMS(HORIZONTAL)*RESCALE_FACTOR/2);
            topLeft(VERTICAL) = CENTER(VERTICAL) - (PICDIMS(VERTICAL)*RESCALE_FACTOR/2);
            Screen('DrawTexture', mainWindow, picIndex, [0 0 PICDIMS],[topLeft topLeft+PICDIMS*RESCALE_FACTOR]);
            timespec = timing.plannedOnsets.cue(n) - SLACK; %first one will be .4 off but it's because it's not expecting there to be an etra .4 seconds
            timing.actualOnsets.cue(n) = Screen('Flip',mainWindow,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.cue(n) - timing.plannedOnsets.cue(n));
            Screen('Close',picIndex);
            
            % free recall judgment
            DrawFormattedText(mainWindow,'+','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            recogEK = easyKeys(recogEK, ...
                'onset', timing.actualOnsets.cue(n), ...
                'stim', stim.stim{n}, ...
                'cond', stim.cond(n), ...
                'cresp', cresp, ...
                'next_window', mainWindow, 'cresp_map', cresp_map, 'valid_map', recog_map );
            
            
            if n <= 3
                % WaitSecs(1);
                OnFB = GetSecs;
                if isnan(recogEK.trials.resp(end))
                    displayText(mainWindow,['Oops, you either did not respond in time, or did not ' ...
                        'press an appropriate key. Remember to respond with "' recog_scale.inputs{1} '" key for OLD IMAGE and "' recog_scale.inputs{2} '" ' ...
                        'key for a NEW IMAGE. Try to respond as quickly and accurately as ' ...
                        'possible.\n\n-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,...
                        'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                else
                    displayText(mainWindow,['Good work! Your response was detected.\n\n-- Press ' ...
                        PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR, ...
                        WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                end
                if n==3
                    displayText(mainWindow,['We will now move onto the main experiment. Remember to respond with "' recog_scale.inputs{1} '" key for OLD IMAGE and "' recog_scale.inputs{2} '" ' ...
                        'key for a NEW IMAGE. Try to respond as quickly and accurately as ' ...
                        'possible.\n\n-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,...
                        'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                    waitForKeyboard(kbTrig_keycode,DEVICE);
                end
                OffFB = GetSecs;
                timing.plannedOnsets.preITI(n+1:end) = timing.plannedOnsets.preITI(n+1:end) + OffFB-OnFB;
                timing.plannedOnsets.cue(n+1:end) = timing.plannedOnsets.cue(n+1:end) + OffFB-OnFB;
                timing.plannedOnsets.lastITI = timing.plannedOnsets.lastITI + OffFB-OnFB;
            end
            
            % rinse, save and repeat
            save(MATLAB_SAVE_FILE,'stim','recog', 'timing', 'config');
        end
        timespec = timing.plannedOnsets.lastITI - SLACK;
        timing.actualOnsets.lastITI = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.lastITI-timing.plannedOnsets.lastITI);
        WaitSecs(stim.isiDuration);
        
        save(MATLAB_SAVE_FILE,'stim','recog', 'timing', 'config');
        % clean up
        
        printlog(LOG_NAME,['\n\nSESSION ' int2str(SESSION) ' ended ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        printlog(LOG_NAME,'\n\n\n******************************************************************************\n');
        endSession(recogEK, 'Congratulations, you''ve finished the experiment! Please contact your experimenter.');
        sca;
        %% MOT
    case [MOT_PREP MOT MOT_PRACTICE MOT_PRACTICE2 MOT_LOCALIZER]
        
        
        if SESSION == MOT_PRACTICE2
            displayText(mainWindow,['Welcome to your fMRI scanning session!\n\nOnce you''re all the way inside the scanner and can read this text, please reach up to your eyes and ' ...
                'fine-tune the position of your mirror. You want to set it so you can see as much of the screen as comfortably as possible. This will be your last chance to adjust ' ...
                'your mirror, so be sure to set it just right.\n\nOnce you''ve adjusted the mirror to your satisfaction, please press the index finger button to test your button pad.'] ...
                ,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,'Great. I detected that button press, which means at least one button works. Now let''s try the rest of them. Please press the middle finger button.' ...
                ,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(keys.code(3,:),DEVICE);
            displayText(mainWindow,'And now the ring finger button...',minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(keys.code(4,:),DEVICE);
            displayText(mainWindow,'And now the pinky finger button...',minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(keys.code(5,:),DEVICE);
            displayText(mainWindow,'And now the thumb button...',minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(keys.code(1,:),DEVICE);
            displayText(mainWindow,['Good news! It looks like the button pad is working just fine.\n\nJust a reminder that we can hear your voice when the scanner is at rest, so ' ...
                'just speak up to ask or tell us something. During a scan, we will abort right away if you use the squeeze ball, but please do so only if there''s something urgent ' ...
                'we need to address immediately.\n\n-- please press the index finger button to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,['During the scan today, it is crucial that you keep your head still. Even a tiny head movement, e.g., caused by stretching your legs, will blur ' ...
                'your brain scan. This is for the same reason that moving objects appear blurry in a photograph.\n\nAs it can be uncomfortable to stay still for a long time, please ' ...
                'go ahead and take the opportunity to stretch or scratch whenever the scanner is silent. Just try your best to keep your head in the same place when you do so.\n\n' ...
                '-- please press the index finger button to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,['We''re now going to start with a five-minute anatomical scan as well as a couple other example scans while you complete some training tasks in preparation for later. Please work through these and we''ll get in ' ...
                'touch with you when you finish.\n\n-- please press the index finger button to continue --'],INSTANT,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
        end

        
        % stimulus presentation parameters
        instant = 0;
        num_dots = 5; %where set total dots
        stim.square_dims = round([20 20] ./ mean(WINDOWSIZE.degrees_per_pixel)); % 20ï¿½ visual angle in pixels
        stim.dot_diameter = 1.5 / mean(WINDOWSIZE.degrees_per_pixel); % 1.5 visual angle in pixels
%         if SESSION < MOT{1}
%             stim.trialDur = 20*SPEED;
%             promptTRs = 3:2:9;
%         else
        stim.trialDur = 30*SPEED; %chaning length from 20 to 30 s or 15 TRs
        promptTRs = 3:3:13; %which TR's to make the prompt active- %TO BE CHANGED-used to be 3:3:13; with 15 2 s TRs--now will have 30 TR's
  %      end
        stim.inter_prompt_interval = 4*SPEED;
        stim.maxspeed = 30;
        stim.minspeed = -30; %0.3; %changed to -30 to go into function--now will allow negative speeds and result in them being faded to black
        stim.targ_dur = 2*SPEED;
        stim.fixBlock = 20; %time in the beginnning they're staring--20s, 10 TRs 
        %stim.betweenPrompt = 2;
        %stim.beforePrompt = 4;
        vis_promptDur = 2*SPEED;
        digits_promptDur = 1.9*SPEED;
        digits_isi = 0.1*SPEED;
        PROGRESS_TEXT = 'INDEX finger';
        probe_promptDur = 2*SPEED;
        probe_listenDur = 0;
        num_digit_qs = 2;
        stim.isiDuration = 2*SPEED;%10*SPEED - stim.TRlength*2 - num_digit_qs*(digits_promptDur + digits_isi);
        stim.mathDur = (digits_promptDur + digits_isi) * num_digit_qs;
        stim.fbDuration = 2*SPEED; %for dot probe choice
        dot_map = keys.map([1 5],:);
        
        %rt parameters 
        OptimalForget = 0.15;
        maxIncrement = 1.25; %will also have to check this
        Kp = 5;
        Ki = .0; %changed after RT worksop from 0.01
        Kd = .5;
        config.initFeedback = OptimalForget; %make it so there's change in speed for the first 4 TR's
        config.initFunction = PID(config.initFeedback,Kp,Ki,Kd,OptimalForget,maxIncrement);
        
        switch SESSION
            case {MOT_PREP,MOT_PRACTICE}
                mot_conds = {'1_targ'};
                show_words = true;
                minimal_format = false;
                day_2 = false;
                realtime = false;
            case {MOT_PRACTICE2,MOT_LOCALIZER}
                mot_conds = {'targ-hard','targ-easy','lure-hard','lure-easy'};
                show_words = true;
                minimal_format = true;
                day_2 = true;
                realtime = false;
            otherwise
                mot_conds = {'realtime'};
                show_words = true;
                minimal_format = true;
                day_2 = true;
                realtime = true; %use this to make other conditions below!
        end
        
        if SESSION == MOT_PRACTICE || SESSION == MOT_PRACTICE2
            stim.header = 'MULTI-TASKING -- PRACTICE';
        elseif SESSION < MOT_PRACTICE2
            stim.header = 'MULTI-TASKING';
        else stim.header = 'MULTI-TASKING';
        end
        instruct_continue = ['\n\n-- Press ' PROGRESS_TEXT ' to continue once you understand these instructions --'];
        stim.instruct1 = [stim.header '\n\nWe will now do a "multi-tasking" twist: we would like you to ' ...
            'try visualizing named scenes while also keeping track of moving dots. One word will cue an image to visualize in the center of the screen, with five dots surrounding it. \n\nDot-tracking works as follows: first, one dot ' ...
            'will appear in red (target) and others in green (non-targets). After two seconds, all dots will turn ' ...
            'green and and move around the screen. Your job is to track the target dot until a "question" dot turns white. ' ...
            'While tracking dots, it is VITAL that you keep your eyes at the center fixation dot for the entire '...
            'time that you are tracking dot motion!!! When the ''question'' dot turns white, you '...
            'can then take your eyes off the center and indicate whether that dot was originally a target (PINKY) or not (THUMB). ' ...
            instruct_continue];
        stim.instruct2 = ['While all this happens, we want you to "multi-task" by visualizing the scene ' ...
            'named by the word in the middle. Every few seconds, the smaller central dot inside the word will turn red: when it does '...
            'you should indicate your CURRENT visualization detail using one of your fingers (THUMB,INDEX,MIDDLE,RING,PINKY). It''s possible that the clarity of your ' ...
            'mental image will change throughout the trial and/or be fuzzier than if you were not tracking dots. As your rating should reflect your CURRENT ' ...
            'rather than best possible image, this would be reflected in your responses. \n\nRemember that mental visualizing is only your second priority: if you lose track ' ...
            'of the dots for even a second, you will get the trial wrong, and we will have to throw out the trial; ' ...
            'and we need as many trials as possible. Getting the target dot correct is your most important task.' ...
            'You should "squeeze in" visualizing when it''s possible. \n\n---- Press ' PROGRESS_TEXT ' once you understand these instructions \nthen PRESS IT AGAIN when you are done viewing the rating scale ---'];
        stim.instruct_summary = [stim.header '\n\nTo summarize this task: you will keep track of target dots moving around the screen while keeping your ' ...
            'eyes fixed on the central dot. The dot-tracking task is your top priority, but you should also try to visualize ' ...
            'the named scene, and rate your mental picture whenever the center dot turns red (with the same keys THUMB-PINKY). \n\nThe speed of the dots may change in '...
            'different trials and also within a trial. Just try to do your best and keep doing the task no matter the dot speed. '...
            'At the end of the trial, when a dot ' ...
            'turns white, you will first move your eyes to the white dot, then press PINKY (if it''s a target) or THUMB (if it''s ' ...
            'not). (Think YES TARGET = PINKY, NOT TARGET = THUMB.) \n\nThere will also be several even/odd questions after each trial that you should try to complete (using ' ...
            'your THUMB for even and PINKY for odd).' final_instruct_continue];
        stim.fMRI_refresher = ['MULTI-TASKING-- fMRI PRACTICE\n\nDot tracking today will be similar to last time, with two ' ...
            'changes. Firstly, dots in half of trials will be moving slow, while dots in the other half will be moving faster.\n\n' ...
            'Secondly, some of the trials will involve familiar words that are not scene names. On these trials, ' ...
            'because the word is not the name of a scene, there is nothing for you to visualize, so press to indicate low visualization.\n\nWe will now review the ' ...
            'instructions for the task (with which you are already familiar).' instruct_continue];
        stim.fMRI_instruct = ['MULTI-TASKING-- fMRI, Stage 1 \n\nDot tracking today will involve two changes from yesterday. ' ...
            'Firstly, dots in half of trials will be moving faster, while dots in the other half will be moving slower. Additionally, the dots that are not the target '...
            'may become dimmer to help you track the target.\n' ...
            'Secondly, some of the trials will involve familiar words that are not scene names. On these trials, ' ...
            'because the word is not the name of a scene, there is nothing for you to visualize.\n\nWe will now review the ' ...
            'instructions for the task (with which you are already familiar).' instruct_continue];
        stim.RT_instruct = ['MULTI-TASKING-- fMRI, Stage 2 \n\n You will now do the same dot tracking task as earlier. However, '...
            'now every word will have been paired with a scene. Additionally, the dot speed may change within '...
            'a trial or the dots that are not the target may become dimmer to help you track the target. Just try your best to keep tracking the target dot while visualizing the scene as best as you can! You will complete '...
            'three runs of this task.' instruct_continue];
        % stimulus data fields
        stim.triggerCounter = 1;
        stim.missedTriggers = 0;
        stim.normalSize = 36;
        stim.smallSize = 24;
        stim.bigSize = 100;
        stim.keys = keyCell([1 5]);
        queueCheck = -1;
        
        % initilize staircasing
        if SESSION == MOT_PREP
            stair = 1;
        else
            stair = 0;
        end
        tGuess = 15; % guess as to the blend that will give us the our pThreshold
        pThreshold = .85; % probability of a resp=1 that we are aiming for
        beta = 3; % steepness of psychometric function, typically 3
        delta = .1; % fraction of trials on which observer presses blindly
        gamma = .5; % response rate when intensity = 0, typically .5 for forced-choiced questions w/ 2 possible answers
        grain = .1; %changed from 0.02 to .1--too little too step size
        range = 30; % changed from 5-- supposed to be a generous range of possible intensities. centered around initial tGuess
        % during prep, we free up the staircasing algorithm--here we're
        % setting parameters
        if stair
            tGuessSd = 5;
            %on day 2, use params from prior session
        elseif ~stair && SESSION > MOT_PREP %take speed and stop updating afterwards
            matlabOpenFile = []; timeout = true;
            %find_str = [ppt_dir 'stable_subj_' num2str(SUBJECT) '_' num2str(MOT_PREP) '*.mat']; %only look in MOT prep--no more updating
            relative_paths = false;
            %             while timeout
            try
                %will need to add this file to the folder to get it to work
                fileCandidates = dir([ppt_dir 'mot_realtime05_' num2str(SUBJECT) '_' num2str(MOT_PREP)  '*.mat']);
                dates = [fileCandidates.datenum];
                names = {fileCandidates.name};
                [~,newest] = max(dates);
                matlabOpenFile = [ppt_dir names{newest}];
                lastRun = load(matlabOpenFile);
                
            catch
                fprintf(['Could not find a prior run file. Using speed = 15.']);
                pause(1);
                lastRun.stim.tGuess = 15;
            end
            %             end
            finalSpeed = stim.maxspeed - lastRun.stim.tGuess(end);
            %finalSpeed = 15;
        end % stair
        if stair
            questStruct = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range); %Quest function for face blocks
        end
        keymap_image = imread(KEY_MAPPING);
        keymap_prompt = Screen('MakeTexture', mainWindow, keymap_image);
        
        % allocate stimuli
        if ~day_2
            stim.lureWords = [];
            if SESSION == MOT_PRACTICE
                [stim.cond stim.condString stim.stim] = counterbalance_items({cues{STIMULI}{LEARN}},{MOTSTRINGS{LEARN}},0);
            else %MOT_PREP
                [stim.cond stim.condString stim.stim] = counterbalance_items({[cues{STIMULI}{LEARN} cues{STIMULI}{LEARN} cues{STIMULI}{LEARN} cues{STIMULI}{LEARN} cues{STIMULI}{LEARN} cues{STIMULI}{LEARN}]},{MOTSTRINGS{LEARN}},0);
            end
            condmap = makeMap({'target'});
        elseif SESSION==MOT_PRACTICE2
            stim.lureWords = lureWords(6:7);
            [stim.cond stim.condString stim.stim] = counterbalance_items({cues{STIMULI}{LEARN}(1) cues{STIMULI}{LEARN}(2), stim.lureWords(1), stim.lureWords(2)},MOTSTRINGS,1);
            condmap = makeMap({'targ_hard','targ_easy','lure_hard','lure_easy'});
            square_bounds = [CENTER-(stim.square_dims/2) CENTER+(stim.square_dims/2)-1];
            stim.condString = condmap.descriptors(stim.cond);
        elseif SESSION == MOT_LOCALIZER
            stim.lureWords = lureWords(8:23);
            [stim.cond stim.condString stim.stim] = counterbalance_items({cues{STIMULI}{LOC}(1:8), cues{STIMULI}{LOC}(9:16), stim.lureWords(1:8), stim.lureWords(9:16)},MOTSTRINGS,1); %looks like he wanted to separate easy/hard conditions
            condmap = makeMap({'targ_hard','targ_easy','lure_hard','lure_easy'});
            square_bounds = [CENTER-(stim.square_dims/2) CENTER+(stim.square_dims/2)-1];
            stim.condString = condmap.descriptors(stim.cond);
        else
            %stim.lureWords = lureWords(1:5); %don't need these anymore but oh well
            [stim.cond stim.condString stim.stim] = counterbalance_items({cues{STIMULI}{REALTIME}}, MOT_RT_STRINGS,1);
            condmap = makeMap({'rt-targ'});
        end
        square_bounds = [CENTER-(stim.square_dims/2) CENTER+(stim.square_dims/2)-1];
        stim.condString = condmap.descriptors(stim.cond);
        stim.repulse = 2;
        % assign parameters based on condition
        lureCounter = 0;
        stim.num_targets = 1;
        repulsor_force_small = 1;
        
        %generate stimulus ID's first so can add them easily
        for i = 1:length(stim.cond)
            pos = find(strcmp(preparedCues,stim.stim{i}));
            if ~isempty(pos)
                stim.id(i) = pos;
            else
                stim.id(i) = -1; %so this should never go during MOT
            end
        end
        
        
        for i=1:length(stim.cond)
            
            if ~day_2 && ~stair
                stim.speed(i) = stim.maxspeed - tGuess; %setting practice speed to initial tGuess
                %stim.speed(i) = -3*i; this was to check about the dots
                %being dark
                %stim.repulse(i) = 5/3;
            else
                if SESSION > 5 && SESSION < MOT{1} %change dot speeds for practice and localizer, but then we want to change dot speed!
                    switch stim.cond(i)
                        case {1,3} %for either of the hard cases
                            stim.speed(i) = finalSpeed; %will have to load last speed and find the speed here
                            % repulsor_force(i) = repulsor_force_small * finalSpeed/0.5;
                        case {2,4} %for either of the easy cases
                            stim.speed(i) = -7.5;%5;%finalSpeed*.5; %again load last speed found here, change to accept max speed-see what first person has for this to decide
                            %repulsor_force(i) = repulsor_force_small;
                    end
                    if stim.cond(i) > 2
                        lureCounter = lureCounter + 1;
                    end
                    stim.condString{i} = mot_conds{stim.cond(i)};
                else %for MOT real-time trials
                    % repulsor_force(i) = stim.speed;
                    stim.speed(i) = 2; %initialize each trial here!!!
                end
            end
            
        end
        stim.num_lures = num_dots - stim.num_targets;
        
        % present instructions: change instructions from MOT PRACTICE TO
        % MOT PREP
        stim.stairInstruct = ['We will now continue with the same task you just did, but for more trials. '...
            'We are going to be adjusting the dot speed in this session, so it is very important for this session that you do not guess randomly. '...
            'If you do fall asleep or weren''t paying attention, do NOT guess randomly; instead '...
            'don''t respond at all! Random guessing will throw off our algorithm.\n\nWe will now review the ' ...
            'instructions for the task (with which you are already familiar).' instruct_continue];
        if SESSION == MOT_PRACTICE
            displayText(mainWindow,stim.instruct1,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,stim.instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            keymap_image = imread(KEY_MAPPING);
            keymap_prompt = Screen('MakeTexture', mainWindow, keymap_image);
            Screen('DrawTexture',mainWindow,keymap_prompt,[],[],[]); %[0 0 keymap_dims],[topLeft topLeft+keymap_dims]);
            Screen('Flip',mainWindow);
            waitForKeyboard(kbTrig_keycode,DEVICE);
        elseif SESSION == MOT_PREP
            displayText(mainWindow,stim.stairInstruct,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            
        elseif SESSION == MOT_PRACTICE2
            displayText(mainWindow,stim.fMRI_refresher,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,stim.instruct1,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            displayText(mainWindow,stim.instruct2,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
            keymap_image = imread(KEY_MAPPING);
            keymap_prompt = Screen('MakeTexture', mainWindow, keymap_image);
            Screen('DrawTexture',mainWindow,keymap_prompt,[],[],[]); %[0 0 keymap_dims],[topLeft topLeft+keymap_dims]);
            Screen('Flip',mainWindow);
            waitForKeyboard(kbTrig_keycode,DEVICE);
        elseif SESSION == MOT_LOCALIZER 
            displayText(mainWindow,stim.fMRI_instruct,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
        elseif SESSION == MOT{1}
            displayText(mainWindow,stim.RT_instruct, minimumDisplay, 'center', COLORS.MAINFONTCOLOR,WRAPCHARS);
            waitForKeyboard(kbTrig_keycode,DEVICE);
        end
        if SESSION == MOT{2} || SESSION == MOT{3}
            stim.instruct_nextMOT = ['You will now continue with the same multitasking task.' final_instruct_continue];
            displayText(mainWindow,stim.instruct_nextMOT,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            %also load in last session information here
            allLast = findNewestFile(ppt_dir,[ppt_dir 'mot_realtime05_' num2str(SUBJECT) '_' num2str(SESSION-1) '*']);
            last = load(allLast);
            lastSpeed = last.stim.lastSpeed; %matrix of motRun (1-3), stimID
            lastDecoding = last.stim.lastRTDecoding;
            lastDecodingFunction = last.stim.lastRTDecodingFunction;
            fprintf(['Loaded speed and classification information from ' allLast '\n']);
        else
            displayText(mainWindow,stim.instruct_summary,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        end
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);      %make sure the instructions are always to continue
        
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            DrawFormattedText(mainWindow,'Waiting for scanner start, hold tight!','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            Screen('Flip', mainWindow);
        end
        
        % initialize mot test
        dotEK = initEasyKeys([exp_string_long '_DOT'], SUBJ_NAME,ppt_dir, ...
            'default_respmap', target_scale, ...
            'stimmap', stimmap, ...
            'condmap', condmap, ...
            'prompt_dur', probe_promptDur, ...
            'listen_dur', probe_listenDur, ...
            'exp_onset', stim.subjStartTime, ...
            'trigger_next', false, ...
            'device', DEVICE);
        subjectiveEK = initEasyKeys([exp_string_long '_SUB'], SUBJ_NAME,ppt_dir, ...
            'default_respmap', subj_scale, ...
            'stimmap', stimmap, ...
            'condmap', condmap, ...
            'prompt_dur', 0.001, ...
            'exp_onset', stim.subjStartTime, ...
            'trigger_next', false, ...
            'device', DEVICE);
        [dotEK,subjectiveEK] = startSession(dotEK,subjectiveEK);
        
        digits_scale = makeMap({'even','odd'},[0 1],keyCell([1 5]));
        condmap = makeMap({'even','odd'});
        digitsEK = initEasyKeys('odd_even', SUBJ_NAME,ppt_dir, ...
            'default_respmap', digits_scale, ...
            'condmap', condmap, ...
            'trigger_next', false, ...
            'prompt_dur', digits_promptDur, ...
            'device', DEVICE);
        
        % present a fixation for the duration of stim.fixBlock (20 s)
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            [timing.trig.wait timing.trig.waitSuccess] = WaitTRPulse(TRIGGER_keycode,DEVICE);
            runStart = timing.trig.wait;
            displayText(mainWindow,STILLREMINDER,STILLDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            DrawFormattedText(mainWindow,'+','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            Screen('Flip', mainWindow)
            config.wait = stim.fixBlock; % stim.fixation - 20 s
        else
            runStart = GetSecs;
            timing.trig.wait = runStart; %for debugging purposes
            config.wait = 0;
        end
        
        config.TR = stim.TRlength;
        config.nTRs.ISI = stim.isiDuration/stim.TRlength;
        %config.nPrompts = 4;
        %config.nTRs.betweenPrompt = stim.betweenPrompt/stim.TRlength;
        config.nTRs.target = stim.targ_dur/stim.TRlength;
        %config.nTRs.prompt = vis_promptDur/stim.TRlength;
        %config.nTRs.beforePrompt = stim.beforePrompt/stim.TRlength;
        config.nTRs.motion = stim.trialDur/stim.TRlength;
        config.nTRs.probe = probe_promptDur/stim.TRlength;
        config.nTRs.feedback = stim.fbDuration/stim.TRlength;
        config.nTRs.math = (num_digit_qs*(digits_promptDur + digits_isi))/stim.TRlength;
        config.nTrials = length(stim.stim);
        
        config.nTRs.perTrial = (config.nTRs.ISI + config.nTRs.target + config.nTRs.motion ...
            + config.nTRs.probe + config.nTRs.feedback + config.nTRs.math);
        config.nTRs.perBlock = config.wait/stim.TRlength + (config.nTRs.perTrial)*config.nTrials+ config.nTRs.ISI; %includes the last ISI and 20 s fixation in the beginning but ...
        % does NOT include that 10 s at the end (last TRs)
        
        % calculate all future onsets
        timing.plannedOnsets.preITI(1:config.nTrials) = runStart + config.wait + ((0:config.nTrials-1)*config.nTRs.perTrial)*config.TR;
        timing.plannedOnsets.target(1:config.nTrials) = timing.plannedOnsets.preITI + config.nTRs.ISI*config.TR;
        %timing.plannedOnsets.motionStart(1:config.nTrials) = timing.plannedOnsets.target + config.nTRs.target*config.TR;
        for mTr = 1:config.nTRs.motion
            if mTr == 1
                timing.plannedOnsets.motion(mTr,1:config.nTrials) = timing.plannedOnsets.target + config.nTRs.target*config.TR;
            else
                timing.plannedOnsets.motion(mTr,1:config.nTrials) = timing.plannedOnsets.motion(mTr-1,:) + config.TR;
            end
        end
        timing.plannedOnsets.probe(1:config.nTrials) = timing.plannedOnsets.motion(1,:) + config.nTRs.motion*config.TR;
        timing.plannedOnsets.feedback(1:config.nTrials) = timing.plannedOnsets.probe + config.nTRs.probe*config.TR;
        timing.plannedOnsets.math(1:config.nTrials) = timing.plannedOnsets.feedback + config.nTRs.feedback*config.TR;
        timing.plannedOnsets.lastITI = timing.plannedOnsets.math(end) + config.nTRs.math*config.TR;%%make sure it pauses for this one
        
        allMotionTRs = convertTR(runStart,timing.plannedOnsets.motion,config.TR); %row,col = mTR,trialnumber
        addTR = 0;
        %showFiles = 1;
        if ~CURRENTLY_ONLINE %we have to add the 10 TR into back to go with prev data
            allMotionTRs = allMotionTRs + stim.fixBlock/stim.TRlength;
            addTR = stim.fixBlock/stim.TRlength;
        end
        rtData.classOutputFileLoad = nan(1,config.nTRs.perBlock + addTR);
        rtData.classOutputFile = cell(1,config.nTRs.perBlock + addTR);
        rtData.rtDecoding = nan(1,config.nTRs.perBlock+ addTR);
        %rtData.smoothRTDecoding = nan(1,config.nTRs.perBlock+ addTR);
        rtData.rtDecodingFunction = nan(1,config.nTRs.perBlock+ addTR);
        rtData.fileList = cell(1,config.nTRs.perBlock + addTR);
        rtData.newestFile = cell(1,config.nTRs.perBlock + addTR);
        % repeat
        stim.lastSpeed = nan(1,stim.num_realtime);%going to save it in a matrix of run,stimID
        stim.lastRTDecoding = nan(1,stim.num_realtime); %file 9 that's applied now
        stim.lastRTDecodingFunction = nan(1,stim.num_realtime);
        stim.changeSpeed = nan(mTr,length(stim.cond));
        stim.motionSpeed = nan(mTr,length(stim.cond));
        
        
        %save the timing, stim ID, and stim conditions here!
        stimID = stim.id;
        stimCond = stim.cond;
        sessionInfoFile = fullfile(ppt_dir, ['SessionInfo' '_' num2str(SESSION) '.mat']);
        save(sessionInfoFile, 'stimCond','stimID', 'timing', 'config'); 
                
        rtData.foundFn = [];
        rtData.RTVEC = {};
        rtData.allSpeeds = {};
        rtData.allTimeChanges = {};
        fileNumber = 1;
        for n=1:length(stim.cond)
            stim.trial = n;
            train = [];
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.preITI(n), timing.trig.preITI_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.preITI(n));
            end
            timespec = timing.plannedOnsets.preITI(n) - SLACK;
            timing.actualOnsets.preITI(n) = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.preITI(n) - timing.plannedOnsets.preITI(n));
            
            % reset quest to prevent us from getting stuck if we had an error in the first few trials
            if stair
                suggestion = QuestQuantile(questStruct);
                if stim.trial <= 4 && (suggestion > stim.maxspeed-8) %this restarts it if they get first 3 wrong
                    questStruct = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range); %Quest function for face blocks
                    suggestion = QuestQuantile(questStruct);
                end
                % initialize trial
                blend = min(suggestion,stim.maxspeed-0.5);% 30-0.5 is as slow as we're willing to go
                blend = max(blend, 0.5); %29.5 is as fast as we're willing to go
                stim.speed(stim.trial) = stim.maxspeed - blend;
            end
            actualspeed = max([stim.speed(stim.trial) 0.3]);
            repulsor_force(stim.trial) = repulsor_force_small * actualspeed;
            %fprintf('For trial %i: speed = %.2d\n', stim.trial,stim.speed(stim.trial));
            dotTarg = []; dotPos = []; dotTargPos = [];
            
            cue = stim.stim{stim.trial};
            embedded_keys = keyCell;
            embedded_scale = subj_scale;
            if stim.cond > 3 % lure condition
                embedded_cresp = keyCell(1);
            else
                embedded_cresp = keyCell(3:5);
            end
            
            % initialize KbQueue
            [embedded_keys, valid_keycode, embedded_cresp] = keyCheck(embedded_keys,embedded_cresp);
            KbQueueCreate(DEVICE,valid_keycode);
            KbQueueStart(DEVICE);
            % initialize dots
            [dots phantom_dots] = initialize_dots(num_dots,stim.num_targets,stim.square_dims,stim.dot_diameter);
            
            % choose a probe
            is_new = logical(randi(2)-1); %either a 1 or 0, if 1; 1 = Z no, 0 = Y
            if is_new, stim.cresp(stim.trial) = 1; else stim.cresp(stim.trial) = 2; end
            if is_new || ~stim.num_targets
                dots(stim.num_targets + randi(num_dots - stim.num_targets)).is_probe = true; %makes it so the first dot can't be the probe
            else dots(randi(stim.num_targets)).is_probe = true; %this forces the first dot index to be the probe, 50% chances
            end
            
            % reveal targets
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.target(n), timing.trig.target_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.target(n));
            end
            timespec = timing.plannedOnsets.target(n) - SLACK;
            show_targs = true; show_probe = false; prompt_active = false;
            targetBin = [];
            shade = 0;
            initFeedback = config.initFeedback;
            initFunction = config.initFunction;
            current_speed = stim.speed(stim.trial);

            if SESSION > MOT{1}
                %load the previous last speed for this stimulus and set
                %this as the speed
                initSpeed = lastSpeed(stim.id(stim.trial)); %check the indexing is right!!
                current_speed = initSpeed;
                initFeedback = lastDecoding(stim.id(stim.trial));
                initFunction = lastDecodingFunction(stim.id(stim.trial));
            end
            rtData.allSpeeds{n} = [];
            rtData.allSpeeds{n}(1) = current_speed;
            rtData.allTimeChanges{n} = [];
            if current_speed < 0
                shade = abs(current_speed);
            end
            [~,timing.actualOnsets.target(n)] = dot_refresh(mainWindow,targetBin,dots,shade,square_bounds,stim.dot_diameter,COLORS,cue,show_targs,show_probe,prompt_active,timespec); %flips but not moving yet
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.target(n) - timing.plannedOnsets.target(n));
            
            % initialize dot following
            show_targs = false;
            stim.frame_counter(stim.trial) = 0;
            
            TRcounter = 0;
            prompt_counter = 0;
            fps = 30;
            repulse = stim.repulse;
            remainingTR = ones(mTr,1);
            
            prompt_active = 0;
            fs = [];
            
            printTR = ones(1,mTr);
            
            fileTR = 1;
            waitForPulse = false;
            printlog(LOG_NAME,'trial\tTR\tprompt active\tspeed\tds\tflip error\tFound File #\tThisTR\tCategSep\n');
            check = [];
            rtData.RTVEC{n} = [];
           
            while abs(GetSecs - timing.plannedOnsets.probe(n)) > 0.050;%SLACK*2 %so this is just constatnly running, stops when it's within a flip
               
                stim.frame_counter(stim.trial) = stim.frame_counter(stim.trial) + 1;
                % here we go!
                
                % this whole thing will happen ONLY at the start of every
                % TR
                nextTRPos = find(remainingTR,1,'first');
                if ~isempty(nextTRPos)
                    nextTRTime = timing.plannedOnsets.motion(nextTRPos,stim.trial);
                    if abs(GetSecs - nextTRTime) <= 0.050
                        %look for speed update here
                        TRcounter = nextTRPos; %update TR count (initialized at 0): so it's the TR that we're currently ON
                        thisTR = allMotionTRs(TRcounter,n);
                        waitForPulse = true;
                        if ismember(TRcounter,promptTRs)
                            prompt_active = 1;
                            KbQueueFlush(DEVICE);
                            prompt_counter = prompt_counter + 1;
                            train.onset(prompt_counter) = GetSecs;
                            check(prompt_counter) = 1;
                        elseif ismember(TRcounter-1,promptTRs) && check(prompt_counter)
                            [train.acc(prompt_counter), train.resp{prompt_counter}, ~, train.rt(prompt_counter), ~, train.resp_str{prompt_counter}] = ...
                                multiChoice(queueCheck, embedded_keys, embedded_scale, embedded_cresp, GetSecs, DEVICE, [],sum(keys.map(3:5,:)),subj_map);
                            
                            if isempty(train.resp{prompt_counter}), train.resp{prompt_counter} = nan; end % timeout
                            if isnan(train.resp{prompt_counter}),
                                train.resp_str{prompt_counter} = nan;
                            else train.resp_str{prompt_counter} = embedded_keys{train.resp{prompt_counter}};
                            end
                            subjectiveEK = easyKeys(subjectiveEK, ...
                                'stim', stim.stim{stim.trial}, ...
                                'onset', train.onset(prompt_counter), ...
                                'cond', stim.cond(stim.trial), ...
                                'nesting', [SESSION stim.trial prompt_counter], ...
                                'cresp', embedded_cresp, ...
                                'simulated_key', train.resp_str{prompt_counter}, ...
                                'cresp_map', sum(keys.map(3:5,:)), 'valid_map', subj_map);
                            check(prompt_counter) = false;
                            fprintf('The response for prompt %i was %i\n', prompt_counter, train.resp{prompt_counter})
                            prompt_active = false;
                        end
                        % the realtime part could help because you know it
                        % only goes into that loop once
                        remainingTR(nextTRPos) = 0;
                        stim.motionSpeed(TRcounter,n) = current_speed; %this is the speed at the onset of the TR
                    elseif GetSecs - nextTRTime > 0.050
                        % then we actually missed the current TR
                        % counter--look for the next
                        % check that we didn't just miss the trigger
                        diffT = (GetSecs - timing.plannedOnsets.motion(:,stim.trial));
                        %[~,correctTR] = min(abs(GetSecs - timing.plannedOnsets.motion(:,stim.trial)));
                        % this will go towards the last TR that we passed
                        correctTR = find(floor(diffT)==0);
                        TRcounter = correctTR;
                        thisTR = allMotionTRs(TRcounter,n);
                        remainingTR(1:correctTR) = 0;
                    end
                end
                
              
   
                waitFrames = 15; % looking at half second rate--start to begin % ACM changed this from 15 to 2 after subject 2
                % but we don't want this to look for the same file more
                % than once or go through the same loop
                if realtime
                    if TRcounter >=4
                        thisTR = allMotionTRs(TRcounter,n); %this is the TR we're actually on KEEP THIS WAY--starts on 4, ends on 10
                        fileTR = thisTR - 1; 
                        if scanNum
                            if ~mod(stim.frame_counter(n),waitFrames) %check only evert so often to spare memory
                                % first check if there's a new file
                                rtData.fileList{thisTR} = ls(classOutputDir);
                                allFn = dir([classOutputDir 'vol' '*']);
                                dates = [allFn.datenum];
                                names = {allFn.name};
                                [~,newestIndex] = max(dates);
                                rtData.newestFile{thisTR} = names{newestIndex};
                                % figure out the NUMBER file that it came from
                                filename = rtData.newestFile{thisTR};
                                startI = 5;
                                fileNumber = str2double(filename(startI:(length(filename)-4)));
                                % let's see if we've loaded this before
                                if ~ismember(fileNumber,rtData.foundFn) && ismember(fileNumber,allMotionTRs(3:end,n)) % so only accept it if it's one of the TR's for that trial
                                    % now add this to the found file
                                    rtData.foundFn(end+1) = fileNumber;

                                    rtData.classOutputFileLoad(fileNumber) = 1;
                                    tempStruct = load(fullfile(classOutputDir,filename));
                                    rtData.rtDecoding(fileNumber) = tempStruct.classOutput;
                                    rtData.RTVEC{n}(end+1) = rtData.rtDecoding(fileNumber);
                                    rtData.rtDecodingFunction(fileNumber) = PID(rtData.RTVEC{n},Kp,Ki,Kd,OptimalForget,maxIncrement);
                                    
                                    % update speed here
                                    current_speed = current_speed + rtData.rtDecodingFunction(fileNumber); % apply in THIS TR what was from 2 TR's ago (indexed by what file it is) so file 3 will be applied at TR5!
                                    stim.changeSpeed(TRcounter,n) = rtData.rtDecodingFunction(fileNumber); %speed changed ON that TR
                                    % so this will update every time a new
                                    % file comes through!
                                    % you can check rtData.rtDecoding to
                                    % see which files are caught here
                                    
                                    
                                    % make sure speed is between [stim.minspeed
                                    % stim.maxspeed] (0,30) right now
                                    current_speed = min([stim.maxspeed current_speed]);
                                    current_speed = max([stim.minspeed current_speed]);
                                    rtData.allSpeeds{n}(end+1) = current_speed;
                                    rtData.allTimeChanges{n}(end+1) = GetSecs;
                                end
                            end
                            
                        else
                            rtData.rtDecoding(fileTR) = rand(1)*2-1;
                            rtData.RTVEC{n}(end+1) = rtData.rtDecoding(fileTR);
                            rtData.rtDecodingFunction(fileTR) = PID(rtData.RTVEC{n},Kp,Ki,Kd,OptimalForget,maxIncrement);
                            % update speed here
                            current_speed = current_speed + rtData.rtDecodingFunction(fileNumber); % apply in THIS TR what was from 2 TR's ago (indexed by what file it is) so file 3 will be applied at TR5!
                            stim.changeSpeed(TRcounter,n) = rtData.rtDecodingFunction(fileNumber); %speed changed ON that TR
                            
                            % make sure speed is between [stim.minspeed
                            % stim.maxspeed] (0,30) right now
                            current_speed = min([stim.maxspeed current_speed]);
                            current_speed = max([stim.minspeed current_speed]);
                        end
                    end
                end
                
                if TRcounter %only do this once the motion should begin
                    [dots, fs,shade] = dot_compute(dots,current_speed,shade,stim.square_dims,stim.dot_diameter,phantom_dots,WINDOWSIZE,repulse,fs,repulsor_force(n));
                    if waitForPulse
                        if CURRENTLY_ONLINE && SESSION >TOCRITERION3 %localizer and up
                            [timing.trig.motion(TRcounter,n), timing.trig.motion_Success(TRcounter,n)] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.motion(TRcounter,n)); %minus 1 because first TR is motionStart
                        end
                        timespec = nextTRTime - SLACK;
                        [targetBin, timing.actualOnsets.motion(TRcounter,n)] = dot_refresh(mainWindow,targetBin,dots,shade,square_bounds,stim.dot_diameter,COLORS,cue,show_targs,show_probe,prompt_active,timespec);
                        %fprintf('Flip time error = %.4f\n', timing.actualOnsets.motion(motion_counter,stim.trial) - timing.plannedOnsets.motion(motion_counter,stim.trial));
                        %fprintf(['mTR ' num2str(motion_counter) '; Speed = ' num2str(current_speed) '\n']);
                        waitForPulse = false;
                        %display report at the end of the TR
                    else
                        [targetBin] = dot_refresh(mainWindow,targetBin,dots,shade,square_bounds,stim.dot_diameter,COLORS,cue,show_targs,show_probe,prompt_active);
                    end
                end
                if TRcounter > 1 && (GetSecs >= timing.plannedOnsets.motion(TRcounter,n) + config.TR-.1) && printTR(TRcounter) %after when should have found file
                    %z = GetSecs - timing.plannedOnsets.motion(TRcounter,n);
                    printlog(LOG_NAME,'%d\t%d\t%d\t\t%5.3f\t%5.3f\t%5.4f\t\t%i\t\t%d\t\t%5.3f\n',n,TRcounter,prompt_active,current_speed,rtData.rtDecodingFunction(fileNumber),timing.actualOnsets.motion(TRcounter,stim.trial) - timing.plannedOnsets.motion(TRcounter,stim.trial),fileNumber,thisTR,rtData.rtDecoding(fileNumber));
                    printTR(TRcounter) = 0;
                elseif TRcounter ==1 && (GetSecs >= timing.plannedOnsets.motion(TRcounter,n) + config.TR-.1) && printTR(TRcounter)
                    printlog(LOG_NAME,'%d\t%d\t%d\t\t%5.3f\t%5.3f\t%5.4f\t\t%i\t\t%d\t\t%5.3f\n',n,TRcounter,prompt_active,current_speed,rtData.rtDecodingFunction(fileNumber),timing.actualOnsets.motion(TRcounter,stim.trial) - timing.plannedOnsets.motion(TRcounter,stim.trial),fileNumber,thisTR,rtData.rtDecoding(fileNumber));
                    printTR(TRcounter) = 0;
                end

                
            end  %30 s trial ends here THEN probe
            % update parameters at the end of the trial
            if realtime
                stim.lastSpeed(stim.id(stim.trial)) = current_speed; %going to save it in a matrix of run,stimID
                stim.lastRTDecoding(stim.id(stim.trial)) = rtData.rtDecoding(fileNumber); %file 9 that's applied now
                stim.lastRTDecodingFunction(stim.id(stim.trial)) = rtData.rtDecodingFunction(fileNumber);
            end
            % present targetness probe
            KbQueueRelease;
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.probe(n), timing.trig.probe_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.probe(n));
            end
            show_probe = true; prompt_active = false;
            timespec = timing.plannedOnsets.probe(n) - SLACK;
            [~,timing.actualOnsets.probe(n)] = dot_refresh(mainWindow,[],dots,shade,square_bounds,stim.dot_diameter,COLORS,cue,show_targs,show_probe,prompt_active,timespec);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.probe(n) - timing.plannedOnsets.probe(n));
            dotEK = easyKeys(dotEK, ...
                'nesting', [SESSION stim.trial], ...
                'stim', stim.stim{stim.trial}, ...
                'cond', stim.cond(stim.trial), ...
                'onset', timing.actualOnsets.probe(n), ...
                'cresp', stim.keys(stim.cresp(stim.trial)), ...
                'cresp_map', dot_map(stim.cresp(stim.trial),:), 'valid_map', target_map );
            
            % log trial
            stim.dur(stim.trial) = GetSecs - timing.actualOnsets.target(n);
            stim.dotLog{n} = dots;
            % stim.vis_train{n} = train;
            clear dots
            
            
            % feedback
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.feedback(n), timing.trig.feedback_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.feedback(n));
            end
            Screen('TextSize',mainWindow,stim.smallSize);
            Screen('TextFont', mainWindow,'Arial');
            timespec = timing.plannedOnsets.feedback(n) - SLACK;
            if dotEK.trials.acc(end)
                timing.actualOnsets.feedback(n) = displayText_specific(mainWindow,'!!!','center',COLORS.GREEN,WRAPCHARS,timespec);
            else
                timing.actualOnsets.feedback(n) = displayText_specific(mainWindow,'X','center',COLORS.RED,WRAPCHARS,timespec);
            end
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.feedback(n) - timing.plannedOnsets.feedback(n));
            
            % PRACTICE FEEDBACK
            
            if SESSION == MOT_PRACTICE && (stim.trial <= 3)
                pause(stim.fbDuration)
                onFB = GetSecs;
                Screen('TextSize',mainWindow,stim.smallSize);
                % mot performance (top priority)
                if dotEK.trials.acc(end)
                    displayText(mainWindow,['Great work! You successfully tracked the target dots and responded correctly.\n\n' ...
                        'Did you remember to keep your eyes fixed on the word in the middle?\n\nAlso, did you remember to move ' ...
                        'your eyes to the white dot at the end of the trial?\n\n' ...
                        '-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                elseif isnan(dotEK.trials.resp(end))
                    displayText(mainWindow,['Oops, you didn''t respond to the dot-following question in time.\n\n' ...
                        'Press PINKY when the white dot is a target and THUMB ' ...
                        'when the white dot is not a target.(Think YES TARGET = PINKY, NO NOT TARGET = THUMB.)\n\nAlso, remember that although ' ...
                        'you have other jobs to do in this multi-tasking task, keeping track of the target dots is the ' ...
                        'most important thing you have to do.\n\n' ...
                        '-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                else
                    displayText(mainWindow,['Oops, you guessed incorrectly about whether the question dot was a target.\n\n' ...
                        'Press PINKY when the white dot is a target and THUMB ' ...
                        'when the white dot is not a target. (Think YES TARGET = PINKY, NO NOT TARGET = THUMB.)\n\nAlso, remember that although ' ...
                        'you have other jobs to do in this multi-tasking task, keeping track of the target dots is the ' ...
                        'most important thing you have to do.\n\n' ...
                        '-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                end
                waitForKeyboard(kbTrig_keycode,DEVICE);
                
                
                % subjective responding
                displayText(mainWindow,['Also, you responded to ' num2str(round(100*(length(train.resp) - sum(cellfun(@isnan,train.resp)))/100)) ' out of ' ...
                    num2str(length(train.resp)) ' "rate detail" prompts (i.e., when the fixation circle turns red). Even when you ' ...
                    'have no free mental resources for mental imagery, you should still always try to respond to these ' ...
                    'prompts. Also, remember that your rating should reflect your CURRENT clarity, even if ' ...
                    'that is not as good as it CAN be.\n\n Just to be sure you remember the scale, we''ll give you another ' ...
                    'look at the hand map.\n\n' ...
                    '-- Press ' PROGRESS_TEXT ' to continue, then press it again to clear the hand map --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                waitForKeyboard(kbTrig_keycode,DEVICE);
                Screen('DrawTexture',mainWindow,keymap_prompt,[],[],[]); %[0 0 keymap_dims],[topLeft topLeft+keymap_dims]);
                Screen('Flip',mainWindow);
                waitForKeyboard(kbTrig_keycode,DEVICE);
                
                
                % prepare for odd-even
                displayText(mainWindow,['Now get ready for the rapid odd/even calculations: indicate "even" sums with your THUMB and "odd" sums with your PINKY finger.\n\n' ...
                    '-- Press ' PROGRESS_TEXT ' to continue --'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                waitForKeyboard(kbTrig_keycode,DEVICE);
                offFB = GetSecs;
                timing.plannedOnsets.preITI(n+1:end) = timing.plannedOnsets.preITI(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.target(n+1:end) = timing.plannedOnsets.target(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.motion(:,n+1:end) = timing.plannedOnsets.motion(:,n+1:end) + (offFB - onFB);
                %timing.plannedOnsets.prompt(:,n+1:end) = timing.plannedOnsets.prompt(:,n+1:end) + (offFB - onFB);
                timing.plannedOnsets.probe(n+1:end) = timing.plannedOnsets.probe(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.feedback(n+1:end)= timing.plannedOnsets.feedback(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.math(n:end) = timing.plannedOnsets.math(n:end) + (offFB - onFB);
                timing.plannedOnsets.lastITI = timing.plannedOnsets.lastITI + (offFB - onFB);
            end
            
            
            %math
            if CURRENTLY_ONLINE && SESSION > TOCRITERION3
                [timing.trig.math(n), timing.trig.math_Success(n)] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.math(n));            end
            timespec = timing.plannedOnsets.math(n) - SLACK;
            [stim.digitAcc(stim.trial), stim.digitRT(stim.trial) timing.actualOnsets.math(n)] = odd_even(digitsEK,num_digit_qs,digits_promptDur,digits_isi,minimal_format,mainWindow,keyCell([1 5]),COLORS,DEVICE,SUBJ_NAME,[SESSION stim.trial], SLACK, timespec, keys);
            fprintf('Flip time error = %.4f\n', timing.actualOnsets.math(n) - timing.plannedOnsets.math(n));
            
            % MOT prep feedback: fight drowsiness
            if SESSION == MOT_PREP && ~mod(n,5)
                pause(digits_isi)
                onFB = GetSecs;
                displayText(mainWindow,['BREAK\n\nYou''re doing great! Feel free to use this opportunity to get up and stretch if you need to. '...
                    'Before continuing, we want to make sure that you''re focused because if you guess randomly '...
                    'you''ll throw off our algorithm and we may not be able to continue with the study.\n Instead, if you''re feeling drowsy, '...
                    'please take a short rest now and resume once you''re ready. \n\n -- Press ' PROGRESS_TEXT ' when you''re ready to continue--'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                waitForKeyboard(kbTrig_keycode,DEVICE);
                
                displayText(mainWindow,['Trial ' num2str(stim.trial) ' of ' num2str(length(stim.stim)) ' complete.'],minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                pause(0.05);
                
                offFB = GetSecs;
                timing.plannedOnsets.preITI(n+1:end) = timing.plannedOnsets.preITI(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.target(n+1:end) = timing.plannedOnsets.target(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.motion(:,n+1:end) = timing.plannedOnsets.motion(:,n+1:end) + (offFB - onFB);
                %  timing.plannedOnsets.prompt(:,n+1:end) = timing.plannedOnsets.prompt(:,n+1:end) + (offFB - onFB);
                timing.plannedOnsets.probe(n+1:end) = timing.plannedOnsets.probe(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.feedback(n+1:end)= timing.plannedOnsets.feedback(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.math(n+1:end) = timing.plannedOnsets.math(n+1:end) + (offFB - onFB);
                timing.plannedOnsets.lastITI = timing.plannedOnsets.lastITI + (offFB - onFB);
            end
            
            % report
            stim.expDuration = (GetSecs - stim.subjStartTime) / 60; % experiment time in mins
            if ~isnan(dotEK.trials.resp(end)) &&  stair
                % update questStruct only if we got a response
                questStruct=QuestUpdate(questStruct,blend,dotEK.trials.acc(end));
            end
            if stair
                stim.tGuess(stim.trial) = QuestMean(questStruct);
                stim.tGuess_sd(stim.trial) = QuestSd(questStruct);
            end
            
            clear train;
            save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        end
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            [timing.trig.lastITI, timing.trig.lastITI_Success] = WaitTRPulse(TRIGGER_keycode,DEVICE, timing.plannedOnsets.lastITI);
        end
        timespec = timing.plannedOnsets.lastITI - SLACK;
        timing.actualOnsets.lastITI = isi_specific(mainWindow,COLORS.MAINFONTCOLOR,timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.lastITI - timing.plannedOnsets.lastITI);
        if GetSecs - timing.actualOnsets.lastITI < 2
            WaitSecs(2 - (GetSecs - timing.actualOnsets.lastITI));
        end
        
        %wait in scanner at end of run
        if CURRENTLY_ONLINE && SESSION > TOCRITERION3
            WaitSecs(10);
        end

        if SESSION < MOT{1}
            save(MATLAB_SAVE_FILE,'stim', 'timing', 'config');
        else
            save(MATLAB_SAVE_FILE, 'stim', 'timing', 'config', 'rtData');
        end
        
        % wrap up
        if SESSION == MOT_PREP
            
            %displayText(mainWindow,CONGRATS,CONGRATSDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            endSession(dotEK, subjectiveEK, CONGRATS);
            load(MATLAB_SAVE_FILE);
            %subplot(1,2,1)
            figure;
            plot(stim.maxspeed-stim.tGuess);
            %subplot(1,2,2)
            %plot(stim.avg_vis_resp);
            sca
        else
            %displayText(mainWindow,CONGRATS,CONGRATSDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
            endSession(dotEK, subjectiveEK,CONGRATS);
            if SESSION < MOT_LOCALIZER
                mot_realtime05(SUBJECT,SESSION+1,SET_SPEED,scanNum,scanNow);
            end
        end
        sca;
        
        %% FRUIT HARVEST
    case {RSVP,RSVP2}
        
        % stimulus presentation parameters
        secs_per_item = 8*SPEED; % secs per item
        stim.targetLatencyMean = 12*SPEED; % time between target exposures for use in distribution. count on this being about 0.5s longer than specified, since filler items will spill over
        stim.targetLatencySd = 0*SPEED; % no jitter
        stim.shortest_expos = 0.300*SPEED;
        stim.longest_expos = 0.750*SPEED;
        stim.nminus1 = 0.5*SPEED; % trial before cue presentation
        stim.fruit_extradelay = 0;
        stim.detectKey = {INDEXFINGER};
        stim.isiDuration = 2*SPEED;
        % session-based declarations
        instruct = ['THE GREAT FRUIT HARVEST\n\nIn this task, words will flash up on the screen very quickly, one after another. If you notice a word that is a ' ...
            'type of fruit, please press the INDEX finger.\n\nNote: there are very few fruit, so make sure to catch them!\n\n-- Press INDEX to begin --'];
        stim.condDuration(REALTIME) = 0; stim.condDuration(OMIT)= 0;
        exposure = 1;
        
        % load or initialize exclusion words and cue durations
        relative_paths = 1; % for use with file search functions
        stim.exclusionList = [preparedCues];
        
        % during practice, filler is going to be a small set of repeated words to familiarize
        if SESSION == RSVP
            stim.fillerCues = lureWords(1:7);
        else stim.fillerCues = lureWords(8:23);
        end
        stim.num_short = 0; stim.num_long = 0; stim.num_omit = 0;
        PROGRESS = INDEXFINGER;
        % final initialization
        fillerCueTargets = readStimulusFile(CUETARGETFILE,ALLMATERIALS);
        triggerNext = false;
        condmap = makeMap({'realtime','omit','lure','fruit'});
        fruitHarvestEK = initEasyKeys([exp_string_long '_FH'], SUBJ_NAME, ppt_dir,...
            'default_respmap', rsvp_scale, ...
            'condmap', condmap, ...
            'trigger_next', triggerNext, ...
            'device', DEVICE);
        fruitHarvestEK = startSession(fruitHarvestEK);
        stim.trial = 0;
        countdown = 0;
        
        stim.availableFruit = round(length(stim.fillerCues) * 0.5);
        stim.scanLength = secs_per_item * length(stim.fillerCues) + stim.availableFruit;
        
        % display or skip instructions, depending on session
        displayText(mainWindow,instruct,minimumDisplay,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        stim.subjStartTime = waitForKeyboard(kbTrig_keycode,DEVICE);
        
        stim.sessionStartTime = GetSecs();
        %[~,ision] = isi(mainWindow,stim.TRlength,COLORS.MAINFONTCOLOR);
        isi(mainWindow,stim.isiDuration,COLORS.MAINFONTCOLOR); %cange this to be actual ISI
        lastCue = GetSecs() - 4; %initialize so it will present maybe 8 s later
        idealLag = normrnd(stim.targetLatencyMean,stim.targetLatencySd);
        stim.enterLoop = GetSecs;
        stim.expDuration = 0;
        % keep pumping out filler words until the time is up
        while stim.expDuration < ((stim.scanLength)/60) %4 seconds less
            
            % initialize trial
            stim.trial = stim.trial + 1;
            cueDistance = GetSecs() - lastCue; % rear-view mirror
            if countdown, countdown = countdown - 1; end
            
            % figure out if there are any cues left
            if ~stim.availableFruit %if there are no more fruit
                timeToCue = inf;
            else % if there are cues left, check whether it's nearly time to present one
                timeToCue =  idealLag - cueDistance;
                buffer_time = unifrnd(stim.shortest_expos,stim.longest_expos);
                if ~countdown && timeToCue < stim.nminus1+buffer_time
                    countdown = 2; % step 5 is sync to TR; step 4 is 0.5s constant filler; step 3 is cue; then two fillers and back to normal
                end
            end
            
            % get stimulus, make sure it wasn't just used
            stim.stim{stim.trial} = [];
            while isempty(stim.stim{stim.trial})
                candidate = [];
                while isempty(candidate)
                    candidate = stim.fillerCues{randi(length(stim.fillerCues))};
                    if stim.trial > 1
                        if strcmp(stim.stim{stim.trial-1},candidate)
                            candidate = [];
                        end
                    end
                end
                stim.stim{stim.trial} = candidate;
            end
            % by default, words are filler with random duration
            stim.cond(stim.trial) = LUREWORD;
            cresp = {nan};
            cresp_map = zeros(1,256);
            valid_map = keys.map(2,:);
            stim.promptDur(stim.trial) = unifrnd(stim.shortest_expos,stim.longest_expos);
            duration = stim.promptDur(stim.trial);
            
            % now override details based on real item type
            switch countdown
                case 1 % target--if the countdown is 1, make it a fruit with 1 s duration (and 33% of the time)
                    % what condition will the cue come from?
                    if stim.availableFruit && (rand() < 0.33) %1/3 of the time take a fruit (stim.availableFruit / (stim.num_short + stim.num_long + stim.availableFruit - cueIndex(EASY) - cueIndex(HARD))))
                        cresp = {INDEXFINGER};
                        cresp_map = keys.map(2,:);
                        stim.cond(stim.trial) = FRUIT;
                        stim.stim{stim.trial} = fillerCueTargets{randi(length(fillerCueTargets))};
                        stim.promptDur(stim.trial) = 1;
                        lastCue = GetSecs() + duration + stim.fruit_extradelay; % assures us a reasonable buffer before memory cue is presented
                        stim.availableFruit = stim.availableFruit - 1;
                        idealLag = normrnd(stim.targetLatencyMean,stim.targetLatencySd);
                        
                    end
                case 2 % filler with n-1 duration
                    stim.promptDur(stim.trial) = stim.nminus1;
            end
            stim.condString{stim.trial} = CONDSTRINGS{stim.cond(stim.trial)};
            
            trial_message = ['sess' num2str(SESSION) '_trial' num2str(stim.trial) '_cond' num2str(stim.cond(stim.trial))];
            
            
            % present stimulus
            if  stim.promptDur(stim.trial) > 0.025 % don't bother drawing less than 25ms--this is the only time we'll draw
                DrawFormattedText(mainWindow,stim.stim{stim.trial},'center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
                if stim.trial-1 > 0 %if it's not the first trial, check to see when th elast trial happened to determine best flip time
                    if stim.promptDur(stim.trial-1)>0.025 %if the last trial occured
                        timing.plannedOnsets.cue(stim.trial) = timing.actualOnsets.cue(stim.trial-1) + stim.promptDur(stim.trial-1);
                        timespec = timing.plannedOnsets.cue(stim.trial)-SLACK; %the onset should be the last onset + the duration
                    elseif promptDur(stim.trial-1)<0.025 && stim.trial-2>0
                        timing.plannedOnsets.cue(stim.trial) = timing.actualOnsets.cue(stim.trial-2) + stim.promptDur(stim.trial-2);
                        timespec = timing.plannedOnsets.cue(stim.trial)-SLACK; %if last trial didn't occur, send it at the last onset
                    else %if it's not the first trial, the last trial didn't happen, and two trials ago didn't happen (unlikely)
                        timing.plannedOnsets.cue(stim.trial) = GetSecs;
                        timespec = timing.plannedOnsets.cue(stim.trial) - SLACK;
                    end
                else %if it's the first trial, flip now
                    timing.plannedOnsets.cue(stim.trial) = GetSecs;
                    timespec = timing.plannedOnsets.cue(stim.trial) - SLACK;
                end
                timing.actualOnsets.cue(stim.trial) = Screen('Flip',mainWindow,timespec);
                fprintf('Flip time error = %.4f\n', timing.actualOnsets.cue(stim.trial) - timing.plannedOnsets.cue(stim.trial));
                
                fruitHarvestEK = easyKeys(fruitHarvestEK, ...
                    'onset', timing.actualOnsets.cue(stim.trial), ...
                    'stim', stim.stim{stim.trial}, ...
                    'cond', stim.cond(stim.trial), ...
                    'cresp', cresp, ...
                    'nesting', [SESSION stim.trial], ...
                    'prompt_dur', stim.promptDur(stim.trial), ...
                    'cresp_map', cresp_map, 'valid_map', valid_map);
                
            else
                fprintf('skip')
            end
            
            % report
            % log the stimulus and compute prior exposures
            stim.expDuration = (GetSecs() - stim.sessionStartTime) / 60;
            stim.harvest_rate = mean(fruitHarvestEK.trials.acc(stim.cond == FRUIT));
            stim.false_fruit = sum(1-fruitHarvestEK.trials.acc(stim.cond ~= FRUIT));
            
            % save to file
            if mod(stim.trial,10)==0 || (stim.enterLoop >= (stim.scanLength-(STABILIZATIONTIME)))
                stim.expDuration = (GetSecs - stim.enterLoop) / 60; % experiment time in mins
                save(MATLAB_SAVE_FILE,'stim','timing');
            end 
        end
        isi(mainWindow,stim.isiDuration,COLORS.MAINFONTCOLOR); %cange this to be actual ISI
        
        % present participant with feedback
        taskSummary = ['Fruit harvest rate: ' num2str(round(stim.harvest_rate*100)) '%\nFalse fruit: ' num2str(stim.false_fruit) ' items'];
        feedbackString = ['Run complete! Your performance in this run:\n\n' taskSummary ]; % Your overall performance:\n\n' overallSummary '\n\n
        displayText(mainWindow,feedbackString,CONGRATSDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        printlog(LOG_NAME,['\n\nSESSION ' int2str(SESSION) ' ended ' datestr(now) ' for SUBJECT number ' int2str(SUBJECT) '\n\n']);
        printlog(LOG_NAME,'\n\n\n******************************************************************************\n');
        
        save(MATLAB_SAVE_FILE,'stim', 'timing');
        
        % wrap up
        if (SESSION == RSVP)
            endSession(fruitHarvestEK, NOTIFY);
            sca
        else
            endSession(fruitHarvestEK, CONGRATS);
            mot_realtime05(SUBJECT,SESSION+1,SET_SPEED,scanNum,scanNow);
        end
        
        %% SCAN PREP
        
    case SCAN_PREP
        % instructions
        displayText(mainWindow,['Great job! Now, we''re now going to take a sequence of short scans before you complete various tasks. Please work through these and we''ll get in ' ...
             'touch with you when you finish.\n\n-- please press the index finger button to continue --'],INSTANT,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        waitForKeyboard(kbTrig_keycode,DEVICE);
        DrawFormattedText(mainWindow,'Waiting for scanner start, hold tight!','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        Screen('Flip', mainWindow);
        
        if CURRENTLY_ONLINE
        [timing.trig.wait timing.trig.waitSuccess] = WaitTRPulse(TRIGGER_keycode,DEVICE);
        runStart = timing.trig.wait;
        displayText(mainWindow,STILLREMINDER,STILLDURATION,'center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        DrawFormattedText(mainWindow,'+','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        Screen('Flip', mainWindow)
        else
            runStart = GetSecs;
        end
       % runStart = timing.trig.wait;
        config.wait = 150; % we want this to be up for 8 seconds to collect sample TR's - this will run for 5 minutes so just stop whenever it's done!
        config.TR = 2;
        timing.plannedOnsets.offEx = runStart + config.wait;
        DrawFormattedText(mainWindow,'Done!','center','center',COLORS.MAINFONTCOLOR,WRAPCHARS);
        timespec = timing.plannedOnsets.offEx-SLACK;
        timing.actualOnsets.offEx = Screen('Flip',mainWindow,timespec);
        fprintf('Flip time error = %.4f\n', timing.actualOnsets.offEx-timing.plannedOnsets.offEx);
        
        save(MATLAB_SAVE_FILE,'timing');
        WaitSecs(2) %wait a little before closing
        sca
        
end
return

%% get image ready for presentation (filename should include path info)
function imageHandle = prepImage(imageFilename,window,scramble)

imageData = imread(imageFilename);
dims = size(imageData);
if exist('scramble','var') && ~isempty(scramble) && scramble %only if you're scrambling images!!
    imageData = reshape(imageData(randperm(numel(imageData))),dims);
end
imageHandle = Screen('MakeTexture', window, imageData);


return

%% quick computation of TR and experiment time relative to experiment start (note: time=0s -> TR1
function [TR timePassed] = calcOnsetTR(trialStartSecs,expStartSecs,trlengthSecs)
timePassed = trialStartSecs-expStartSecs;
TR = ((round((((timePassed)*1000) / (trlengthSecs*1000))*10)/10))+1;
return


function [targetBin,timeon] = dot_refresh(window,targetBin,dots,shade,square_bounds,dot_dia,COLORS,cue,show_targs,show_probe,prompt_active,timespec)
maxdim = 15;
shade = min([maxdim shade]); % make it so max is 30 (or completely black)

% prepare time-gating
target_frame_rate = 30;
time_bins = 0:1000/target_frame_rate:1000;

% define colors for drawing
square_col = COLORS.BLACK;
targ_col = COLORS.RED;
normal_col = COLORS.GREEN;
probe_col = COLORS.WHITE;
font_col = COLORS.WHITE;
centerdot_size = 5;

% determine if we want a picture or black background
if isempty(cue) || ischar(cue)
    picture_mode = false;
    fix_col = COLORS.GREY;
    bumper_size = 0;
else
    picture_mode = true;
    fix_col = COLORS.BLACK;
    bumper_size = 2;
end
if prompt_active
    fix_col = COLORS.RED;
end

% draw square
offset = square_bounds(1:2);
if picture_mode
    bg_texture = Screen('MakeTexture', window, cue);
    Screen('DrawTexture', window, bg_texture, [0 0 size(cue)], square_bounds, []);
else
    Screen('FillRect', window, square_col, square_bounds)
end

% collect dots to draw
drawTargs = []; drawProbe = []; drawNormal = []; drawShadow = [];
trial = find(~isnan(dots(1).pos(:,1)),1,'last');
for i=1:length(dots)
    if picture_mode
        drawShadow = [drawShadow; [round(dots(i).pos(trial,1)+offset(1)) round(dots(i).pos(trial,2)+offset(2))]];
    end
    if show_targs && dots(i).is_target
        drawTargs = [drawTargs; [round(dots(i).pos(trial,1)+offset(1)) round(dots(i).pos(trial,2)+offset(2))]];
    elseif show_probe && dots(i).is_probe
        drawProbe = [drawProbe; [round(dots(i).pos(trial,1)+offset(1)) round(dots(i).pos(trial,2)+offset(2))]];
    else
        drawNormal = [drawNormal; [round(dots(i).pos(trial,1)+offset(1)) round(dots(i).pos(trial,2)+offset(2))]];
    end
end

% draw the dots
%right now going to change for all stages
if ~isempty(drawShadow)
    Screen('FillOval',window,COLORS.BLACK,[drawShadow(:,1:2)-(dot_dia/2)-bumper_size drawShadow(:,1:2)+dot_dia/2 + bumper_size]',dot_dia+ (bumper_size*2))
end
if ~isempty(drawTargs)
    Screen('FillOval',window,targ_col,[drawTargs(:,1:2)-(dot_dia/2) drawTargs(:,1:2) + dot_dia/2]',dot_dia)
end
if ~isempty(drawProbe)
    Screen('FillOval',window,probe_col,[drawProbe(:,1:2)-(dot_dia/2) drawProbe(:,1:2) + dot_dia/2]',dot_dia)
end
if ~isempty(drawNormal)
    if size(drawNormal,1) == length(dots) %search through all dots
        for i = 1:length(dots)
            if ~dots(i).is_target
                col = normal_col;
                col(2) = normal_col(2) - (shade/maxdim)*200;
                Screen('FillOval',window,col,[drawNormal(i,1:2)-(dot_dia/2) drawNormal(i,1:2) + dot_dia/2]',dot_dia)
            else
                Screen('FillOval',window,normal_col,[drawNormal(i,1:2)-(dot_dia/2) drawNormal(i,1:2) + dot_dia/2]',dot_dia)
            end
        end
    else
        %less than that many anyway--fade them
        col = normal_col;
        col(2) = normal_col(2) - (shade/maxdim)*200;
        Screen('FillOval',window,col,[drawNormal(:,1:2)-(dot_dia/2) drawNormal(:,1:2) + dot_dia/2]',dot_dia)
    end
end

% cue or fixation
Screen('TextSize', window,36);
Screen('TextFont', window,'Arial');
if ~isempty(cue) && ~picture_mode && ischar(cue)
    DrawFormattedText(window,cue,'center','center',font_col);
end

% draw fixation point
center = [(square_bounds(1)+square_bounds(3))/2 (square_bounds(2)+square_bounds(4))/2 + 4]';
center_fix_pos = [center - centerdot_size; center + centerdot_size];
Screen('FillOval',window,fix_col,center_fix_pos,dot_dia)

% show the result
Screen('DrawingFinished', window);
if picture_mode, Screen('Close', bg_texture); end

% control dot-drawing refresh rate
if isempty(targetBin)
    now = round(mod(GetSecs(),1)*1000);
    targetBin = histc(now,time_bins) + 1;
else
    bin = 0;
    while bin ~= targetBin
        now = round(mod(GetSecs(),1)*1000);
        bin = histc(now,time_bins);
        if bin ~= targetBin
            pause(0.01);
        end
    end
    targetBin = bin + 1;
end
if ~exist('timespec', 'var')
    timeon = Screen('Flip', window);
else
    timeon = Screen('Flip',window,timespec);
    %    fprintf('Flip time error = %.4f\n', timeon - timespec);
end

return

function [dots,framestart,shade] = dot_compute(dots,speed,shade,square_dims,dot_dia,phantom_dots,windowSize,repulse,framestart,repulsor_force)

% parameters
fps = 30;
%speed = 25;
%speed = mean(windowSize.pixels)*2;
%repfulse = 7/3;
if speed <= 0
    shade = abs(speed);
    speed = 0.3; %make it so the dots never actually stop moving
    %take gray amount as distance from 3;
end
rate = speed * (1 / fps) / mean(windowSize.degrees_per_pixel);
%rate = mean(windowSize.pixels)*2/fps; this is the eqn to get 1 Hz
%motion****
%speed_limit = speed * (2 / fps) / mean(windowSize.degrees_per_pixel); % frame-wise speed limit in degrees (0.119 for 33.33ms draw rate)
repulsor_focus = 2;
%repulsor_force = 0;
%repulsor_force = speed * 0.1;
%repulsor_force = 2;
repulsor_distance = dot_dia*repulse; %4/3 first
bumper_limit = dot_dia*1.5;
%bumper_limit = dot_dia*proximity;
% initialize
num_dots = length(dots);
prev_trial = find(~isnan(dots(1).pos(:,1)),1,'last');
trial = prev_trial + 1;
dotList = zeros(num_dots,2);

% randomly update all dot trajectories
for i = 1:num_dots
    
    %dots(i).trajectory(trial,:) = dots(i).trajectory(prev_trial,:) + (([rand() rand()]-0.5)*rate); %negative .5 to +.5
    %dots(i).trajectory(trial,:) = dots(i).trajectory(prev_trial,:);
    vx = (rand-0.5)*(rate);
    vy = sign(rand-0.5)*sqrt(rate^2 - vx^2);
    if prev_trial == 1
        dots(i).trajectory(trial,:) = [vx vy];
    else
        dots(i).trajectory(trial,:) = dots(i).trajectory(prev_trial,:); %keep it in the same direction unless acted upon by outside force!!
    end
    dots(i).pos(trial,:) = dots(i).pos(prev_trial,:) + dots(i).trajectory(trial,:) ;
    dotList(i,:) = dots(i).pos(trial,:); % accumulate dot positions
end
%bigDotList = [dotList; phantom_dots];
bigDotList = [dotList];
% imminent collision detector
distances = ipdm(dotList) + (eye(size(dotList,1)) + tril(ones(size(dotList,1))))*bumper_limit;
collisions = distances < bumper_limit;
if any(any(collisions))
    [dot1 dot2] = find(collisions); %do this for every collision
    for i=1:length(dot1)
        %first we're going to find the direction that they collided with on
        %the past trial
        % to catch bad behavior, we're assuming they collided in prev_trial
        %PASS = 0;
        if trial > 2 %check for bad behavior
            
            D1 = dots(dot1(i)).pos(trial-1:trial,:);
            D2 = dots(dot2(i)).pos(trial-1:trial,:);
            
            lastDistance = pdist([D1(1,:); D2(1,:)]);
            currentDistance = pdist([D1(2,:); D2(2,:)]);
            D1pos = D1(2,:);
            D2pos = D2(2,:);
            subtractpositions = abs(D1 - D2);
            %find which direction they're colliding in
            
            %velocities = [dots(dot1(i)).trajectory(trial,:); dots(dot2(i)).trajectory(trial,:)];
            %x_collision = sign(prod(velocities(:,1))) < 0;
            %y_collision = sign(prod(velocities(:,2))) < 0;
            x_collision = subtractpositions(1,1) - subtractpositions(2,1) > 0;
            y_collision = subtractpositions(1,2) - subtractpositions(2,2) > 0;
            %collision on the last refresh
            %xvel = [dots(dot1(i)).trajectory(trial-2,1) dots(dot1(i)).trajectory(trial-1,1)];
            % yvel = [dots(dot1(i)).trajectory(trial-2,1) dots(dot1(i)).trajectory(trial-1,1)];
            
            if currentDistance < lastDistance %|| (sign(prod(xvel)) > 0) || (sign(prod(yvel)) > 0) %this means they're getting closer
                % 3/20/16: try to make it so dots ALWAYS move away from
                % each other instead of just reversing direction
                
                
                if x_collision %&&  sign(prod(xvel(:,1))) > 0;%check if they didn't also flip signs in x velocity 2 trials ago
                    if D1pos(1) < D2pos(1) %if first ball is to the left, move it to the left
                        dots(dot1(i)).trajectory(trial,1) = abs(dots(dot1(i)).trajectory(prev_trial,1)) * -1;
                        dots(dot2(i)).trajectory(trial,1) = abs(dots(dot2(i)).trajectory(prev_trial,1));
                    else %if the first dot is to the right
                        dots(dot1(i)).trajectory(trial,1) = abs(dots(dot1(i)).trajectory(prev_trial,1));
                        dots(dot2(i)).trajectory(trial,1) = abs(dots(dot2(i)).trajectory(prev_trial,1)) * -1;
                        %dots(dot1(i)).trajectory(trial,1) = dots(dot1(i)).trajectory(prev_trial,1) * -1;
                        %dots(dot2(i)).trajectory(trial,1) = dots(dot2(i)).trajectory(prev_trial,1) * -1;
                    end
                end
                if y_collision %&& sign(prod(yvel(:,1))) > 0;
                    if D1pos(2) < D2pos(2) %if the first dot is higher
                        dots(dot1(i)).trajectory(trial,2) = abs(dots(dot1(i)).trajectory(prev_trial,2)) * -1;
                        dots(dot2(i)).trajectory(trial,2) = abs(dots(dot2(i)).trajectory(prev_trial,2));
                    else
                        dots(dot1(i)).trajectory(trial,2) = abs(dots(dot1(i)).trajectory(prev_trial,2));
                        dots(dot2(i)).trajectory(trial,2) = abs(dots(dot2(i)).trajectory(prev_trial,2)) * -1;
                    end
                    
                    %dots(dot1(i)).trajectory(trial,2) = dots(dot1(i)).trajectory(prev_trial,2) * -1;
                    %dots(dot2(i)).trajectory(trial,2) = dots(dot2(i)).trajectory(prev_trial,2) * -1;
                    
                end
                %               if ~x_collision && ~y_collision
                %                   %stop here and figure it out
                %                   %see why this would be the case--look for clues
                %                   daldsaljsld
                %               end
                
            else %make sure if you're not reversing that the dots aren't moving towards one another
                fprintf('saved?!?!?!?')
                % framestart(i,1) = trial;
                % framestart(i,2) = dot1;
                % framestart(i,3) = dot2;
                %dots(dot1(i)).trajectory(trial,:) = dots(dot1(i)).trajectory(prev_trial,:);
                %dots(dot2(i)).trajectory(trial,:) = dots(dot2(i)).trajectory(prev_trial,:);
            end
            % dots(dot1(i)).pos(trial,:) = dots(i).pos(prev_trial,:) + dots(dot1(i)).trajectory(trial,:);
        end
    end
end
% if  ~isempty(framestart) && trial >= framestart(1,1) + 10
%    % sldlad
% elseif isempty(framestart)
%     framestart = [];
% end
% reflect at the edges

%repulsor_force = 40;
% apply force fields based on dot distances
distances = ipdm(bigDotList) - repulsor_distance;
distances(distances<0) = NaN;
for i = 1:num_dots
    my_pos = dots(i).pos(trial,:);
    pos_diff = -1 * bsxfun(@minus,bigDotList,my_pos);
    for j = 1:length(pos_diff)
        pos_diff_norm(j,:) = pos_diff(j,:)/norm(pos_diff(j,:));
    end
    dist_factor = distances(i,:) .^ -repulsor_focus; %magnitude difference
    vecs = bsxfun(@times,pos_diff_norm,dist_factor');
    %they hit and then here tellin?
    %repulsor_force*sum(vecs(~isnan(vecs(:,1)),:),1);
    newroute = dots(i).trajectory(trial,:) + repulsor_force*sum(vecs(~isnan(vecs(:,1)),:),1);
    
    C = sqrt(rate^2/(newroute(1)^2 + newroute(2)^2));
    dots(i).trajectory(trial,:) = C*newroute;
    %dots(i).trajectory(trial,:) = ratio*newroute;
end

edge_dist = 2/3;

for i = 1:num_dots
    reflect =  [(dots(i).pos(trial,:) < dot_dia*edge_dist) ; (dots(i).pos(trial,:) > square_dims - dot_dia*edge_dist)];
    %reflect = -1 * ((dots(i).pos(trial,:) < dot_dia/1.5) + (dots(i).pos(trial,:) > square_dims - dot_dia/1.5));
    if sum(sum(reflect))~=0
        [rows cols] = find(reflect~=0);
        if sum(reflect(1,:)) ~=0 %this means too small so set trajectory positive in that direction
            coords = find(reflect(1,:)~=0);
            dots(i).trajectory(trial,coords) = abs(dots(i).trajectory(trial,coords)); %make positive
        elseif sum(reflect(2,:)~=0) %here too big so make negative
            coords = find(reflect(2,:)~=0);
            dots(i).trajectory(trial,coords) = -1*abs(dots(i).trajectory(trial,coords));
        end
        %reflect(reflect == 0) = 1;
        %dots(i).trajectory(trial,:) = dots(i).trajectory(prev_trial,:) .* reflect;
        % dots(i).pos(trial,:) = dots(i).pos(prev_trial,:) + dots(i).trajectory(trial,:);
    end
end

%     % don't break the speed limit
%     for i = 1:num_dots
%         speeding_ticket = abs(dots(i).trajectory(trial,:)) > speed_limit;
%         if any(speeding_ticket);
%             bigger_val = max(abs(dots(i).trajectory(trial,:)));
%             dots(i).trajectory(trial,:) = dots(i).trajectory(trial,:) / bigger_val * speed_limit;
%             dots(i).pos(trial,:) = dots(i).pos(prev_trial,:) + dots(i).trajectory(trial,:);
%         end
%     end
%
%
return

function [dots phantom_dots] = initialize_dots(num_dots,num_targets,square_dims,dot_dia)

% params
min_boundary = dot_dia;
min_distance = dot_dia*2;
safe_dims = square_dims - (min_boundary*2);

% phantom dots
phantom_dots(1,:) = square_dims / 2;
phantom_dots(end+1,:) = [1 1];
phantom_dots(end+1,:) = square_dims;
phantom_dots(end+1,:) = [square_dims(1) 1];
phantom_dots(end+1,:) = [1 square_dims(2)];

% initialize
is_good = false;

% keep trying to generate dots until we get it right
while ~is_good
    dotPosList = phantom_dots; targetList = false(length(phantom_dots),1);
    
    % create new dots iteratively
    for dot = 1:num_dots
        % begin will null metadata
        if dot <= num_targets
            dots(dot).is_target = true;
        else
            dots(dot).is_target = false;
        end
        dots(dot).is_probe = false;
        dots(dot).pos = nan(600,2);
        dots(dot).trajectory = nan(600,2);
        dots(dot).trajectory(1,:) = 0;
        
        % initialize position
        good_dot = false;
        while ~good_dot
            for d = 1:2
                dots(dot).pos(1,d) = rand()*safe_dims(d)+min_boundary;
            end
            distance = ipdm(dots(dot).pos(1,:),dotPosList);
            if all(distance > min_distance)
                good_dot = true;
            end
        end
        
        % log the new dot
        dots(dot).log = dots(dot).pos(1,:);
        dotPosList = [dotPosList; dots(dot).pos(1,:)];
        targetList = [targetList; dots(dot).is_target];
    end % dot
    
    % on-multi-target trials, reject the selection if there is bunching of targets
    if sum(targetList) > 1
        % assess general dot distances in ten 5-subset samples of dots
        for i = 1:10
            for j = 1:5
                choice(j) = randi(10);
            end
            general_distance(i) = mean(mean(ipdm(dotPosList(choice,:))));
        end
        general_distance = mean(general_distance);
        
        % compare this to the distance in our target subsample
        inter_targ_distance = mean(mean(ipdm(dotPosList(targetList,:))));
        if abs(inter_targ_distance - general_distance) < (general_distance * 0.2)
            is_good = true;
        end
    else
        is_good = true;
    end
end % is_good

return


function [acc rt timefirst] = odd_even(digitsEK,num_qs,q_dur,isi_dur,min_format,window,keys,COLORS,DEVICE,SUBJ_NAME, nest_info,SLACK,timespecFirst, keysObject)

% initialize
old_keys = [1 5];
triggerNext = false; %changed because presentation has to be the same!!!!
console = false;
files = false;

acc = nan; rt = nan;
CORRECT_COLOR = COLORS.GREEN;
INCORRECT_COLOR = COLORS.RED;
% loop over digit probes
for digit_q = 1:num_qs
    
    % prepare trial data
    digit_picks = [0 0];
    while sum(digit_picks) <= 10
        digit_picks = [randi(9) randi(9)];
    end
    is_odd = mod(sum(digit_picks),2);
    if min_format
        digit_prompt = [num2str(digit_picks(1)) '  +  ' num2str(digit_picks(2))];
    else
        digit_prompt = [num2str(digit_picks(1)) '  +  ' num2str(digit_picks(2)) '\n\n'];
        prompt_white = '\n\nEVEN                              ODD';
    end
    cresp = is_odd + 1;
    
    % execute
    DrawFormattedText(window,digit_prompt,'center','center',COLORS.MAINFONTCOLOR,0);
    if ~min_format
        DrawFormattedText(window,prompt_white,'center','center',COLORS.MAINFONTCOLOR,0);
    end
    if digit_q == 1
        onset = Screen('Flip',window,timespecFirst);
        timefirst = onset;
    else
        timespec = ision(digit_q-1) + isi_dur - SLACK;
        onset = Screen('Flip',window,timespec);
    end
    digitsEK = easyKeys(digitsEK, ...
        'onset', onset, ...
        'stim', digit_prompt, ...
        'nesting', [nest_info digit_q], ...
        'cresp', keys(cresp), 'cresp_map', keysObject.map(old_keys(cresp),:), 'valid_map', sum(keysObject.map(old_keys,:)));
    timespec = onset + q_dur - SLACK;
    ision(digit_q) = isi_specific(window,COLORS.MAINFONTCOLOR,timespec);
end % loop over q's

% summary stats
acc = mean(digitsEK.trials.acc);
valid_rts = ~isnan(digitsEK.trials.rt);
if sum(valid_rts) > 0
    rt = mean(digitsEK.trials.rt(valid_rts));
else rt = nan;
end


return