% analyze behav

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


%% look at descriptive ratings
subjectVec = [4];
for s = 1:length(subjectVec)
    subjectNum = subjectVec(s);
    behavioral_dir = [base_path 'BehavioralData/' num2str(subjectNum) '/'];
    
    recallSession = [RECALL1 RECALL2];
    
    for i = 1:2
        r = dir(fullfile(behavioral_dir, ['EK' num2str(recallSession(i)) '_' 'SUB'  '*.mat']));
        r = load(fullfile(behavioral_dir,r(end).name));
        trials = table2cell(r.datastruct.trials);
        stimID = cell2mat(trials(:,8));
        cond = cell2mat(trials(:,9));
        rating = cell2mat(trials(:,12));
        [~,indSort] = sort(stimID);
        ratingOrdered = rating(indSort);
        condOrdered = cond(indSort);
        easy = find(condOrdered==2);
        hard = find(condOrdered==1);
        
        easyScores(i,1:10) = ratingOrdered(easy);
        hardScores(i,1:10) = ratingOrdered(hard);
        %easy_score(i) = nanmean(rating(easy));
        %hard_score(i) = nanmean(rating(hard));
    end
    %easyRatingDiff(s) = diff(easy_score);
    %hardRatingDiff(s) = diff(hard_score);
    save(fullfile(behavioral_dir, 'ratings.mat'), 'easyScores', 'hardScores')
end
%% now graph
figure;
dataH = [nanmedian(HTF) nanmedian(HTS) nanmedian(HLF) nanmedian(HLS)];
dataE = [nanmedian(ETF) nanmedian(ETS) nanmedian(ELF) nanmedian(ELS)];
plot(1:4,dataH, 'r.', 1:4, dataE, 'k.', 'MarkerSize', 12)
firstgroup = [HTF; ETF];
secondgroup = [HTS; ETS];
avgratio = [nanmean(firstgroup,2) nanmean(secondgroup,2)];
eavgratio = [nanstd(firstgroup,[],2)/sqrt(length(subjectVec)-1) nanstd(secondgroup,[],2)/sqrt(length(subjectVec)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['MOT ';'OMIT']);
legend('Target First', 'Target Second')
xlabel('Stimulus Type')
ylabel('Recognition RT')
title('Recognition Accuracy')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylim([0 1])
%%
%%figure;
dataH = [nanmedian(HLF) nanmedian(HLS) nanmedian(HLF) nanmedian(HLS)];
dataE = [nanmedian(ELF) nanmedian(ELS) nanmedian(ELF) nanmedian(ELS)];
plot(1:4,dataH, 'r.', 1:4, dataE, 'k.', 'MarkerSize', 12)
firstgroup = [HLF; ELF];
secondgroup = [HLS; ELS];
avgratio = [nanmean(firstgroup,2) nanmean(secondgroup,2)];
eavgratio = [nanstd(firstgroup,[],2)/sqrt(length(subjectVec)-1) nanstd(secondgroup,[],2)/sqrt(length(subjectVec)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['MOT ';'OMIT']);
legend('Lure First', 'Lure Second')
xlabel('Stimulus Type')
ylabel('Recognition RT')
title('Recognition Accuracy')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylim([0 1])

%print(thisfig, sprintf('%sallrecogAcc.pdf', allplotDir), '-dpdf')

%%
%% look at descriptive ratings
subjectVec = [3];
nresp = 3*10;
allresp = zeros(10,4*3);
for s = 1:length(subjectVec)
    subjectNum = subjectVec(s);
    behavioral_dir = [base_path 'BehavioralData/' num2str(subjectNum) '/'];
    for i = 1:3
        SESSION = MOT{i};
        r = dir(fullfile(behavioral_dir, ['EK' num2str(SESSION) '_' 'SUB'  '*.mat']));
        r = load(fullfile(behavioral_dir,r(end).name));
        trials = table2cell(r.datastruct.trials);
        stimID = cell2mat(trials(:,8));
        resp = cell2mat(trials(:,12));


        IDmat = reshape(stimID,4,length(stimID)/4)';
        [~,indSort] = sort(IDmat(:,1));
        IDorder = IDmat(indSort,:);
        RESPmat = reshape(resp,4,length(resp)/4)';
        RESPorder = RESPmat(indSort,:);
        allresp(:, (i-1)*4 + 1: i*4 ) = RESPorder;
    end
    save(fullfile(behavioral_dir, 'RTresponses.mat'), 'allresp')
end
%%
%% now convert recog to cell
subjects = [4];
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