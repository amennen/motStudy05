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
RECALL2 = MOT{end} + 1; % post-scan rsvp memory test
DESCRIPTION = RECALL2 + 1; %26
ASSOCIATES = DESCRIPTION + 1; %27
base_path = [fileparts(which('mot_realtime04MB.m')) filesep];


%% look at descriptive ratings
subjectVec = [4:6];
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
        easy = find(cond==2);
        hard = find(cond==1);
        easy_score(i) = nanmean(rating(easy));
        hard_score(i) = nanmean(rating(hard));
    end
    easyRatingDiff(s) = diff(easy_score);
    hardRatingDiff(s) = diff(hard_score);
    %% look at AB' d' ratings
    % YS suggestion: look at the difference between presented first and
    % presented after if related--so if getting right once will mess up getting
    % right the other time
    stimFile = dir(fullfile(behavioral_dir, ['mot_*' num2str(ASSOCIATES) '*.mat']));
    sf = load(fullfile(behavioral_dir, stimFile(end).name));
    Afirst = sf.stim.AAPID;
    id = sf.stim.id;
    nPractice = 3;
    for i = 1:length(sf.stim.id)
        if ismember(sf.stim.id(i),sf.stim.AAPID)
            if i <=20+nPractice
                c(i) = 1; % target
            else
                c(i) = 2; % lure
            end
        else
            if i<=20+nPractice
                c(i) = 2; % lure
            else
                c(i) = 1; % target
            end
        end
    end
    r = dir(fullfile(behavioral_dir, ['_RECOG' '*.mat']));
    r = load(fullfile(behavioral_dir,r(end).name));
    trials = table2cell(r.datastruct.trials);
    stimID = cell2mat(trials(:,8));
    cond = cell2mat(trials(:,9));
    acc = cell2mat(trials(:,11));
    rt = cell2mat(trials(:,13));
    cresp = cell2mat(trials(:,22));
    
    easy = find(cond==2);
    hard = find(cond==1);
    target = find(cresp==1);
    lure = find(cresp==2);
    
    hardTargets = intersect(hard,target);
    hardLures = intersect(hard,lure);
    easyTargets = intersect(easy,target);
    easyLures = intersect(easy,lure);
    % out of hard trials, find rt to familiar items
    HT_RT(s) = nanmedian(rt(hardTargets));
    ET_RT(s) = nanmedian(rt(easyTargets));
    HL_RT(s) = nanmedian(rt(hardLures));
    EL_RT(s) = nanmedian(rt(easyLures));
    
    % should make sure it's not influenced by order for RT******
    % if pair appears in first half, the index would be less than 23
    hardTargetFirst = hardTargets(hardTargets<=23);
    hardTargetSecond = hardTargets(hardTargets>23);
    hardLureFirst = hardLures(hardLures<=23);
    hardLureSecond = hardLures(hardLures>23);
    easyTargetFirst = easyTargets(easyTargets<=23);
    easyTargetSecond = easyTargets(easyTargets>23);
    easyLureFirst = easyLures(easyLures<=23);
    easyLureSecond = easyLures(easyLures>23);
    
    % average over all HTF
    HTF(s) = nanmean(rt(hardTargetFirst));
    HTS(s) = nanmean(rt(hardTargetSecond));
    HLF(s) = nanmean(rt(hardLureFirst));
    HLS(s) = nanmean(rt(hardLureSecond));
    
    ETF(s) = nanmean(rt(easyTargetFirst));
    ETS(s) = nanmean(rt(easyTargetSecond));
    ELF(s) = nanmean(rt(easyLureFirst));
    ELS(s) = nanmean(rt(easyLureSecond));
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
subjectVec = [4:8];
nresp = 3*10;
allresp = zeros(nresp,4);
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

        if s == 1 && length(stimID) ~= 40 % then we have to add another row
            beg = stimID(1:7,:);
            part2 = stimID(8:end,:);
            beg(8,:) = stimID(7,:);
            stimID = [beg;part2];
            begresp = resp(1:7,:);
            part2resp = resp(8:end,:);
            begresp(8,:) = NaN;
            resp = [begresp;part2resp];
        end
        IDmat = reshape(stimID,4,length(stimID)/4)';
        [~,indSort] = sort(IDmat(:,1));
        IDorder = IDmat(indSort,:);
        RESPmat = reshape(resp,4,length(resp)/4)';
        RESPorder = RESPmat(indSort,:);
        %allresp((i-1)*10 + 1: i*10,:,s ) = RESPorder;
    end

   
    
end