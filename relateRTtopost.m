% purpose: relate MOT trials to how the stimuli did afterwards

% for each subject
% for each stimuli:
%   look a number of TR's spent high (0.1+), medium vals (-0.1-+0.1), neg
%   (-0.1-)
% then look at the post score (and pre-post score)
%(this should be for all TR's--not just the FB ones

%want as an output:
% ACTUALLY: what would be most helpful: for each stimuli: average
% activation for recall1, then average activation in MOT, then average
% activation in recall2 [then can maybe understand how it's changing]

nstim = 10;
nTRs = 15;
projectName = 'motStudy03';
allplotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' ];
cd('/Data1/code/motStudy03/code/')
nblock = 3;
shiftTR = 2;
svec = [3:9 11 12 13];
nsub = length(svec);
for s = 1:nsub
    allsep = [];
    subjectNum = svec(s);
    for iblock = 1:nblock
        blockNum = iblock;
        SESSION = 20 + blockNum;
        behavioral_dir = ['/Data1/code/' projectName '/code/BehavioralData/' num2str(subjectNum) '/'];
        save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/'];
        classOutputDir = fullfile(save_dir,['motRun' num2str(blockNum)], 'classOutput/');
        runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
        fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime02_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
        plotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' num2str(subjectNum) '/'];
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end
        matlabOpenFile = [behavioral_dir '/' fileSpeed(end).name];
        d = load(matlabOpenFile);
        
        % now find the category separation for that block
        allMotionTRs = convertTR(d.timing.trig.wait,d.timing.plannedOnsets.motion,d.config.TR); %row,col = mTR,trialnumber
        TRvector = reshape(allMotionTRs,1,numel(allMotionTRs));
        for i=1:length(TRvector)
            fileTR = TRvector(i) + 2;
            [~, tempfn{fileTR}] = GetSpecificClassOutputFile(classOutputDir,fileTR);
            tempStruct = load(fullfile(classOutputDir, tempfn{fileTR}));
            categsep(i) = tempStruct.classOutput;
        end
        sepbytrial = reshape(categsep,15,10);
        sepbytrial = sepbytrial'; %results by trial number, TR number
        sepvec = reshape(sepbytrial,1,numel(sepbytrial));
        [~,indSort] = sort(d.stim.id);
        sepinorder = sepbytrial(indSort,:);
        %test if fb only
        sepbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = sepinorder;
    end
    %average separation over all trials
    avgsepbystim(:,s) = mean(sepbystim,2);
    [stimL,j ] = find(sepbystim<-0.1);
    [stimM,j] = ind2sub(size(sepbystim),intersect(find(sepbystim>=-0.1) ,find(sepbystim<=0.1)));
    %for z = 1:length(stimM)
    %    sepbystim(stimM(z),j(z))
    %end
    [stimH,j] = find(sepbystim>0.1);

    for p = 1:10
        nlow(p,s) = length(find(stimL==p));
        nmed(p,s) = length(find(stimM==p));
        nhigh(p,s) = length(find(stimH==p));
    end
    %% now for the same stimuli load the pre/post patterns: change recall file process and just save the output?
    % code that just runs recall file process for each subject so you can just
    % look it up later with the .mat file
    recallSessions = [20 24];
    nTRptrial = 4;
    for i = 1:2
        session = recallSessions(i);
        rfn = dir(fullfile(save_dir, ['recallpatternsdata_' num2str(session)  '*.mat']));
        rp = load(fullfile(save_dir,rfn(end).name));
        rpat = rp.patterns;
        [expat,trials,stimOrder] = GetSessionInfoRT(subjectNum,session,behavioral_dir);
        testTrials = find(any(expat.regressor.allCond));
        allcond = rpat.regressor.allCond(:,testTrials);
        categSep = rpat.categSep(:,testTrials);
        %categSep = rpat.categSep(:,union(testTrials,testTrials+shiftTR)); %all testTR's plus 2 before
        
        z = reshape(categSep,nTRptrial,20); %for 20 trials --make sure this works here!
        byTrial = z';
        RTtrials = byTrial(trials.hard,:);
        RTtrials = RTtrials(stimOrder.hard,:);
        OMITtrials = byTrial(trials.easy,:);
        OMITtrials = OMITtrials(stimOrder.easy,:);
        
        RTevidence(:,:,i) = RTtrials;
        OMITevidence(:,:,i) = OMITtrials;
    end
    PrePostRT = RTevidence(:,:,2) - RTevidence(:,:,1);
    avgPreRT(:,s) = mean(RTevidence(:,:,1),2);
    PrePostOMIT = OMITevidence(:,:,2) - OMITevidence(:,:,1);
    PrePostRTbyT(s,:) = mean(PrePostRT);
    avgPrePostOM(:,s) = mean(PrePostOMIT,2);
    %average PrePostRT over all 4 TR's
    avgPrePostRT(:,s) = mean(PrePostRT,2);
    avgPostRT(:,s) = mean(RTevidence(:,:,2),2);
    PrePostOMbyT(s,:) = mean(PrePostOMIT);
end
%%
%now save data and send to jukebox
folder= '/jukebox/norman/amennen/PythonMot3';
fn = 'plottInfo.mat';
save(fn,'avgsepbystim','PrePostRT','avgPrePostOM', 'avgPrePostRT','avgPostRT', 'diff_stim_hard','all_hard','all_easy','avgPreRT', 'hardacc_ordered', 'hardrt_ordered', 'PrePostRTbyT', 'PrePostOMbyT')
unix(['scp ' fn ' amennen@apps.pni.princeton.edu:' folder '/'])
newfn = 'newPrePost.mat';
save(newfn, 'PrePostRTbyT', 'PrePostOMbyT')
unix(['scp ' newfn ' amennen@apps.pni.princeton.edu:' folder '/'])

%% now plot the whole thing
h = figure;
for s = 1:nsub
%scatter(avgsepbystim(:,s),avgPrePostRT(:,s),  'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
p = scatter(avgsepbystim(:,s),avgPrePostRT(:,s),50, 'fill');
alpha(p,0.1)
p = patch(avgsepbystim(:,s),avgPrePostRT(:,s), 'c');
set(p, 'EdgeColor', none)
hold on;
end
xlim([-0.25 0.25])
ylim([-.5 0.8])
x = reshape(avgsepbystim,nsub*10,1);
y = reshape(avgPrePostRT,nsub*10,1);
[rho,pval] = corrcoef([x y]);
text(0.1,.6,['corr = ' num2str(rho(1,2))]);
text(0.1,.45, ['p = ' num2str(pval(1,2))]);
title('PostPre Recall vs. MOT Evidence')
xlabel('Average Category Evidence During RT')
ylabel('Average Evidence Difference Score Post-Pre')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepPrePost.pdf', allplotDir), '-dpdf')


h = figure;
for s = 1:nsub
    scatter(avgsepbystim(:,s),avgPostRT(:,s), 'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
    hold on;
end
xlim([-0.25 0.25])
ylim([-.5 0.8])
x = reshape(avgsepbystim,nsub*10,1);
y = reshape(avgPostRT,nsub*10,1);
[rho,pval] = corrcoef([x y]);
text(0.1,.6,['corr = ' num2str(rho(1,2))]);
text(0.1,.45, ['p = ' num2str(pval(1,2))]);
title('Post Recall vs. MOT Evidence')
xlabel('Average Category Evidence During RT')
ylabel('Average Post Evidence After')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepPost.pdf', allplotDir), '-dpdf')

%% see if relationship with number of TR's
h = figure;
for s = 1:nsub
scatter(nlow(:,s),avgPrePostRT(:,s), 'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
hold on;
end
xlim([0 30])
ylim([-.5 0.8])
x = reshape(nlow,nsub*10,1);
y = reshape(avgPrePostRT,nsub*10,1);
[rho,pval] = corrcoef([x y]);
text(0.1,.6,['corr = ' num2str(rho(1,2))]);
text(0.1,.45, ['p = ' num2str(pval(1,2))]);
title('PostPre Recall vs. nLow Evidence RT')
xlabel('n TRs < -0.1')
ylabel('Average Evidence Difference Score Post-Pre')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepnLow.pdf', allplotDir), '-dpdf')

h = figure;
for s = 1:nsub
scatter(nmed(:,s),avgPrePostRT(:,s),'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
hold on;
end
xlim([0 30])
ylim([-.5 0.8])
x = reshape(nmed,nsub*10,1);
y = reshape(avgPrePostRT,nsub*10,1);
[rho,pval] = corrcoef([x y]);
text(0.1,.6,['corr = ' num2str(rho(1,2))]);
text(0.1,.45, ['p = ' num2str(pval(1,2))]);
title('PostPre Recall vs. nMed Evidence RT')
xlabel('n TRs between -0.1:0.1')
ylabel('Average Evidence Difference Score Post-Pre')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepnMed.pdf', allplotDir), '-dpdf')

h = figure;
for s = 1:nsub
scatter(nhigh(:,s),avgPrePostRT(:,s), 'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
hold on;
end
xlim([0 30])
ylim([-.5 0.8])
x = reshape(nhigh,nsub*10,1);
y = reshape(avgPrePostRT,nsub*10,1);
[rho,pval] = corrcoef([x y]);
text(0.1,.6,['corr = ' num2str(rho(1,2))]);
text(0.1,.45, ['p = ' num2str(pval(1,2))]);
title('PostPre Recall vs. High Evidence RT')
xlabel('n TRs > 0.1')
ylabel('Average Evidence Difference Score Post-Pre')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepnHigh.pdf', allplotDir), '-dpdf')


%%
h = figure;
for s = 1:nsub
scatter(avgsepbystim(:,s),diff_stim_hard(:,s), 'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
hold on;
end
xlim([-.2 .2])
ylim([-3 3])
y = reshape(diff_stim_hard,nsub*10,1);
x = reshape(avgsepbystim,nsub*10,1);
[rho,pval] = corrcoef([x(~isnan(y)) y(~isnan(y))]);
text(0.1,2,['corr = ' num2str(rho(1,2))]);
text(0.1,1.5, ['p = ' num2str(pval(1,2))]);
title('Difference in Ratings vs. Avg Evidence MOT')
xlabel('Avg MOT Evidence')
ylabel('Detail Rating Score Post-Pre')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepRating.pdf', allplotDir), '-dpdf')

h = figure;
for s = 1:nsub
scatter(avgsepbystim(:,s),hardacc_ordered(:,s), 'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
hold on;
end
xlim([-.2 .2])
ylim([0 1])
y = reshape(hardacc_ordered,nsub*10,1);
x = reshape(avgsepbystim,nsub*10,1);
[rho,pval] = corrcoef([x(~isnan(y)) y(~isnan(y))]);
text(0.1,.8,['corr = ' num2str(rho(1,2))]);
text(0.1,.75, ['p = ' num2str(pval(1,2))]);
title('Recog Acc vs. Avg Evidence MOT')
xlabel('Avg MOT Evidence')
ylabel('Recognition Accuracy')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepAcc.pdf', allplotDir), '-dpdf')

h = figure;
for s = 1:nsub
scatter(avgsepbystim(:,s),hardrt_ordered(:,s), 'jitter','on', 'jitterAmount',0.005,'LineWidth', 2);
hold on;
end
xlim([-.2 .2])
ylim([0 1])
y = reshape(hardrt_ordered,nsub*10,1);
x = reshape(avgsepbystim,nsub*10,1);
[rho,pval] = corrcoef([x(~isnan(y)) y(~isnan(y))]);
text(0.1,.9,['corr = ' num2str(rho(1,2))]);
text(0.1,.85, ['p = ' num2str(pval(1,2))]);
title('Recognition RT vs. Avg Evidence MOT')
xlabel('Avg MOT Evidence')
ylabel('Recognition RT')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%savgsepRT.pdf', allplotDir), '-dpdf')
%%
h = figure;
for s = 1:nsub
scatter(all_easy(:,s),avgPrePostOM(:,s), 'LineWidth', 2.5);
hold on;
end
xlim([1 5])
ylim([-1 1])
x = reshape(all_easy,nsub*10,1);
y = reshape(avgPrePostOM,nsub*10,1);
[rho,pval] = corrcoef([x(~isnan(x)) y(~isnan(x))]);
text(4,.8,['corr = ' num2str(rho(1,2))]);
text(4,.7, ['p = ' num2str(pval(1,2))]);
title('OMIT: PostPreEv vs. Recall 1 Rating')
xlabel('Recall 1 Rating')
ylabel('Average Post-Pre Evidence')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%sOMrating_sep.pdf', allplotDir), '-dpdf')

h = figure;
for s = 1:nsub
scatter(all_hard(:,s),avgPrePostRT(:,s), 'LineWidth', 2.5);
hold on;
end
xlim([1 5])
ylim([-1 1])
x = reshape(all_hard,nsub*10,1);
y = reshape(avgPrePostRT,nsub*10,1);
[rho,pval] = corrcoef([x(~isnan(x)) y(~isnan(x))]);
text(4,.8,['corr = ' num2str(rho(1,2))]);
text(4,.7, ['p = ' num2str(pval(1,2))]);
title('RT: PostPreEv vs. Recall 1 Rating')
xlabel('Avg MOT Evidence')
ylabel('Recognition RT')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
print(h, sprintf('%sRTrating_sep.pdf', allplotDir), '-dpdf')

%% NOW GET DATA WITH OLD EXPERIMENT
oldDir = '/Data1/code/motStudy02/code/';
cd(oldDir)
projectName = 'motStudy02';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;

subjectVec = [8 12 14 15 16 18 20 22 26 27 28 30 31 32]; %now made it so for all subjects and can separate them into RT/YC afterwards in python 2/17

nsub = length(svec);

for s = 1:nsub
    allsep = [];
    subjectNum = svec(s);
    
    
    %% now for the same stimuli load the pre/post patterns: change recall file process and just save the output?
    % code that just runs recall file process for each subject so you can just
    % look it up later with the .mat file
    recallSessions = [19 23];
    nTRptrial = 4;
    for i = 1:2
        session = recallSessions(i);
        rfn = dir(fullfile(save_dir, ['recallpatternsdata_' num2str(session)  '*.mat']));
        rp = load(fullfile(save_dir,rfn(end).name));
        rpat = rp.patterns;
        [expat,trials,stimOrder] = GetSessionInfoRT(subjectNum,session,behavioral_dir);
        testTrials = find(any(expat.regressor.allCond));
        allcond = rpat.regressor.allCond(:,testTrials);
        categSep = rpat.categSep(:,testTrials);
        %categSep = rpat.categSep(:,union(testTrials,testTrials+shiftTR)); %all testTR's plus 2 before
        
        z = reshape(categSep,nTRptrial,20); %for 20 trials --make sure this works here!
        byTrial = z';
        RTtrials = byTrial(trials.hard,:);
        RTtrials = RTtrials(stimOrder.hard,:);
        OMITtrials = byTrial(trials.easy,:);
        OMITtrials = OMITtrials(stimOrder.easy,:);
        
        RTevidence(:,:,i) = RTtrials;
        OMITevidence(:,:,i) = OMITtrials;
    end
    PrePostRT = RTevidence(:,:,2) - RTevidence(:,:,1);
    avgPreRT(:,s) = mean(RTevidence(:,:,1),2);
    PrePostOMIT = OMITevidence(:,:,2) - OMITevidence(:,:,1);
    PrePostRTbyT(s,:) = mean(PrePostRT);
    avgPrePostOM(:,s) = mean(PrePostOMIT,2);
    %average PrePostRT over all 4 TR's
    avgPrePostRT(:,s) = mean(PrePostRT,2);
    avgPostRT(:,s) = mean(RTevidence(:,:,2),2);
    PrePostOMbyT(s,:) = mean(PrePostOMIT);
end

