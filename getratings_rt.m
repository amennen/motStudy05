%calculate subjective details during MOT
%eh a lot of noise here--next: check whether the 
%cd /Volumes/norman/amennen/behav_test_anne/Participant' Data'/1
%number of participants here

%look at the recognition memory at the end and listen to wav files! (use
%recogdata.m to look at the recognition memory)
%clear all;
projectName = 'motStudy03';
base_path = [fileparts(which('mot_realtime02.m')) filesep];

% don't put in 22 until have subject
svec = [3 4 5 6 7 8 9 11 12 13];

nsub = length(svec);
recallSession = [20 24];
nstim = 10;
allplotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' ];


for s = 1:nsub
    behavioral_dir = [base_path 'BehavioralData/' num2str(svec(s)) '/'];
    for i = 1:length(recallSession)
        r = dir(fullfile(behavioral_dir, ['EK' num2str(recallSession(i)) '_' 'SUB'  '*.mat'])); 
        r = load(fullfile(behavioral_dir,r(end).name)); 
        trials = table2cell(r.datastruct.trials);
        stimID = cell2mat(trials(:,8));
%         
%         if find(RT_m == svec(s)) %then it's in the RT group
%             matched = find(RT_m == svec(s));
%         elseif find(YC_m == svec(s)) % then in YC group
%             matched = find(YC_m == svec(s));
%         end
%         goodStim = overlapping{matched};
%         goodTrials = find(ismember(stimID,goodStim));
        cond = cell2mat(trials(:,9));
        rating = cell2mat(trials(:,12));
        easy = find(cond==2);
        hard = find(cond==1);
        rating_easy = rating(easy);
        rating_hard = rating(hard);
        %stimID = stimID(goodTrials,:);
        [~, horder] = sort(stimID(find(cond==1)));
        [~, eorder] = sort(stimID(find(cond==2)));
        easy_ordered(i,:) = rating_easy(eorder);
        hard_ordered(i,:) = rating_hard(horder);
        %hAvg(s,i) = nanmean(rating(hard));
        %eAvg(s,i) = nanmean(rating(easy));
        
    end
    all_easy(:,s) = easy_ordered(1,:)';
    all_hard(:,s) = hard_ordered(1,:)';
    diff_stim_hard(:,s) = hard_ordered(2,:)' - hard_ordered(1,:)';
    diff_easy(s) = nanmean(easy_ordered(2,:) - easy_ordered(1,:));
    diff_hard(s) = nanmean(hard_ordered(2,:) - hard_ordered(1,:));
    %clear easy_ordered hard_ordered
    %[allRem] = findRememberedStim(svec(s));
    if svec(s) >=13 %use diff filename
        r = dir(fullfile(behavioral_dir, ['EK25_' 'RECOG'  '*.mat']));
        r = load(fullfile(behavioral_dir,r(end).name));
    else
        r = dir(fullfile(behavioral_dir, ['_' 'RECOG'  '*.mat']));
        r = load(fullfile(behavioral_dir,r(end).name));
    end
    trials = table2cell(r.datastruct.trials);
    stimID = cell2mat(trials(:,8));
    %goodTrials = find(ismember(stimID,goodStim));
    cond = cell2mat(trials(:,9));
    [~, horder] = sort(stimID(find(cond==1)));
    acc = cell2mat(trials(:,11));
    rt = cell2mat(trials(:,13));
    easy = find(cond==2);
    hard = find(cond==1);
    hardacc_ordered(:,s) = acc(hard(horder));
    hardrt_ordered(:,s) = rt(hard(horder));
    easy_score(s) = nanmean(acc(easy));
    hard_score(s) = nanmean(acc(hard));
    easy_rt(s) = nanmedian(rt(easy));
    hard_rt(s) = nanmedian(rt(hard));
    
end

% diff_easy = easy_ordered(:,easy_rem{s},2) - easy_ordered(:,easy_rem{s},1);
% diff_hard = hard_ordered(:,hard_rem{s},2) - hard_ordered(:,hard_rem{s},1);
e_ALL = nanmean(diff_easy);
h_ALL = nanmean(diff_hard);

ee_ALL = nanstd(diff_easy)/sqrt(nsub-1);
eh_ALL = nanstd(diff_hard)/sqrt(nsub-1);

eALLD = [eh_ALL; ee_ALL];
ALLD = [h_ALL; e_ALL];
h = figure;
barwitherr(eALLD,ALLD)
set(gca,'XTickLabel' , ['RT  '; 'Omit']);
title('Average Post - Pre Difference')
xlabel('Stim Type')
ylabel('Level of Detail Difference')
%ylim([1 5.5])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)
%legend('Pre MOT', 'Post MOT')
%print(h, sprintf('%sYCONLY_ratings.pdf', allplotDir), '-dpdf')


e_ALL = nanmean(easy_score);
h_ALL = nanmean(hard_score);

ee_ALL = nanstd(easy_score)/sqrt(nsub-1);
eh_ALL = nanstd(hard_score)/sqrt(nsub-1);

eALLD = [eh_ALL; ee_ALL];
ALLD = [h_ALL; e_ALL];
h = figure;
barwitherr(eALLD,ALLD)
set(gca,'XTickLabel' , ['RT  '; 'Omit']);
title('Recog Acc')
xlabel('Stim Type')
ylabel('Accuracy')
%ylim([1 5.5])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)


e_ALL = nanmedian(easy_rt);
h_ALL = nanmedian(hard_rt);

ee_ALL = nanstd(easy_rt)/sqrt(nsub-1);
eh_ALL = nanstd(hard_rt)/sqrt(nsub-1);

eALLD = [eh_ALL; ee_ALL];
ALLD = [h_ALL; e_ALL];
h = figure;
barwitherr(eALLD,ALLD)
set(gca,'XTickLabel' , ['RT  '; 'Omit']);
title('RT')
xlabel('Stim Type')
ylabel('RT Median')
%ylim([1 5.5])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)



%% level of detail differences across groups
firstgroup = [diff_hard(iRT); diff_easy(iRT)];
secondgroup = [diff_hard(iYC); diff_easy(iYC)];
avgratio = [nanmean(firstgroup,2) nanmean(secondgroup,2)];
eavgratio = [nanstd(firstgroup,[],2)/sqrt(length(firstgroup)-1) nanstd(secondgroup,[],2)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['MOT ';'OMIT']);
legend('Realtime', 'Yoked')
xlabel('Stimulus Type')
ylabel('Level of Detail Difference')
title('Average Post - Pre Detail Difference')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylim([-.4 .8])
print(thisfig, sprintf('%salldetailRatings.pdf', allplotDir), '-dpdf')



cats = {'RT-MOT', 'YC-MOT', 'RT-OMIT', 'YC-OMIT'};
pl = {diff_hard(iRT), diff_hard(iYC), diff_easy(iRT), diff_easy(iYC)};
%clear mp;
%[~,mp] = ttest2(distDec2(iRT),distDec2(iYC));
%ps = [mp];
yl='Change in Imagery Detail'; %y-axis label
h = figure;plotSpread(pl,'xNames',cats,'showMM',2,'yLabel',yl); %this plots the beeswarm
h=gcf;set(h,'PaperOrientation','landscape'); %these two lines grabs some attributes important for plotting significance
%xt = get(gca, 'XTick');yt = get(gca, 'YTick');
%hold on;plotSig([1 2],yt,ps,0);hold off; %keep hold on and do plotSig.
%pn=[picd num2str(vers) '-' num2str(scramm) 'ChangeInPrecision'];%print(pn,'-depsc'); %print fig
%ylim([-1.25 1.25])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title('Post - Pre MOT Change in Detail');
print(h, sprintf('%sbeesdetails.pdf', allplotDir), '-dpdf')




%% recognition accuracy
firstgroup = [hard_score(iRT); easy_score(iRT)];
secondgroup = [hard_score(iYC); easy_score(iYC)];
avgratio = [nanmean(firstgroup,2) nanmean(secondgroup,2)];
eavgratio = [nanstd(firstgroup,[],2)/sqrt(length(firstgroup)-1) nanstd(secondgroup,[],2)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['MOT ';'OMIT']);
legend('Realtime', 'Yoked')
xlabel('Stimulus Type')
ylabel('Recognition Rate')
title('Recognition Accuracy')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylim([0 1])
%print(thisfig, sprintf('%sallrecogAcc.pdf', allplotDir), '-dpdf')
%% recognition RT
firstgroup = [hard_rt(iRT); easy_rt(iRT)];
secondgroup = [hard_rt(iYC); easy_rt(iYC)];
avgratio = [nanmedian(firstgroup,2) nanmedian(secondgroup,2)];
eavgratio = [nanstd(firstgroup,[],2)/sqrt(length(firstgroup)-1) nanstd(secondgroup,[],2)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['MOT ';'OMIT']);
legend('Realtime', 'Yoked')
xlabel('Stimulus Type')
ylabel('RT (s)')
title('Recognition RT')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylim([0 1])
%print(thisfig, sprintf('%sallrecogRT.pdf', allplotDir), '-dpdf')