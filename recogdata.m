%recognition task
recog = 24; %associates task
projectName = 'motStudy02';
allplotDir = ['/Data1/code/' projectName '/' 'Plots' '/' ];
svec = [8 12:15 18 21 22];
RTgroup = 0;
YCgroup = 1;
if RTgroup
svec = [8 12:15 18 21 22];
elseif YCgroup
    svec = [16 20];
end
NSUB = length(svec);

for s = 1:NSUB
    [allRem] = findRememberedStim(svec(s));
    
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(svec(s)) '/'];
        r = dir(fullfile(behavioral_dir, ['_' 'RECOG'  '*.mat'])); 
        r = load(fullfile(behavioral_dir,r(end).name)); 
        trials = table2cell(r.datastruct.trials);
        stimID = cell2mat(trials(:,8));
        goodTrials = find(ismember(stimID,allRem));
        cond = cell2mat(trials(goodTrials,9));
        acc = cell2mat(trials(goodTrials,11));
        rt = cell2mat(trials(goodTrials,13));
        easy = find(cond==2);
        hard = find(cond==1);
        easy_score(s) = mean(acc(easy));
        hard_score(s) = mean(acc(hard));
        easy_rt(s) = nanmedian(rt(easy));
        hard_rt(s) = nanmedian(rt(hard));
end
eALLD = [std(hard_score)/sqrt(NSUB-1);std(easy_score)/sqrt(NSUB-1)];
ALLD = [mean(hard_score); mean(easy_score)];
h = figure;
barwitherr(eALLD,ALLD)
set(gca,'XTickLabel' , ['RT  '; 'Omit']);
title('Recognition Accuracy')
xlabel('Stim Type')
ylabel('Recognition Rate (%)')
%ylim([1 5.5])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)
%print(h, sprintf('%srecogacc.pdf', allplotDir), '-dpdf')

eALLD = [nanstd(hard_rt)/sqrt(NSUB-1);nanstd(easy_rt)/sqrt(NSUB-1)];
ALLD = [nanmedian(hard_rt); nanmedian(easy_rt)];
h = figure;
barwitherr(eALLD,ALLD)
set(gca,'XTickLabel' , ['RT  '; 'Omit']);
title('Recognition RT')
xlabel('Stim Type')
ylabel('RT (s)')
%ylim([1 5.5])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)
%print(h, sprintf('%srecogrt.pdf', allplotDir), '-dpdf')
