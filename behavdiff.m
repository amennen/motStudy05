        

close all;
clear all;
projectName = 'motStudy02';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;
nblock = 3;
svec = [3:5 7];
nsub = length(svec);
sepbystim = zeros(nstim,nTRs*3);
speedbystim = zeros(nstim,nTRs*3);

for s = 1:nsub
    subjectNum = svec(s);
behavioral_dir = [fileparts(which('mot_realtime01.m')) '/BehavioralData/' num2str(subjectNum) '/'];

    recallFile = dir(fullfile(behavioral_dir, ['EK19_SUB' '*mat']));
    r1 = load([behavioral_dir '/' recallFile(end).name]);
    z = table2cell(r1.datastruct.trials(:,16));
    z(cellfun(@(x) any(isnan(x)),z)) = {'00'};

    resp1 = cell2mat(z);
    resp1 = resp1(:,1);
    
    stimOrder = table2array(r1.datastruct.trials(:,8));
    RTorder = stimOrder(stimOrder<11);
    RTonly = resp1(stimOrder<11);
    [~,sortedID] = sort(RTorder);
    r1Sort = RTonly(sortedID);
    
     %now all recall changes
    [~,allsort] = sort(stimOrder);
    R1recall = resp1(allsort);
    
    recallFile = dir(fullfile(behavioral_dir, ['EK23_SUB' '*mat']));
    r2 = load([behavioral_dir '/' recallFile(end).name]);
    z = table2cell(r2.datastruct.trials(:,16));
    z(cellfun(@(x) any(isnan(x)),z)) = {'00'}; %for nan's!
    resp2 = cell2mat(z);
    resp2 = resp2(:,1);
    
    stimOrder = table2array(r2.datastruct.trials(:,8));
    RTorder = stimOrder(stimOrder<11);
    RTonly = resp2(stimOrder<11);
    [~,sortedID] = sort(RTorder);
    r2Sort = RTonly(sortedID);
    
     %now all recall changes
    [~,allsort] = sort(stimOrder);
    R2recall = resp2(allsort);
    
    recalldiff = R2recall - R1recall;
   rtdiff(s) = mean(recalldiff(1:10)); 
   omitdiff(s) = mean(recalldiff(11:end));
end

allavg = mean([rtdiff;omitdiff],2);
allstd = std([rtdiff;omitdiff], [],2)/sqrt(nsub-1);
figure;
h = barwitherr(allstd,allavg)
set(gca,'XTickLabel' , ['Realtime'; '  Omit  ']);
title('Detail Change Post-Pre MOT')
xlabel('Stimulus Category')
ylabel('Detail Change Post-Pre MOT')
set(h,'FaceColor', [73 91 252]/255);
%ylim([1 5.5])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)
%legend('Realtime', 'Omit')


        