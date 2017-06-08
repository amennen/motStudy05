% here: going to quantify how long it takes to drop down classifier
% evidence
% can compare between the two versions and with realtime and non realtime
% versions

% need: speed (from behavioral file)
% category separation (can also pull from behavioral file)
close all;
clear all;thisDir = '/Data1/code/motStudy03/code/';
cd(thisDir)
projectName = 'motStudy03';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;
sepTRs = 17;
FBTRs = 11;
nblock = 3;
svec = [3:9 11 12 13];

nsub = length(svec);
sepbystim = zeros(nstim,nTRs*3);
speedbystim = zeros(nstim,nTRs*3);
allplotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' ];
iRT = 1:nsub;
for s = 1:nsub
    subjectNum = svec(s);
    for iblock = 1:nblock
        blockNum = iblock;
        SESSION = 20 + blockNum;
        %blockNum = SESSION - 20 + 1;
        
        %behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
        behavioral_dir = [fileparts(which('mot_realtime02.m')) '/BehavioralData/' num2str(subjectNum) '/'];
        save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
        runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
        fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime02_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
        names = {fileSpeed.name};
        dates = [fileSpeed.datenum];
        [~,newest] = max(dates);
        plotDir = ['/Data1/code/' projectName '/' 'Plots' '/' num2str(subjectNum) '/'];
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end
        matlabOpenFile = [behavioral_dir '/' names{newest}];
        d = load(matlabOpenFile);
        allSpeed = d.stim.motionSpeed; %matrix of TR's
        speedVector = reshape(allSpeed,1,numel(allSpeed));
        allMotionTRs = convertTR(d.timing.trig.wait,d.timing.plannedOnsets.motion,d.config.TR); %row,col = mTR,trialnumber
        allMotionTRs = allMotionTRs + 2;%[allMotionTRs; allMotionTRs(end,:)+1; allMotionTRs(end,:) + 2]; %add in the next 2 TR's for HDF
        onlyFbTRs = allMotionTRs(5:end,:);
        TRvector = reshape(allMotionTRs,1,numel(allMotionTRs));
        FBTRVector = reshape(onlyFbTRs,1,numel(onlyFbTRs));
        run = dir([runHeader 'motpatternsdata_' num2str(SESSION) '*']);
        names = {run.name};
        dates = [run.datenum];
        [~,newest] = max(dates);
        run = load(fullfile(runHeader,run(end).name));
        categsep = run.patterns.categsep(TRvector - 10); %minus 10 because we take out those 10
        sepbytrial = reshape(categsep,nTRs,10);
        sepbytrial = sepbytrial'; %results by trial number, TR number
        %want to match speed changes so change back to 15 TR's
        %sepbytrial = sepbytrial(:,1:nTRs);
        speedbytrial = reshape(speedVector,nTRs,nstim);
        speedbytrial = speedbytrial';
        [~,indSort] = sort(d.stim.id);
        sepinorder = sepbytrial(indSort,:);
        speedinorder = speedbytrial(indSort,:);
        sepbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = sepinorder;
        FBsepbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = sepinorder(:,5:end);
        speedbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = speedinorder;
        FBspeedbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = speedinorder(:,5:end);
        allspeedchanges = d.stim.changeSpeed;
        speedchange = allspeedchanges'; % now in trial, TR order
        speedchangeinorder = speedchange(indSort,:);
        dsbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = speedchangeinorder;
        meandiffbyblock(iblock) = mean(mean(abs(diff(sepbystim,1,2))));
        FBmeandiffbyblock(iblock) = mean(mean(abs(diff(FBsepbystim,1,2))));
        
        % look how separation evidence is changing by TR
        allDiff = diff(diff(sepinorder(:,5:end),1,2),1,2)/4;
        secondDiff = reshape(allDiff,1,numel(allDiff));
        zTR = length(secondDiff);
        allSecondDiff(s,(iblock-1)*zTR + 1: iblock*zTR) = secondDiff;
    end
    meandiff = mean(meandiffbyblock);
    FBmeandiff = mean(FBmeandiffbyblock);
  
    postime = [];
    negtime = [];
    zerotime = [];
    posinc = [];
    negdec = [];
%     for i = 1:nstim
%         %thisSep = smoothedsep(i,:); %can choose to look at smoothed
%         %version or raw separation
%         thisSep = sepbystim(i,:);
%         thisdS = dsbystim(i,:);
%         thisSpeed = speedbystim(i,:);
%         highPts = find(thisSep>0.15);
%         %acceptable = intersect(find(mod(1:45,15)~=1),find(mod(1:45,15)~=2)); % don't take if in the first 2 TR's because weren't looking at anything 
%         %keep = find(ismember(highPts,acceptable));
%         %highPts = highPts(keep);
%         dist = nan(1,length(highPts));
%         clear avgds;
%         for h = 1:length(highPts)
%             newvals = thisSep(highPts(h)+1:end);
%             %nextdown = find(newvals<0);
%             nextdown = find(thisSep(highPts(h)) - newvals >= meandiff*2 );
%             if ~isempty(nextdown)
%                 dist(h) = nextdown(1);
%                 %avgds(h) = sum(thisdS(highPts(h)+1:highPts(h)+nextdown(1)));
%                 avgds(h) = thisSpeed(highPts(h)+nextdown(1)) - thisSpeed(highPts(h));
%                 if avgds(h) > 0
%                     postime = [postime dist(h)];
%                     posinc = [posinc avgds(h)];
%                 elseif avgds(h) < 0
%                     negtime = [negtime dist(h)];
%                     negdec = [negdec avgds(h)];
%                 end
%             end
%         end
%         distDec1(s,i) = nanmean(dist); 
%     end
    
    % now do for FB TRs only
    for i = 1:nstim
        %thisSep = smoothedsep(i,:); %can choose to look at smoothed
        %version or raw separation
        thisSep = FBsepbystim(i,:);
        thisSpeed = FBspeedbystim(i,:);
        highPts = find(thisSep>0.15);
        %acceptable = intersect(find(mod(1:45,15)~=1),find(mod(1:45,15)~=2)); % don't take if in the first 2 TR's because weren't looking at anything 
        %keep = find(ismember(highPts,acceptable));
        %highPts = highPts(keep);
        dist = nan(1,length(highPts));
        clear avgds;
        for h = 1:length(highPts)
            newvals = thisSep(highPts(h)+1:end);
            %nextdown = find(newvals<0);
            nextdown = find(thisSep(highPts(h)) - newvals >= FBmeandiff );
            if ~isempty(nextdown)
                actualIndex = highPts(h)+nextdown(1);
                if (~mod(actualIndex,FBTRs) || mod(actualIndex,FBTRs) > mod(highPts(h),FBTRs)) && (mod(highPts(h),FBTRs)) %only accept if in the same block
                    dist(h) = nextdown(1);
                    %avgds(h) = sum(thisdS(highPts(h)+1:highPts(h)+nextdown(1)));
                    avgds(h) = thisSpeed(highPts(h)+nextdown(1)) - thisSpeed(highPts(h));
                    if avgds(h) > 0
                        postime = [postime dist(h)];
                        posinc = [posinc avgds(h)];
                    elseif avgds(h) < 0
                        negtime = [negtime dist(h)];
                        negdec = [negdec avgds(h)];
                    end
                end
            end
        end
        distDec1(s,i) = nanmean(dist); 
    end
    
    npos(s) = length(postime);
    nneg(s) = length(negtime);
    posavg(s) = mean(posinc);
    negavg(s) = mean(negdec);
    distPos(s) = mean(postime);
    distNeg(s) = mean(negtime);
    distDec2 = nanmean(distDec1,2);
    
    %take the average of nTRs in range
    goodRange = [0.05 0.15];
    %remove feedback and just look at all datapoints
    z1 =find(sepbystim>=goodRange(1));
    z2 = find(sepbystim<=goodRange(2));
    nGoodRange(s) = length(intersect(z1,z2))/numel(sepbystim);
    nLow(s) = length(find(sepbystim<=-.1))/numel(sepbystim);
    nHigh(s) = length(find(sepbystim>=.3))/numel(sepbystim);
    allSepMean(s) = mean(mean(sepbystim));
    vectorSep(s,:) = reshape(sepbystim,1,numel(sepbystim));
end

%% now compare decrease of retrieval evidence across both groups
%if (~updated && ~oldonly)
    cats = {'RT', 'YC'};
    pl = {distDec2(iRT), distDec2(iYC)};
    clear mp;
    [~,mp] = ttest2(distDec2(iRT),distDec2(iYC));
    ps = [mp];
    yl='nTRs to Decrease Evidence'; %y-axis label
    h = figure;plotSpread(pl,'xNames',cats,'showMM',2,'yLabel',yl); %this plots the beeswarm
    h=gcf;set(h,'PaperOrientation','landscape'); %these two lines grabs some attributes important for plotting significance
    xt = get(gca, 'XTick');yt = get(gca, 'YTick');
    hold on;plotSig([1 2],yt,ps,0);hold off; %keep hold on and do plotSig.
    %pn=[picd num2str(vers) '-' num2str(scramm) 'ChangeInPrecision'];%print(pn,'-depsc'); %print fig
    %ylim([-1.25 1.25])
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    title('Time to Decrease Evidence by Group');
    print(h, sprintf('%sbeesTIMETODEC.pdf', allplotDir), '-dpdf')
    
    
    firstgroup = distDec2(iRT);
    secondgroup = distDec2(iYC);
    avgratio = [mean(firstgroup) mean(secondgroup)];
    eavgratio = [std(firstgroup)/sqrt(length(firstgroup)-1) std(secondgroup)/sqrt(length(secondgroup)-1)];
    thisfig = figure;
    barwitherr(eavgratio,avgratio)
    set(gca,'XTickLabel' , ['RT';'YC']);
    xlabel('Subject Group')
    ylabel('TR''s to Decrease')
    title('Time to Decrease Evidence by Group')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    %ylim([-.2 0.2])
    print(thisfig, sprintf('%sTIMETODEC.pdf', allplotDir), '-dpdf')
%end

%% now compare time to decrease with types of feedback
%     firstPos = distPos(iRT);
%     secondPos = distPos(iYC);
%     firstNeg = distNeg(iRT);
%     secondNeg = distNeg(iYC);
%     
%     avgratio = [mean(firstPos)  mean(secondPos); mean(firstNeg) mean(secondNeg)];
%     eavgratio = [std(firstPos)/sqrt(length(firstPos-1)) std(secondPos)/sqrt(length(secondPos-1));std(firstNeg)/sqrt(length(firstPos-1)) std(secondNeg)/sqrt(length(secondPos-1))];
%     thisfig = figure;
%     barwitherr(eavgratio,avgratio)
%     set(gca,'XTickLabel' , ['dS > 0';'dS < 0']);
%     legend('Old 4', 'New 6')
%     ylabel('TR''s to Decrease')
%     title('Evidence Response Time, Separated by \DeltaS')
%     set(findall(gcf,'-property','FontSize'),'FontSize',20)
%     %print(thisfig, sprintf('%sweakeningbygroupbysign.pdf', allplotDir), '-dpdf')

%% and then do time spent in optimal weakening zone-just fb TR's only
cats = {'RT', 'YC'};
pl = {nGoodRange(iRT), nGoodRange(iYC)};
clear mp;
[~,mp] = ttest2(nGoodRange(iRT),nGoodRange(iYC));
ps = [mp];
yl='Ratio in Good Range'; %y-axis label
h = figure;plotSpread(pl,'xNames',cats,'showMM',2,'yLabel',yl); %this plots the beeswarm
h=gcf;set(h,'PaperOrientation','landscape'); %these two lines grabs some attributes important for plotting significance
xt = get(gca, 'XTick');yt = get(gca, 'YTick');
hold on;plotSig([1 2],yt,ps,0);hold off; %keep hold on and do plotSig.
%pn=[picd num2str(vers) '-' num2str(scramm) 'ChangeInPrecision'];%print(pn,'-depsc'); %print fig
%ylim([-1.25 1.25])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title('Ratio Spent in Optimal Range');
print(h, sprintf('%sbeesTIMEINOPTIMAL.pdf', allplotDir), '-dpdf')


firstgroup = nGoodRange(iRT);
secondgroup = nGoodRange(iYC);
avgratio = [mean(firstgroup) mean(secondgroup)];
eavgratio = [std(firstgroup)/sqrt(length(firstgroup)-1) std(secondgroup)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['RT';'YC']);
xlabel('Subject Group')
ylabel('% FB Time in Optimal Range')
title('Time Spent in Optimal Range')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%ylim([-.2 0.2])
print(thisfig, sprintf('%sTIMEINOPTIMAL.pdf', allplotDir), '-dpdf')
%% and then do time spent BELOW optimal weakening zone-just fb TR's only
cats = {'RT', 'YC'};
pl = {nLow(iRT), nLow(iYC)};
clear mp;
[~,mp] = ttest2(nLow(iRT),nLow(iYC));
ps = [mp];
yl='Ratio in Low Range'; %y-axis label
h = figure;plotSpread(pl,'xNames',cats,'showMM',2,'yLabel',yl); %this plots the beeswarm
h=gcf;set(h,'PaperOrientation','landscape'); %these two lines grabs some attributes important for plotting significance
xt = get(gca, 'XTick');yt = get(gca, 'YTick');
hold on;plotSig([1 2],yt,ps,0);hold off; %keep hold on and do plotSig.
%pn=[picd num2str(vers) '-' num2str(scramm) 'ChangeInPrecision'];%print(pn,'-depsc'); %print fig
%ylim([-1.25 1.25])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title('Ratio Spent in Low Range');
print(h, sprintf('%sbeesTIMELOW.pdf', allplotDir), '-dpdf')


firstgroup = nLow(iRT);
secondgroup = nLow(iYC);
avgratio = [mean(firstgroup) mean(secondgroup)];
eavgratio = [std(firstgroup)/sqrt(length(firstgroup)-1) std(secondgroup)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['RT';'YC']);
xlabel('Subject Group')
ylabel('% FB Time in Low Range')
title('Time Spent in Low')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%ylim([-.2 0.2])
%print(thisfig, sprintf('%sTIMELOW.pdf', allplotDir), '-dpdf')
%% and then do time spent in HIGH zone-just fb TR's only
cats = {'RT', 'YC'};
pl = {nHigh(iRT), nHigh(iYC)};
clear mp;
[~,mp] = ttest2(nHigh(iRT),nHigh(iYC));
ps = [mp];
yl='Ratio in High Range'; %y-axis label
h = figure;plotSpread(pl,'xNames',cats,'showMM',2,'yLabel',yl); %this plots the beeswarm
h=gcf;set(h,'PaperOrientation','landscape'); %these two lines grabs some attributes important for plotting significance
xt = get(gca, 'XTick');yt = get(gca, 'YTick');
hold on;plotSig([1 2],yt,ps,0);hold off; %keep hold on and do plotSig.
%pn=[picd num2str(vers) '-' num2str(scramm) 'ChangeInPrecision'];%print(pn,'-depsc'); %print fig
%ylim([-1.25 1.25])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title('Ratio Spent in High Range');
print(h, sprintf('%sbeesTIMEHIGH.pdf', allplotDir), '-dpdf')



firstgroup = nHigh(iRT);
secondgroup = nHigh(iYC);
avgratio = [mean(firstgroup) mean(secondgroup)];
eavgratio = [std(firstgroup)/sqrt(length(firstgroup)-1) std(secondgroup)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['RT';'YC']);
xlabel('Subject Group')
ylabel('% FB Time in High Range')
title('Time Spent in High Range')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%ylim([-.2 0.2])
%print(thisfig, sprintf('%sTIMEHIGH.pdf', allplotDir), '-dpdf')

%% mean evidence of each group
firstgroup =allSepMean(iRT);
secondgroup = allSepMean(iYC);
avgratio = [mean(firstgroup) mean(secondgroup)];
eavgratio = [std(firstgroup)/sqrt(length(firstgroup)-1) std(secondgroup)/sqrt(length(secondgroup)-1)];
thisfig = figure;
barwitherr(eavgratio,avgratio)
set(gca,'XTickLabel' , ['RT';'YC']);
xlabel('Subject Group')
ylabel('Mean of Classifier Evidence')
title('Classifier Evidence of Retrieval during FB')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
print(thisfig, sprintf('%sMEANEVIDENCE.pdf', allplotDir), '-dpdf')

%% histograms of evidence!
allRT = vectorSep(iRT,:);
RTvec = reshape(allRT,1,numel(allRT));
allYC = vectorSep(iYC,:);
YCvec = reshape(allYC,1,numel(allYC));

[c1,b1] = hist(RTvec,[-1:.1:1]);
hold on;
[c2,b] = hist(YCvec,b1);
RTnorm = c1/sum(c1);
YCnorm = c2/sum(c2);

figure;
bar(b,[RTnorm' YCnorm']);
ylim([0 .3]);
xlim([-1 1]);
legend('RT', 'YC')

%% histograms of second deriv
cats = {'RT', 'YC'};
SDRT = abs(reshape(allSecondDiff(iRT,:),numel(allSecondDiff(iRT,:)),1));
SDYC = abs(reshape(allSecondDiff(iYC,:),numel(allSecondDiff(iYC,:)),1));
pl = {SDRT, SDYC};
clear mp;
[~,mp] = ttest2(SDRT,SDYC);
ps = [mp];
yl='Second Derivative'; %y-axis label
h = figure;plotSpread(pl,'xNames',cats,'showMM',2,'yLabel',yl); %this plots the beeswarm
h=gcf;set(h,'PaperOrientation','landscape'); %these two lines grabs some attributes important for plotting significance
xt = get(gca, 'XTick');yt = get(gca, 'YTick');
hold on;plotSig([1 2],yt,ps,0);hold off; %keep hold on and do plotSig.
%pn=[picd num2str(vers) '-' num2str(scramm) 'ChangeInPrecision'];%print(pn,'-depsc'); %print fig
%ylim([-1.25 1.25])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
title('Second Derivative by Group');
print(h, sprintf('%sbeesSECONDD.pdf', allplotDir), '-dpdf')

thisfig = figure;
distributionPlot(pl, 'showMM', 2, 'xNames', cats, 'ylabel', yl, 'colormap', copper)
xlim([.5 2.5])
ylim([-.05 .4])
title('Second Difference of Evidence During MOT')
xlabel('Subject Group')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
print(thisfig, sprintf('%sviolinsSECONDD.pdf', allplotDir), '-dpdf')


secondRT = allSecondDiff(iRT,:);
secondYC = allSecondDiff(iYC,:);
RTvec = reshape(secondRT,1,numel(secondRT));
YCvec = reshape(secondYC,1,numel(secondYC));

[c1,b1] = hist(RTvec,[-1.5:.15:1.5]);
hold on;
[c2,b] = hist(YCvec,b1);
RTnorm = c1/sum(c1);
YCnorm = c2/sum(c2);

figure;
bar(b,[RTnorm' YCnorm']);
%ylim([0 .3]);
xlim([-1.5 1.5]);
legend('RT', 'YC')


