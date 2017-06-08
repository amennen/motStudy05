% compare feedback
% run compareFeedback in motStudy02 folder first!
%close all;
%clear all;
projectName = 'motStudy03';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;
sepTRs = 17;
FBTRs = 11;
nblock = 3;
highB = 0.15;
lowB = -.15;
avgRange = 12; % divided by 2 is Number of TR's
svec = [3:9 11 12 13];
onlyRem = 1;
RT = svec;

nsub = length(svec);
sepbystim = zeros(nstim,nTRs*3);
speedbystim = zeros(nstim,nTRs*3);
allplotDir = ['/Data1/code/' projectName '/' 'Plots' '/' ];

for s = 1:nsub
    timecourse_high = [];
    timecourse_low = [];
    
    subjectNum = svec(s);
    remStim = findRememberedStim(subjectNum);
   
    for iblock = 1:nblock
        blockNum = iblock;
        SESSION = 20 + blockNum;
        
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
        if onlyRem
            goodTrials = find(ismember(d.stim.id,remStim));
        else
            goodTrials = 1:nstim;
        end
        
        nTotalTR = 120; %120 fb TR's per block        
        allSpeed = d.stim.motionSpeed; %matrix of TR's
        if onlyRem
            allSpeed = allSpeed(:,goodTrials);
        end
        speedVector = reshape(allSpeed,1,numel(allSpeed));
        allMotionTRs = convertTR(d.timing.trig.wait,d.timing.plannedOnsets.motion,d.config.TR); %row,col = mTR,trialnumber
        allMotionTRs = allMotionTRs + 2;%[allMotionTRs; allMotionTRs(end,:)+1; allMotionTRs(end,:) + 2]; %add in the next 2 TR's for HDF
        onlyFbTRs = allMotionTRs(4:end,:);
        FBTR2 = allMotionTRs(5:end,end);
        TRvector = reshape(allMotionTRs,1,numel(allMotionTRs));
        FBTRVector = reshape(onlyFbTRs,1,numel(onlyFbTRs));
        FBTRVector2 = reshape(FBTR2,1,numel(FBTR2));
        run = dir([runHeader 'motpatternsdata_' num2str(SESSION) '*']);
        names = {run.name};
        dates = [run.datenum];
        [~,newest] = max(dates);
        run = load(fullfile(runHeader,run(end).name));
        categsep = run.patterns.categsep(TRvector - 10); %minus 10 because we take out those 10
        sepbytrial = reshape(categsep,nTRs,10);
        if onlyRem
            sepbytrial = sepbytrial(:,goodTrials);
        end
        allsepchange = diff(sepbytrial,1,1);
        FBsepchange = reshape(allsepchange(4:end,:),1,numel(allsepchange(4:end,:)));
        allsep = reshape(sepbytrial(5:end,:),1,numel(sepbytrial(5:end,:)));
        allspeedchanges = diff(d.stim.motionSpeed,1,1);
        if onlyRem
            allspeedchanges = allspeedchanges(:,goodTrials);
        end
        FBspeed = reshape(allSpeed(5:end,:),1,numel(allSpeed(5:end,:)));
        FBspeedchange = reshape(allspeedchanges(4:end,:),1,numel(allspeedchanges(4:end,:)));
        FBTRs = length(FBspeedchange);
        secondDspeed = diff(allspeedchanges,1,1);
        d2speed = reshape(secondDspeed(3:end,:),1,numel(secondDspeed(3:end,:)));
      
        
        %look get avg time when above .1
        %first try only dot tracking conditions
        for i = 1:goodTrials
            tcourse = sepbytrial(:,i);
            x1 = 1:length(tcourse);
            x2 = 1:.5:length(tcourse);
            y2 = interp1q(x1',tcourse,x2');
            for j = 3:length(y2) - avgRange
                if y2(j-1) < highB && y2(j) > highB %then this is a POSITIVE CROSSING POINT
                    timecourse_high(end+1,:) = y2(j-2:j+avgRange);
                elseif y2(j-1) > lowB && y2(j) < lowB %% then this is a NEGATIVE crossing point
                    timecourse_low(end+1,:) = y2(j-2:j+avgRange);
                end
            end
        end
        
    alld2speed((iblock-1)*(FBTRs) + 1: iblock*(FBTRs) ,s) = d2speed;
    end
    avg_high(s,:) = nanmean(timecourse_high);
    avg_low(s,:) = nanmean(timecourse_low);
end
%% prove that more over and undershooting is because of feedback (relate previous dot speed to evidence max or min)

% between one peak and the next min-- ask what the change in speed was
% between those points?
timeHigh = [ nanmean(avg_high); nanmean(RTcompare_high)];
nsub2 = size(RTcompare_high,1);
eHigh = [nanstd(avg_high,[],1)/sqrt(nsub-1) ;nanstd(RTcompare_high,[],1)/sqrt(nsub2-1) ]; 
npts = size(avg_high,2);
figure;
mseb(1:npts,timeHigh, eHigh);
title(sprintf('High Timecourse'))
xlim([1 npts])

set(gca, 'XTick', [1:npts])
set(gca,'XTickLabel',['-2';'-1'; ' 0'; ' 1'; ' 2'; ' 3'; ' 4'; ' 5'; ' 6'; ' 7'; ' 8'; ' 9'; '10'; '11'; '12']);
ylabel('Retrieval - Control Evidence')
xlabel('Time Points')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('Real-time', 'Orig RT')
%print(h, sprintf('%sHighTimecourseUPDATED.pdf', allplotDir), '-dpdf') 
%% do the same for low points
timelow = [ nanmean(avg_low); nanmean(RTcompare_low)];
elow = [nanstd(avg_low,[],1)/sqrt(nsub-1);nanstd(RTcompare_high,[],1)/sqrt(nsub2-1) ] ;
npts = size(avg_low,2);

h = figure;
mseb(1:npts,timelow, elow);
xlim([1 npts])
title(sprintf('Low Timecourse'))
set(gca, 'XTick', [1:npts])
%set(gca,'XTickLabel',['0'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8']);
set(gca,'XTickLabel',['-2';'-1'; ' 0'; ' 1'; ' 2'; ' 3'; ' 4'; ' 5'; ' 6'; ' 7'; ' 8'; ' 9'; '10'; '11'; '12']);

ylabel('Retrieval - Control Evidence')
xlabel('Time Points')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('Real-time', 'Orig RT')
%print(h, sprintf('%sLowTimecourseUDPATEDp25.pdf', allplotDir), '-dpdf') 


%% next check that this holds for only remembered stimuli (maybe just for each person)
% and then also plot second derivative
