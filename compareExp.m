% compare old and new feedback systems
% want to compare: distribution of evidence during feedback and second
% derivative of evidence during feedback
% also compare number of TR's spent in optimal zone!

oldDir = '/Data1/code/motStudy02/code/';
cd(oldDir)
projectName = 'motStudy02';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;
sepTRs = 17;
FBTRs = 11;
nblock = 3;
goodRange = [0 0.2]; % change back to 0.5 0.15 from wide

highB = 0.2;
lowB = highB * -1;
speedH = 1;
speedL = speedH * -1; %now change code for loop
avgRange = 12; %divided by 2 is the number of TR's

svec = [8 12 14 15 18 22 31];%take out 22 to make even by group];

nsub = length(svec);
sepbystimA = zeros(nstim,nTRs*3,nsub);


for s = 1:nsub
    subjectNum = svec(s);
    timecourse_high = [];
    timecourse_low = [];
    for iblock = 1:nblock
        blockNum = iblock;
        SESSION = 19 + blockNum;
        %blockNum = SESSION - 20 + 1;
        save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
        behavioral_dir = [fileparts(which('mot_realtime01.m')) '/BehavioralData/' num2str(subjectNum) '/'];
        runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
        fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime01_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
        names = {fileSpeed.name};
        dates = [fileSpeed.datenum];
        [~,newest] = max(dates);
        plotDir = ['/Data1/code/' projectName '/' 'Plots' '/' num2str(subjectNum) '/'];
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end
        matlabOpenFile = [behavioral_dir '/' names{newest}];
        d = load(matlabOpenFile);
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
        sepbytrial = reshape(categsep,nTRs,10); %right now TR x Trial
        sepbytrial = sepbytrial'; %results by trial number x TR number so 10 x 15
        [~,indSort] = sort(d.stim.id);
        sepinorder = sepbytrial(indSort,:); % so this is nTrials x 15 TR's
        sepbystimA(:,(iblock-1)*nTRs + 1: iblock*nTRs,s ) = sepinorder;
        FBsepbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = sepinorder(:,5:end);
        
        % now get all the speed changes
        allSpeed = d.stim.motionSpeed;
        speedbytrial = allSpeed'; %now it's trial x TR's
        speedinorder = speedbytrial(indSort,:);
       %subtract to get the difference in speed
        FBds = diff(speedinorder,[],2);
        FBdsbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = FBds(:,4:end);
        timecourse_high = nan(1,15);
        timecourse_low = nan(1,15);
        %TAKE OUT FOR SPEED
        % now see how timecourses change
%         for i = 1:nstim % because we're using all the stimuli now
%             tcourse = sepbytrial(i,:);
%             x1 = 1:length(tcourse);
%             x2 = 1:.5:length(tcourse);
%             y2 = interp1q(x1',tcourse',x2');
%             for j = 3:length(y2) - avgRange
%                 if y2(j-1) < highB && y2(j) > highB %then this is a POSITIVE CROSSING POINT
%                     timecourse_high(end+1,:) = y2(j-2:j+avgRange);
%                 elseif y2(j-1) > lowB && y2(j) < lowB %% then this is a NEGATIVE crossing point
%                     timecourse_low(end+1,:) = y2(j-2:j+avgRange);
%                 end
%             end
%         end
        
        %would be +-1.25 because that's what the setting for max increase
        %or decrease?
        for i = 1:nstim % because we're using all the stimuli now
            tcourse = sepbytrial(i,:);
            dscourse = FBds;
            x1 = 1:length(dscourse);
            x2 = 1:.5:length(dscourse);
            y2 = interp1q(x1',dscourse',x2');
            
            x1s = 1:length(tcourse);
            x2s = 1:.5:length(tcourse);
            y2s = interp1q(x1',tcourse',x2');
            for j = 3:length(y2) - avgRange -1
                if y2(j-1) < speedH && y2(j) > speedH %then this is a POSITIVE CROSSING POINT
                    timecourse_high(end+1,:) = y2s(j-2+1:j+avgRange+1);
                elseif y2(j-1) > speedL && y2(j) < speedL %% then this is a NEGATIVE crossing point
                    timecourse_low(end+1,:) = y2s(j-2+1:j+avgRange+1); %you go + 1 because of the index difference
                end
            end
        end
        
        % look how separation evidence is changing by TR
        allDiff = diff(diff(sepinorder(:,5:end),1,2),1,2)/4;
        secondDiff = reshape(allDiff,1,numel(allDiff));
        zTR = length(secondDiff);
        allSecondDiffA(s,(iblock-1)*zTR + 1: iblock*zTR) = secondDiff;
    end
    
    
    
    %take the average of nTRs in range
    %remove feedback and just look at all datapoints
    z1 =find(FBsepbystim>=goodRange(1));
    z2 = find(FBsepbystim<=goodRange(2));
    nGoodRangeA(s) = length(intersect(z1,z2))/numel(FBsepbystim);
    nConsecA(s) = sum(diff(intersect(z1,z2))==1)/numel(FBsepbystim);
    nLowA(s) = length(find(FBsepbystim<=goodRange(1)))/numel(FBsepbystim);
    nHighA(s) = length(find(FBsepbystim>=goodRange(2)))/numel(FBsepbystim);
    vectorSepA(s,:) = reshape(FBsepbystim,1,numel(FBsepbystim));
    if isnan(timecourse_high)
        avg_highA(s,:) = nan(1,15);
    else
        avg_highA(s,:) = nanmean(timecourse_high);
    end
    if isnan(timecourse_low)
        avg_lowA(s,:) = nan(1,15);
    else
        avg_lowA(s,:) = nanmean(timecourse_low);
    end
end
%% now do the same thing for new data
oldDir = '/Data1/code/motStudy03/code/';
cd(oldDir)
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
sepbystimB = zeros(nstim,nTRs*3,nsub);


for s = 1:nsub
    subjectNum = svec(s);
    timecourse_high = [];
    timecourse_low = [];
    for iblock = 1:nblock
        blockNum = iblock;
        SESSION = 20 + blockNum;
        %blockNum = SESSION - 20 + 1;
        save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
        behavioral_dir = [fileparts(which('mot_realtime02.m')) '/BehavioralData/' num2str(subjectNum) '/'];
        runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
        fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime02_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
        names = {fileSpeed.name};
        dates = [fileSpeed.datenum];
        [~,newest] = max(dates);
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end
        matlabOpenFile = [behavioral_dir '/' names{newest}];
        d = load(matlabOpenFile);
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
        sepbytrial = reshape(categsep,nTRs,10); %right now TR x Trial
        sepbytrial = sepbytrial'; %results by trial number x TR number so 10 x 15
        [~,indSort] = sort(d.stim.id);
        sepinorder = sepbytrial(indSort,:);
        sepbystimB(:,(iblock-1)*nTRs + 1: iblock*nTRs,s ) = sepinorder;
        FBsepbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = sepinorder(:,5:end);
        
         % now get all the speed changes
        allSpeed = d.stim.motionSpeed;
        speedbytrial = allSpeed'; %now it's trial x TR's
        speedinorder = speedbytrial(indSort,:);
       %subtract to get the difference in speed
        FBds = diff(speedinorder,[],2);
        FBdsbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = FBds(:,4:end);
        
        % look how separation evidence is changing by TR
        allDiff = diff(diff(sepinorder(:,5:end),1,2),1,2)/4;
        secondDiff = reshape(allDiff,1,numel(allDiff));
        zTR = length(secondDiff);
        allSecondDiffB(s,(iblock-1)*zTR + 1: iblock*zTR) = secondDiff;
        % now see how timecourses change
        
        %take this out if going from SPEED and NOT EVIDENCE!
%         for i = 1:nstim % because we're using all the stimuli now
%             tcourse = sepbytrial(i,:);
%             x1 = 1:length(tcourse);
%             x2 = 1:.5:length(tcourse);
%             y2 = interp1q(x1',tcourse',x2');
%             for j = 3:length(y2) - avgRange
%                 if y2(j-1) < highB && y2(j) > highB %then this is a POSITIVE CROSSING POINT
%                     timecourse_high(end+1,:) = y2(j-2:j+avgRange);
%                 elseif y2(j-1) > lowB && y2(j) < lowB %% then this is a NEGATIVE crossing point
%                     timecourse_low(end+1,:) = y2(j-2:j+avgRange);
%                 end
%             end
%         end
    end
    
       %would be +-1.25 because that's what the setting for max increase
        %or decrease?
        timecourse_high = nan(1,15);
        timecourse_low = nan(1,15);
        for i = 1:nstim % because we're using all the stimuli now
            tcourse = sepbytrial(i,:);
            dscourse = FBds;
            x1 = 1:length(dscourse);
            x2 = 1:.5:length(dscourse);
            y2 = interp1q(x1',dscourse',x2');
            
            x1s = 1:length(tcourse);
            x2s = 1:.5:length(tcourse);
            y2s = interp1q(x1',tcourse',x2');
            for j = 3:length(y2) - avgRange -1
                if y2(j-1) < speedH && y2(j) > speedH %then this is a POSITIVE CROSSING POINT
                    timecourse_high(end+1,:) = y2s(j-2+1:j+avgRange+1);
                elseif y2(j-1) > speedL && y2(j) < speedL %% then this is a NEGATIVE crossing point
                    timecourse_low(end+1,:) = y2s(j-2+1:j+avgRange+1); %you go + 1 because of the index difference
                end
            end
        end
    
    %take the average of nTRs in range
    %remove feedback and just look at all datapoints
    z1 =find(FBsepbystim>=goodRange(1));
    z2 = find(FBsepbystim<=goodRange(2));
    nGoodRangeB(s) = length(intersect(z1,z2))/numel(FBsepbystim);
    nConsecB(s) = sum(diff(intersect(z1,z2))==1)/numel(FBsepbystim);
    nLowB(s) = length(find(FBsepbystim<=goodRange(1)))/numel(FBsepbystim);
    nHighB(s) = length(find(FBsepbystim>=goodRange(2)))/numel(FBsepbystim);
    vectorSepB(s,:) = reshape(FBsepbystim,1,numel(FBsepbystim));
     if isnan(timecourse_high)
        avg_highB(s,:) = nan(1,15);
    else
        avg_highB(s,:) = nanmean(timecourse_high);
    end
    if isnan(timecourse_low)
        avg_lowB(s,:) = nan(1,15);
    else
        avg_lowB(s,:) = nanmean(timecourse_low);
    end
    
end
%% now save to plot in python
folder= '/jukebox/norman/amennen/PythonMot3';
save('compareExp.mat','nGoodRangeA', 'nGoodRangeB','avg_highA', 'avg_highB', 'avg_lowA', 'avg_lowB', 'nLowA', 'nLowB', 'nHighA', 'nHighB','vectorSepA', 'vectorSepB', 'allSecondDiffA', 'allSecondDiffB','nConsecA','nConsecB','sepbystimA', 'sepbystimB');
unix(['scp ' 'compareExp.mat' ' amennen@apps.pni.princeton.edu:' folder '/' ])
