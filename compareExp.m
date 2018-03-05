% compare old and new feedback systems
% want to compare: distribution of evidence during feedback and second
% derivative of evidence during feedback
% also compare number of TR's spent in optimal zone!

oldDir = '/Data1/code/motStudy05/code/';
cd(oldDir)
projectName = 'motStudy05';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;
sepTRs = 17;
FBTRs = 11;
nblock = 3;
goodRange = [0.05 0.25]; % change back to 0.5 0.15 from wide
optimal = 0.15;
highB = optimal + .05;
lowB = optimal -.05;
speedH = 1;
speedL = speedH * -1; %now change code for loop
avgRange = 12; %divided by 2 is the number of TR's

svec = [1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40];

nsub = length(svec);
sepbystimD = zeros(nstim,nTRs*3,nsub);
speedbystimD = zeros(nstim,nTRs*3,nsub);


for s = 1:nsub
    subjectNum = svec(s);
    fprintf('now on subject %i\n', s);
    timecourse_high = [];
    timecourse_low = [];
    for iblock = 1:nblock
        blockNum = iblock;
        fprintf('now on block %i\n', iblock);
        SESSION = 20 + blockNum;
        %blockNum = SESSION - 20 + 1;
        save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
        behavioral_dir = [fileparts(which('mot_realtime05.m')) '/BehavioralData/' num2str(subjectNum) '/'];
        runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
        fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime05_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
        names = {fileSpeed.name};
        dates = [fileSpeed.datenum];
        [~,newest] = max(dates);

        matlabOpenFile = [behavioral_dir '/' names{newest}];
        d = load(matlabOpenFile);
        allMotionTRs = convertTR(d.timing.trig.wait,d.timing.plannedOnsets.motion,d.config.TR); %row,col = mTR,trialnumber
        allMotionTRs = allMotionTRs + 2;%[allMotionTRs; allMotionTRs(end,:)+1; allMotionTRs(end,:) + 2]; %add in the next 2 TR's for HDF
        TRvector = reshape(allMotionTRs,1,numel(allMotionTRs));
        fname = findOldestFile(runHeader,fullfile(runHeader,['motpatterns'  '*.mat']));
        run = load(fname);
        categsep = run.patterns.categsep(TRvector - 10); %minus 10 because we take out those 10
        sepbytrial = reshape(categsep,nTRs,10); %right now TR x Trial
        sepbytrial = sepbytrial'; %results by trial number x TR number so 10 x 15
        [~,indSort] = sort(d.stim.id);
        sepinorder = sepbytrial(indSort,:); % so this is nTrials x 15 TR's
        sepbystimD(:,(iblock-1)*nTRs + 1: iblock*nTRs,s ) = sepinorder;
        FBsepbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = sepinorder(:,5:end);
        
        % now get all the speed changes
        allSpeed = d.stim.motionSpeed;
        speedbytrial = allSpeed'; %now it's trial x TR's
        speedinorder = speedbytrial(indSort,:);
        speedbystimD(:,(iblock-1)*nTRs + 1: iblock*nTRs,s ) = speedinorder;
        
       %subtract to get the difference in speed
        FBds = diff(speedinorder,[],2);
        FBdsbystim(:,(iblock-1)*FBTRs + 1: iblock*FBTRs ) = FBds(:,4:end);
        
        timecourse_high = nan(1,15);
        timecourse_low = nan(1,15);
        %TAKE OUT FOR SPEED
        % now see how timecourses change
        for i = 1:nstim % because we're using all the stimuli now
            tcourse = sepbytrial(i,:);
            x1 = 1:length(tcourse);
            x2 = 1:.5:length(tcourse);
            y2 = interp1q(x1',tcourse',x2');
            for j = 3:length(y2) - avgRange
                if y2(j-1) < highB && y2(j) > highB %then this is a POSITIVE CROSSING POINT
                    timecourse_high(end+1,:) = y2(j-2:j+avgRange);
                elseif y2(j-1) > lowB && y2(j) < lowB %% then this is a NEGATIVE crossing point
                    timecourse_low(end+1,:) = y2(j-2:j+avgRange);
                end
            end
        end
        
        %would be +-1.25 because that's what the setting for max increase
        %or decrease?
%         for i = 1:nstim % because we're using all the stimuli now
%             tcourse = sepbytrial(i,:);
%             dscourse = FBds;
%             x1 = 1:length(dscourse);
%             x2 = 1:.5:length(dscourse);
%             y2 = interp1q(x1',dscourse',x2');
%             
%             x1s = 1:length(tcourse);
%             x2s = 1:.5:length(tcourse);
%             y2s = interp1q(x1',tcourse',x2');
%             for j = 3:length(y2) - avgRange -1
%                 if y2(j-1) < speedH && y2(j) > speedH %then this is a POSITIVE CROSSING POINT
%                     timecourse_high(end+1,:) = y2s(j-2+1:j+avgRange+1);
%                 elseif y2(j-1) > speedL && y2(j) < speedL %% then this is a NEGATIVE crossing point
%                     timecourse_low(end+1,:) = y2s(j-2+1:j+avgRange+1); %you go + 1 because of the index difference
%                 end
%             end
%         end
%         
        % look how separation evidence is changing by TR
        allDiff = diff(diff(sepinorder(:,5:end),1,2),1,2)/4;
        secondDiff = reshape(allDiff,1,numel(allDiff));
        zTR = length(secondDiff);
        allSecondDiffD(s,(iblock-1)*zTR + 1: iblock*zTR) = secondDiff;
    end
    
    
    
    %take the average of nTRs in range
    %remove feedback and just look at all datapoints
    z1 =find(FBsepbystim>=goodRange(1));
    z2 = find(FBsepbystim<=goodRange(2));
    nGoodRangeD(s) = length(intersect(z1,z2))/numel(FBsepbystim);
    nConsecD(s) = sum(diff(intersect(z1,z2))==1)/numel(FBsepbystim);
    nLowD(s) = length(find(FBsepbystim<=goodRange(1)))/numel(FBsepbystim);
    nHighD(s) = length(find(FBsepbystim>=goodRange(2)))/numel(FBsepbystim);
    vectorSepD(s,:) = reshape(FBsepbystim,1,numel(FBsepbystim));
    if isnan(timecourse_high)
        avg_highD(s,:) = nan(1,15);
    else
        avg_highD(s,:) = nanmean(timecourse_high);
    end
    if isnan(timecourse_low)
        avg_lowD(s,:) = nan(1,15);
    else
        avg_lowD(s,:) = nanmean(timecourse_low);
    end
end
%% save to plot in python
folder= '/jukebox/norman/amennen/PythonMot5';
save('compareExp5.mat','nGoodRangeD','avg_highD', 'avg_lowD',  'nLowD','nHighD','vectorSepD', 'allSecondDiffD', 'nConsecD','sepbystimD' , 'speedbystimD');
unix(['scp ' 'compareExp5.mat' ' amennen@apps.pni.princeton.edu:' folder '/' ])
