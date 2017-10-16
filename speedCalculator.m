% check that the data is filerred properly
%addpath('/Data1/code/motStudy05/code/')
close all;
%clear all;
s1 = REALTIMESUBJECT;
subject = 12;
run = 2;

% want to get all realtime speeds
% want to save all the speeds under the trial ID and load that information


OptimalForget = 0.15;
maxIncrement = 1.25; %will also have to check this
Kp = 5;
Ki = .0; %changed after RT worksop from 0.01
Kd = .5;
folder = ['../data/' num2str(subject) '/Localizer'];
fname = findNewestFile(folder,fullfile(folder,['locpatterns'  '*.mat']));
p = load(fname);
sigVox = p.patterns.sigVox;

fname = findNewestFile(folder,fullfile(folder,['loctrainedModel'  '*.mat']));
p = load(fname);
weights = p.trainedModel.ridge.betas;

fprintf('now calculating for subject number: %i, run number: %i\n\n\n', subject,run);
%cd(num2str(subject));
r = load(['../data/' num2str(subject) '/reg/retrieval_mask.mat']);
roiInds = find(r.mask_brain);


folder = ['../data/' num2str(subject) '/motRun' num2str(run)];
fname = findOldestFile(folder,fullfile(folder,['motpatterns'  '*.mat']));
z0 = load(fname);
% use the file that was used for actual feedback

nTRs = size(z.patterns.raw,1);
% first check that the ones that it's setting to zero are faulty
badsigvox = [];
guess_orig = [];
guess_new = [];


lv = [];
for t = 1:length(z.patterns.lowVar)
    if ~isempty(z.patterns.lowVar{t})
        lv = [lv z.patterns.lowVar{t}];
    end
end
lv = unique(lv)
highLowVar = find(~ismember(lv,z.patterns.allLow))
if ~isempty(highLowVar)
    % then check if sig and add to bad voxels)
    for i=1:length(highLowVar)
        issig = ismember(highLowVar(i),sigVox);
        if issig
            badsigvox = [badsigvox highLowVar(i)];
        end
    end
end
folder = ['BehavioralData/' num2str(subject)];
if ~isempty(badsigvox)
    SESSION = 20 + run;
    thisrun = findNewestFile(folder,[folder '/mot_realtime05_' num2str(subject) '_' num2str(SESSION) '*']);
    now = load(thisrun);
    stimID = now.stim.id;
    TRrel = find(any(z.patterns.regressor.twoCond,1)) + 2;
    TRmat = reshape(TRrel,15,10);
    speed1 = {};
    speed2 = {};
    for trial = 1:10
        current_speed = 2;
        if run > 1
            %load the previous last speed for this stimulus and set
            %this as the speed
            allLast = findNewestFile(folder,[folder '/mot_realtime05_' num2str(subject) '_' num2str(SESSION-1) '*']);
            last = load(allLast);
            lastSpeed = last.stim.lastSpeed; %matrix of motRun (1-3), stimID
            initSpeed = lastSpeed(stimID(trial)); %check the indexing is right!!
            current_speed = initSpeed;
        end
        speed1{trial} = [];
        speed1{trial}(1) = current_speed;
        speed2{trial} = [];
        speed2{trial}(1) = current_speed;
        for t=1:12
            ti = TRmat(t,trial);
            % now check classifier differences
            %for t =1:length(TRrel);
            %ti = TRrel(t);
            guess_orig(t) = (weights(:,1)' *z0.patterns.raw_sm_filt_z(ti,sigVox)') -  (weights(:,2)' *z0.patterns.raw_sm_filt_z(ti,sigVox)');
            guess_new(t) = (weights(:,1)' *z.patterns.raw_sm_filt_z(ti,sigVox)') -  (weights(:,2)' *z.patterns.raw_sm_filt_z(ti,sigVox)');
            ds1(t) = PID(guess_orig(1:t),Kp,Ki,Kd,OptimalForget,maxIncrement);
            ds2(t) = PID(guess_new(1:t),Kp,Ki,Kd,OptimalForget,maxIncrement);
            current_speed1 = speed1{trial}(end) + ds1(t);
            current_speed2 = speed2{trial}(end) + ds2(t);
            current_speed1 = min([now.stim.maxspeed current_speed1]);
            current_speed1 = max([now.stim.minspeed current_speed1]);
            current_speed2 = min([now.stim.maxspeed current_speed2]);
            current_speed2 = max([now.stim.minspeed current_speed2]);
            speed1{trial}(end+1) = current_speed1;
            speed2{trial}(end+1) = current_speed2;
        end
        
    end
    changeds(:,run) = abs(ds1 - ds2);
    figure;
    plot(abs(guess_orig - guess_new), 'k', 'linewidth', 2)
    title('Abs(difference before - after correction)')
    xlabel('TR #')
    ylabel('Difference in classifier output')
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    figure;
    plot(abs(ds1 - ds2), 'k', 'linewidth', 2)
    title('Abs(change speed before - after correction)')
    xlabel('TR #')
    ylabel('Difference in \Delta speed')
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    d_speed = [];
    for t = 1:10
        d_speed = [d_speed speed1{t} - speed2{t}];
    end
    figure;
    plot(d_speed, 'k', 'linewidth', 2);
    title('Speed before - after correction')
    xlabel('TR #')
    ylabel('Difference in Speed')
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
end



%%
allMotionTRs = convertTR(runStart,timing.plannedOnsets.motion,config.TR); %row,col = mTR,trialnumber
if TRcounter >=4
    thisTR = allMotionTRs(TRcounter,n); %this is the TR we're actually on KEEP THIS WAY--starts on 4, ends on 10
    fileTR = thisTR - 1;
    % first check if there's a new file
    rtData.fileList{thisTR} = ls(classOutputDir);
    allFn = dir([classOutputDir 'vol' '*']);
    dates = [allFn.datenum];
    names = {allFn.name};
    [~,newestIndex] = max(dates);
    rtData.newestFile{thisTR} = names{newestIndex};
    % figure out the NUMBER file that it came from
    filename = rtData.newestFile{thisTR};
    startI = 5;
    fileNumber = str2double(filename(startI:(length(filename)-4)));
    % let's see if we've loaded this before
    if ~ismember(fileNumber,rtData.foundFn) && ismember(fileNumber,allMotionTRs(3:end,n)) % so only accept it if it's one of the TR's for that trial
        % now add this to the found file
        rtData.foundFn(end+1) = fileNumber;
        
        rtData.classOutputFileLoad(fileNumber) = 1;
        tempStruct = load(fullfile(classOutputDir,filename));
        rtData.rtDecoding(fileNumber) = tempStruct.classOutput;
        rtData.RTVEC{n}(end+1) = rtData.rtDecoding(fileNumber);
        rtData.rtDecodingFunction(fileNumber) = PID(rtData.RTVEC{n},Kp,Ki,Kd,OptimalForget,maxIncrement);
        
        % update speed here
        current_speed = current_speed + rtData.rtDecodingFunction(fileNumber); % apply in THIS TR what was from 2 TR's ago (indexed by what file it is) so file 3 will be applied at TR5!
        stim.changeSpeed(TRcounter,n) = rtData.rtDecodingFunction(fileNumber); %speed changed ON that TR
        % so this will update every time a new
        % file comes through!
        % you can check rtData.rtDecoding to
        % see which files are caught here
        
        
        % make sure speed is between [stim.minspeed
        % stim.maxspeed] (0,30) right now
        current_speed = min([stim.maxspeed current_speed]);
        current_speed = max([stim.minspeed current_speed]);
        rtData.allSpeeds{n}(end+1) = current_speed;
        rtData.allTimeChanges{n}(end+1) = GetSecs;
    end
end
