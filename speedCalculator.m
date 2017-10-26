% check that the data is filerred properly
%addpath('/Data1/code/motStudy05/code/')
close all;
%clear all;
yokedSubjects = [20,24,36,34,11,27,17,16,35,30,39,25,37,31,38,33];
for s = 1:length(yokedSubjects)
    subject = yokedSubjects(s); % this is the yoked subject
    %%IMPORTANT: CHECK IF NEED TO MAKE NEW PATTERNS OR NOT!!-don't need to
    %%because scanNow =0 so it doesn't overwrite classifier output!!
    fprintf('calculating for subject %i\n', subject);
    for session = 21:23
        runNum = session - 21 + 1;
        fprintf('calculating for runNumber %i or session %i\n', runNum, session);
        base_path = [fileparts(which('mot_realtime05.m')) filesep];
        subjDir = [base_path 'BehavioralData' filesep num2str(subject) filesep];
        load([subjDir 'matchSubj.mat'])
        s2Dir = [base_path 'BehavioralData' filesep num2str(s2) filesep];
        r = dir(fullfile(subjDir, ['mot_realtime05_' num2str(subject) '_' num2str(session)  '*.mat']));
        r = load(fullfile(subjDir,r(end).name));
        r_s2 = dir(fullfile(s2Dir, ['mot_realtime05_' num2str(s2) '_' num2str(session)  '*.mat']));
        r_s2 = load(fullfile(s2Dir,r_s2(end).name));
        MATLAB_SAVE_FILE = fullfile(subjDir,['mot_realtime05_SIMULATED_' num2str(subject) '_' num2str(session) '_' datestr(now,'ddmmmyy_HHMM') '.mat']);
        
        % set where we're looking
        dicom_dir = fullfile('/Data1/code/motStudy05/', 'data', num2str(subject)); %where all the dicom information is FOR THAT SUBJECT
        if subject==24
            classOutputDir = fullfile(dicom_dir,['motRun' num2str(runNum)], 'classOutput2/');
        else
            classOutputDir = fullfile(dicom_dir,['motRun' num2str(runNum)], 'classOutput/');
        end
        % load previous info if needed
        if runNum > 1
            allLast = findNewestFile(subjDir,[subjDir 'mot_realtime05_SIMULATED_' num2str(subject) '_' num2str(session-1) '*']);
            last = load(allLast);
            lastSpeed = last.stim.lastSpeed; %matrix of motRun (1-3), stimID
            lastDecoding = last.stim.lastRTDecoding;
            lastDecodingFunction = last.stim.lastRTDecodingFunction;
        end
        
        OptimalForget = 0.15;
        maxIncrement = 1.25; %will also have to check this
        Kp = 5;
        Ki = .0; %changed after RT worksop from 0.01
        Kd = .5;
        
        % use the file that was used for actual feedback
        
        
        %%
        rtData.foundFn = [];
        rtData.RTVEC = {};
        rtData.allSpeeds = {};
        rtData.allTimeChanges = {};
        allMotionTRs = convertTR(r.timing.trig.wait,r.timing.plannedOnsets.motion,r.config.TR); %row,col = mTR,trialnumber
        for n = 1:10 % trial we're on
            current_speed = r.stim.speed(n);
            if runNum == 1
                current_speed = 2;
            else
                initSpeed = lastSpeed(r.stim.id(n)); %check the indexing is right!!
                current_speed = initSpeed;
                initFeedback = lastDecoding(r.stim.id(n));
                initFunction = lastDecodingFunction(r.stim.id(n));
            end
            rtData.allSpeeds{n} = [];
            rtData.allSpeeds{n}(1) = current_speed;
            rtData.allTimeChanges{n} = [];
            rtData.RTVEC{n} = [];
            for TRcounter = 1:15 % TR that the display is showing on
                if TRcounter >=4
                    thisTR = allMotionTRs(TRcounter,n); %this is the TR we're actually on KEEP THIS WAY--starts on 4, ends on 10
                    fileTR = thisTR - 1;
                    % figure out the NUMBER file that it came from
                    filename = r_s2.rtData.newestFile{thisTR};
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
                        current_speed = min([r.stim.maxspeed current_speed]);
                        current_speed = max([r.stim.minspeed current_speed]);
                        rtData.allSpeeds{n}(end+1) = current_speed;
                        rtData.allTimeChanges{n}(end+1) = GetSecs;
                    end
                end
            end
            stim.lastSpeed(r.stim.id(n)) = current_speed; %going to save it in a matrix of run,stimID
            stim.lastRTDecoding(r.stim.id(n)) = rtData.rtDecoding(fileNumber); %file 9 that's applied now
            stim.lastRTDecodingFunction(r.stim.id(n)) = rtData.rtDecodingFunction(fileNumber);
        end
        save(MATLAB_SAVE_FILE, 'stim', 'rtData');
    end
end