% purpose: get things ready for Python analysis:
%check if current data has good classification using first 4 vs. all TR's
%need: locpatterns data for each subject
%subject numbers
%what to name them
%final path
folder= '/jukebox/norman/amennen/PythonMot5';
subjectVec = [1 3:6];
%subjectVec = 9; %just put in the newest subject if you've already done for a lot
projectName = 'motStudy05';
for s = 1:length(subjectVec)
    subjectNum = subjectVec(s);
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
    %loc_dir = ['/Data1/code/' projectName '/' 'data' '/' num2str(subjectNum) '/Localizer/'];
    save_dir = ['/Data1/code/' projectName '/' 'data' '/' num2str(subjectNum) '/'];
    %fname = findNewestFile(loc_dir, fullfile(loc_dir, ['locpatterns' '*.mat']));
    fn1 = findNewestFile(save_dir, fullfile(save_dir, ['paraHCG_Recall1' '*.mat']));
    fn2 = findNewestFile(save_dir, fullfile(save_dir, ['paraHCG_Recall2' '*.mat']));
    newname1 = ['r1pat' num2str(s) '.mat'];
    newname2 = ['r2pat' num2str(s) '.mat'];
    
    unix(['scp ' fn1 ' amennen@apps.pni.princeton.edu:' folder '/' newname1])
    unix(['scp ' fn2 ' amennen@apps.pni.princeton.edu:' folder '/' newname2])
end

%% for anatomical recall data
% it may be easier to just load all the recall 1 and 2 files, make big
% matrices with ordered of what stimuli/category they are
folder= '/jukebox/norman/amennen/PythonMot3';

subjectVec = [3 4 5 6 7 8 9 11 12 13];
%subjectVec = 9; %just put in the newest subject if you've already done for a lot
projectName = 'motStudy03';
nTRsperTrial = 4; %4 TRs of imagining image
shiftTR = 2;
recallSession = [20 24];

for s = 1:length(subjectVec)
    allRecall = [];
    RTPAT = [];
    OMITPAT =[];
    for i = 1:2
        subjectNum = subjectVec(s);
        SESSION = recallSession(i);
        behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
        save_dir = ['/Data1/code/' projectName '/' 'data' '/' num2str(subjectNum) '/'];
        fn = findNewestFile(save_dir, fullfile(save_dir, ['paraHCG_Recall' num2str(i) '*.mat']));
        rp = load(fn);
        rpat = rp.patterns;
        [expat,trials,stimOrder] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir);
        testTrials = find(any(expat.regressor.allCond));
        allcond = rpat.regressor.allCond(:,testTrials);
        bigPat = rpat.raw_sm_filt_z(testTrials+shiftTR,:); % go an extra 2 for plotting purposes
        nVox = size(bigPat,2);
        z = reshape(bigPat,nTRsperTrial,20, nVox);
        zp = permute(z,[2 1 3]); %this is now ntrials, 4 trs, nVox so now can take average!
        avgPat = squeeze(mean(zp,2)); % now each trial has an average pattern!
        RTtrials = avgPat(trials.hard,:);
        RTtrials = RTtrials(stimOrder.hard,:);
        OMITtrials = avgPat(trials.easy,:);
        OMITtrials = OMITtrials(stimOrder.easy,:);
        
        RTPAT(:,:,i) = RTtrials;
        OMITPAT(:,:,i) = OMITtrials;
    end
    allRecall.RT = RTPAT;
    allRecall.OMIT = OMITPAT;
    newfn = ['recallPAT' num2str(s) '.mat'];
    save(newfn,'allRecall');
    unix(['scp ' newfn ' amennen@apps.pni.princeton.edu:' folder '/' newfn])
end
%%
folder= '/jukebox/norman/amennen/PythonMot5/Loc/';
folder= '/jukebox/norman/amennen/PythonMot4/Loc/';
subjectVec = [3:6]; %now made it so for all subjects and can separate them into RT/YC afterwards in python 2/17

subjectVec = 9; %just put in the newest subject if you've already done for a lot
projectName = 'motStudy05';
for s = 1:length(subjectVec)
    subjectNum = subjectVec(s);
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
    loc_dir = ['/Data1/code/' projectName '/' 'data' '/' num2str(subjectNum) '/Localizer/'];
    fname = findNewestFile(loc_dir, fullfile(loc_dir, ['retrieval_CORR_locpatterns' '*.mat']));
    fname = findNewestFile(loc_dir, fullfile(loc_dir, ['locpatterns' '*.mat']));
    fname = findNewestFile(loc_dir, fullfile(loc_dir, ['wholebrain_NOFM_locpatterns' '*.mat']));
    newname = ['locpat' num2str(subjectNum) '.mat'];
    unix(['scp ' fname ' amennen@apps.pni.princeton.edu:' folder '/' newname])
end
%%
% go to each subject
% have matrix for each sessionsubject = 20;
ntrials = 10;
nTRs = 12;
nruns = 3;
svec = [20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33];
for i = 1:length(svec)
    subject = svec(i);
    
    YOKED_EV = zeros(ntrials,nTRs,nruns);
    RT_EV = zeros(ntrials,nTRs,nruns);
    
    base_path = [fileparts(which('mot_realtime05.m')) filesep];
    subjDir = [base_path 'BehavioralData' filesep num2str(subject) filesep];
    load([subjDir 'matchSubj.mat'])
    s2Dir = [base_path 'BehavioralData' filesep num2str(s2) filesep];
    for s=1:nruns
        session = 20 + s;
        r = dir(fullfile(subjDir, ['mot_realtime05_SIMULATED_' num2str(subject) '_' num2str(session)  '*.mat']));
        r = load(fullfile(subjDir,r(end).name));
        r_s2 = dir(fullfile(s2Dir, ['mot_realtime05_' num2str(s2) '_' num2str(session)  '*.mat']));
        r_s2 = load(fullfile(s2Dir,r_s2(end).name));
        [~,indSort] = sort(r_s2.stim.id);
        for t = 1:10
            trial = indSort(t);
            RT_EV(t,:,s) = r_s2.rtData.RTVEC{trial};
            YOKED_EV(t,:,s) = r.rtData.RTVEC{trial};
        end
        
    end
    
    
    folder = '/jukebox/norman/amennen/PythonMot5/Evidence';
    fn = ['ev_' num2str(s2) '.mat'];
    savedFile = [ '../data/' fn];
    save(savedFile, 'RT_EV', 'YOKED_EV')
    unix(['scp ' savedFile ' amennen@apps.pni.princeton.edu:' folder '/' fn])
    
    
end