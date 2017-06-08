% purpose: get things ready for Python analysis:
%check if current data has good classification using first 4 vs. all TR's
%need: locpatterns data for each subject
%subject numbers
%what to name them
%final path
folder= '/jukebox/norman/amennen/PythonMot3';
subjectVec = [3 4 5 6 7 8 9 11 12 13];
%subjectVec = 9; %just put in the newest subject if you've already done for a lot
projectName = 'motStudy03';
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
    
    unix(['scp ' fname ' amennen@apps.pni.princeton.edu:' folder '/' newname])
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