% for anatomical recall data
% it may be easier to just load all the recall 1 and 2 files, make big
% matrices with ordered of what stimuli/category they are
folder= '/jukebox/norman/amennen/PythonMot5';

subjectVec = [1 3:6];
%subjectVec = 9; %just put in the newest subject if you've already done for a lot
projectName = 'motStudy05';
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
    newfn = ['recallPAT' num2str(s) '.mat'];
    save(newfn,'RTPAT', 'OMITPAT');
    unix(['scp ' newfn ' amennen@apps.pni.princeton.edu:' folder '/' newfn])
end