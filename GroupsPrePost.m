% GropuPrePost: to analyze data across RT and YC groups

% analyze pre- and post- MOT recall periods

%what we want it to do:
% - open session info (pre and post)
% - open trained model
% - open patterns from pre and post scan
% - classify and subtract to see differences
% - plot
% - eventually have this as a function for every subject where date is an
% input, so is run number

% first set filepaths and information

%variables
%subjectNum = 3;
%runNum = 1;
clear all;
projectName = 'motStudy04';
onlyRem = 1; %if should only look at the stimuli that subject answered >1 for remembering in recall 1
rating_thresh = 4;
onlyForg = 0;
post = 0;
plotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' ]; %should be all
%plot dir?
svec = [4 5 6 7 8];
runvec = ones(1,length(svec));
nTRsperTrial = 8; %because 4 in task, then 2 before 2 after
if length(runvec)~=length(svec)
    error('Enter in the runs AND date numbers!!')
end
%datevec = { '1-11-17', '1-13-17'};
datevec = {  '4-20-17', '4-22-17', '4-23-17', '5-10-17', '5-11-17'};
TR = 1;
shiftTR = 4/TR;
RT = svec;
NSUB = length(svec);
for s = 1:NSUB
    subjectNum = svec(s);
    runNum = runvec(s);
    date = datevec{s};
 
    recallSession = [20 24];
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
    save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/'];
    process_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/' 'reg' '/'];
    roi_dir = ['/Data1/code/' projectName '/data/'];
    code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
    locPatterns_dir = fullfile(save_dir, 'Localizer/');
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
    addpath(genpath(code_dir));
    
    % get recall data from subject
    for i = 1:2
        

        SESSION = recallSession(i);
        rfn = dir(fullfile(save_dir, ['recallpatternsdata_' num2str(SESSION)  '*.mat']));
        rp = load(fullfile(save_dir,rfn(end).name));
        rpat = rp.patterns;
        [expat,trials,stimOrder] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir);
        testTrials = find(any(expat.regressor.allCond));
        allcond = rpat.regressor.allCond(:,testTrials);
        categSep = rpat.categSep(:,testTrials);
        %categSep = rpat.categSep(:,union(testTrials,testTrials+shiftTR+shiftTR)); % go an extra 2 for plotting purposes
        z = reshape(categSep,nTRsperTrial,20); %for 20 trials --make sure this works here!
        byTrial = z';
        RTtrials = byTrial(trials.hard,:);
        %now do in that specific order
        RTtrials = RTtrials(stimOrder.hard,:);
        OMITtrials = byTrial(trials.easy,:);
        OMITtrials = OMITtrials(stimOrder.easy,:);
        
        RTevidence(:,:,i) = RTtrials;
        OMITevidence(:,:,i) = OMITtrials;
        
    end
    
    % now find post - pre difference
    
   
        PrePostRT = RTevidence(:,:,2) - RTevidence(:,:,1);
        PrePostOMIT = OMITevidence(:,:,2) - OMITevidence(:,:,1);
        PostOnlyRT1 = RTevidence(:,:,2);
        PostOnlyOM1 = OMITevidence(:,:,2);

    PostOnlyRT(s,:) = mean(PostOnlyRT1,1);
    PostOnlyOM(s,:) = mean(PostOnlyOM1,1);
    RTavg(s,:) = mean(PrePostRT,1);
    OMITavg(s,:) = mean(PrePostOMIT,1);
    
    
end
%take data only from RT group
RT_i = find(ismember(svec,RT));
nRT = length(RT_i);
if post
    RTgroup_RT = PostOnlyRT;
    RTgroup_OM = PostOnlyOM;
else
    RTgroup_RT = RTavg(RT_i,:);
    RTgroup_OM = OMITavg(RT_i,:);
end

% calculate number of kept trials too!
h1 = figure;
%alldiffmeans = [RTavg;OMITavg];
%alldiffstd = [std(PrePostRT)/sqrt(size(PrePostRT,1)-1);std(PrePostOMIT)/sqrt(size(PrePostRT,1)-1)];
allRT = nanmean(RTgroup_RT,1);
eRT = nanstd(RTgroup_RT,[],1)/sqrt(nRT-1);
allOMIT = nanmean(RTgroup_OM,1);
eOMIT = nanstd(RTgroup_OM,[],1)/sqrt(nRT-1);
alldiffmeans = [allRT;allOMIT];
alldiffstd = [eRT;eOMIT];
mseb(1:nTRsperTrial,alldiffmeans, alldiffstd)
legend('Realtime', 'Omit')
ylim([-.5 .3])
if onlyRem
    avg_kh = mean(num_kh)/10;
    avg_ke = mean(num_ke)/10;
    text(4,-0.3, sprintf('Average kept rt = %.2f',avg_kh));
    text(4,-0.35, sprintf('Average kept omit = %.2f',avg_ke));
    title(sprintf('Post MOT Classifier; n = %i; boundary = %i',nRT,rating_thresh))
else
title(sprintf('Post MOT Classifier; n = %i',nRT))
end
set(gca, 'XTick', [1:nTRsperTrial])
%set(gca,'XTickLabel',['-2'; '-1'; ' 0'; ' 1'; ' 2'; ' 3'; ' 4'; ' 5'; '6'; '7'; '8'; '9'; ']);
ylabel('Target - Lure Evidence')
xlabel('TR (2s)')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%line([3 3], [-1 1], 'Color', 'k', 'LineWidth', 3);
%line([6 6], [-1 1], 'Color', 'k', 'LineWidth', 3);

xlim([1 nTRsperTrial])
%ylim([-.25 .25])
print(h1, sprintf('%sresults21_9sub_post_boundary5.pdf', plotDir), '-dpdf')

%%
h1 = figure;
%alldiffmeans = [RTavg;OMITavg];
%alldiffstd = [std(PrePostRT)/sqrt(size(PrePostRT,1)-1);std(PrePostOMIT)/sqrt(size(PrePostRT,1)-1)];
allRT = nanmean(YCgroup_RT,1);
eRT = nanstd(YCgroup_RT,[],1)/sqrt(nYC-1);
allOMIT = nanmean(YCgroup_OM,1);
eOMIT = nanstd(YCgroup_OM,[],1)/sqrt(nYC-1);
alldiffmeans = [allRT;allOMIT];
alldiffstd = [eRT;eOMIT];
mseb(1:nTRsperTrial,alldiffmeans, alldiffstd)
legend('Realtime', 'Omit')
title(sprintf('Post - Pre MOT Classifier Difference, YC n = %i',nYC))
set(gca, 'XTick', [1:nTRsperTrial])
%set(gca,'XTickLabel',['-2'; '-1'; ' 0'; ' 1'; ' 2'; ' 3'; ' 4'; ' 5'; '6'; '7'; '8'; '9'; ']);
ylabel('Target - Lure Evidence')
xlabel('TR (2s)')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%line([3 3], [-1 1], 'Color', 'k', 'LineWidth', 3);
%line([6 6], [-1 1], 'Color', 'k', 'LineWidth', 3);

xlim([1 nTRsperTrial])
ylim([-.25 .25])
%print(h1, sprintf('%sresults_updated0914_aonlyForg.pdf', plotDir), '-dpdf')
