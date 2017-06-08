% script made on 8/3 to classify localizer patterns to see how much they
% changed from TR to TR-- compared to feedback trials
%run plot function first!! PlotDotSpeedClassifier.m
svec = [3:5 7];
nsub = length(svec);
runvec = [1 1 2 1];

prev = 1;
projectName = 'motStudy02';
SESSION = 18;
keepTR = 4; %should change to 8 maybe???? from 4 because we increased MOT
shiftTR = 2;
ntrials = 32;

for s = 1:nsub
    subjectNum = svec(s);
    runNum = runvec(s);
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
    save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
    process_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/' 'reg' '/'];
    roi_dir = ['/Data1/code/' projectName '/data/'];
    code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
    locPatterns_dir = fullfile(save_dir, 'Localizer/');
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
    addpath(genpath(code_dir));
    
    [newpattern t] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir,keepTR);
    patterns.regressor.allCond = newpattern.regressor.allCond;
    TargHardIdx = find(patterns.regressor.allCond(1,:));
    TargEasyIdx = find(patterns.regressor.allCond(2,:));
    LureHardIdx = find(patterns.regressor.allCond(3,:));
    LureEasyIdx = find(patterns.regressor.allCond(4,:));
    %load trained model
    allfn = dir([locPatterns_dir 'loctrainedModel_' num2str(runNum) '*']); %
    %take the last model saved
    load(fullfile(locPatterns_dir, allfn(end).name));
    
    %load localizer run's standard deviation and voxelInfo
    allLast = dir([locPatterns_dir 'locpatternsdata_' '*']);
    loc = load(fullfile(locPatterns_dir, allLast(end).name));
    
    %now get all trials
    allTRs = find(any(patterns.regressor.allCond(:,:),1));
    allTargs = loc.patterns.regressor.twoCond(:,allTRs);
    allPats = loc.patterns.raw_sm_filt_z(allTRs+shiftTR,loc.patterns.sigVox);
    [allActs,~] = test_ridge(allPats,allTargs,trainedModel);
    relcond = patterns.regressor.allCond(:,allTRs);
    
    categsep = -diff(allActs);
    thIdx = find(relcond(1,:));
    teIdx = find(relcond(2,:));
    lhIdx = find(relcond(3,:));
    leIdx = find(relcond(4,:));
    
    thByTrial2 = (reshape(categsep(thIdx),keepTR,length(thIdx)/keepTR))';
    thAvg(s) = mean(mean(abs(diff(thByTrial2,1,2))));
    
    teByTrial2 = (reshape(categsep(teIdx),keepTR,length(teIdx)/keepTR))';
    teAvg(s) = mean(mean(abs(diff(teByTrial2,1,2))));
    
    lhByTrial2 = (reshape(categsep(lhIdx),keepTR,length(lhIdx)/keepTR))';
    lhAvg(s) = mean(mean(abs(diff(lhByTrial2,1,2))));
    
    leByTrial2 = (reshape(categsep(leIdx),keepTR,length(leIdx)/keepTR))';
    leAvg(s) = mean(mean(abs(diff(leByTrial2,1,2))));
    
    % do for time
    allstart = keepTR:keepTR:ntrials*keepTR-keepTR;
    allend = allstart + 1;
    alldiff = [];
    for j = 1:ntrials -1
        i1 = allstart(j);
        i2 = allend(j);
        alldiff(j) = diff([categsep(i1) categsep(i2)]);
    end
    alldiffavg(s) = mean(abs(alldiff));
end
loc = mean([thAvg' teAvg' lhAvg' leAvg'],2);
avges = [mean(loc) mean(innerdiff) mean(boundarystim) ];
errors = [std(loc) std(innerdiff) std(boundarystim) ]/sqrt(nsub-1);
thisfig = figure;
barwitherr(errors,avges)
set(gca,'XTickLabel' , ['Localizer';'Inner RT '; ' Edge RT ']);
xlabel('Type')
ylabel('Average Abs Difference')
title('Absolute Value Differences by Type')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylim([0 .4])
print(thisfig, sprintf('%sdifferentDist.pdf', plotDir), '-dpdf')