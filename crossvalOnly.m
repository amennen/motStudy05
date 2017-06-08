
%cross-validate only--just for checking how the classifier is doing on
%subject data after the fact! ahhh
projectName = 'motStudy03';
subvec = [3 4 5];
nsub = length(subvec);
featureSelect = 1;
allplotDir = ['/Data1/code/' projectName '/' 'Plots' '/' ];
keepTR = 15;
%training: cross-validation
SESSION = 19; %localizer task- change this if the session number changes!
ktr_vals = 15;%4:15;
nTrainTrials = 28;
useAvg = 1; %whether or not to average TRs for classifier training
penalty = 1;
saveModel = 1;
for s = 1:nsub
    for ktr = 1:length(ktr_vals)
        keepTR = ktr_vals(ktr);
        trialVec = 1:keepTR*nTrainTrials;
        trialVec = ceil(trialVec/keepTR);
    subjectNum = subvec(s);
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
    loc_dir = ['/Data1/code/' projectName '/' 'data' '/' num2str(subjectNum) '/Localizer/'];

    fname = findNewestFile(loc_dir, fullfile(loc_dir, ['locpreprocpatterns' '*.mat']));
    load(fname);
    %first cross-validate
    %print xval results
    fprintf('\n*********************************************\n');
    fprintf(sprintf('beginning model cross-validation...\n for subject%i',s));
    
    %parameters
    shiftTR = 2;
    startXVAL = tic;
    
    %first get session information
    [newpattern t] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir,keepTR); %we're only training on 4 TR's
    [tempPats tempT] = GetSessionInfoRT(subjectNum,SESSION,behavioral_dir,15);
    test_regressor = tempPats.regressor.allCond; % here it has all 15 TR's
    test_selector = tempPats.selector.allxval;
    patterns.regressor.allCond = newpattern.regressor.allCond;
    patterns.regressor.twoCond = newpattern.regressor.twoCond;
    patterns.selector.xval = newpattern.selector.xval;
    patterns.selector.allxval = newpattern.selector.allxval;
    nIter = size(patterns.selector.allxval,1);
    %shift regressor
    nCond = size(patterns.regressor.twoCond,1);
    for j = 1:nIter
        selector = patterns.selector.allxval(j,:);
        thisTestSelector = test_selector(j,:);
        %testTE = find(test_regressor(1,:));
        
        easyIdx = find(patterns.regressor.allCond(2,:));
        hardIdx = find(patterns.regressor.allCond(1,:));
        trainIdx = find(selector == 1);
        %trainIdx = intersect(hardIdx,trainIdx);
        testIdx = find(thisTestSelector == 2);
        
        % now shift indices forward
        %trainIdx = trainIdx + shiftTR;
        %testIdx = testIdx + shiftTR;
        
        trainPats = patterns.raw_sm_filt_z(trainIdx+shiftTR,:);
        testPats = patterns.raw_sm_filt_z(testIdx+shiftTR,:);
        trainTargs = patterns.regressor.twoCond(:,trainIdx);
        %testTargs = patterns.regressor.twoCond(:,testIdx);
        testTargs = tempPats.regressor.twoCond(:,testIdx);
        testAllCond = tempPats.regressor.allCond(:,testIdx);
        if featureSelect
            thr = 0.1;
            p = run_mathworks_anova(trainPats',trainTargs);
            sigVox = find(p<thr);
            trainPats = trainPats(:,sigVox);
            testPats = testPats(:,sigVox);
        end
        
        %average train pats
        trainPats_avg = zeros(nTrainTrials,size(trainPats,2));
        trainTargs_avg = zeros(2,nTrainTrials);
        for t = 1:nTrainTrials
            avgInd = find(trialVec==t);
            trainPats_avg(t,:) = mean(trainPats(avgInd,:),1);
            trainTargs_avg(:,t) = mean(trainTargs(:,avgInd),2);
        end
        if useAvg
            trainTargs = trainTargs_avg;
            trainPats = trainPats_avg;
        else
        end
        scratchpad = train_ridge(trainPats,trainTargs,penalty);
        trainedModel = scratchpad;
        if saveModel
            save(sprintf('%saverageModel',loc_dir),'trainedModel','sigVox');
        end
        [acts scratchpad] = test_ridge(testPats,testTargs,scratchpad);
        %acts is nCond x nVoxels in the mask

        %calculate AUC for JUST TARGET vs. LURE
        clear labels;
        for i = 1:length(acts)
            condition = find(testTargs(:,i));
            if condition == 1
                labels{i} = 'target';
            elseif condition == 2
                labels{i} = 'lure';
            end
        end
      [X,Y,t,AUC(j)] = perfcurve(labels,acts(1,:), 'target');
        
        %calculate AUC SEPARATELY for easy targets vs. lure && hard targets vs.
        %lure
        testTargsFour = tempPats.regressor.allCond(:,testIdx);
        testTH = find(testAllCond(1,:));
        testTE = find(testAllCond(2,:));
        testLH = find(testAllCond(3,:));
        testLE = find(testAllCond(4,:));
        
        actDiff = acts(1,:) - acts(2,:); % targ - lure activation
%         TH_timecourse(j,:) = actDiff(testTH);
%         TE_timecourse(j,:) = actDiff(testTE);
%         LH_timecourse(j,:) = actDiff(testLH);
%         LE_timecourse(j,:) = actDiff(testLE);
        
        hardIdx = find(testTargsFour(1,:)==1);
        easyIdx = find(testTargsFour(2,:)==1);
        lureIdx = find(testTargs(2,:)==1);
        
        actsHard = acts(1,[hardIdx lureIdx]);
        actsEasy = acts(1,[easyIdx lureIdx]);
        clear labelsHard;
        for i = 1:length(actsHard)
            if i <= length(hardIdx)
                labelsHard{i} = 'target';
            else
                labelsHard{i} = 'lure';
            end
        end
        [X,Y,t,AUC_hard(j)] = perfcurve(labelsHard,actsHard, 'target');
        [X,Y,t,AUC_easy(j)] = perfcurve(labelsHard,actsEasy, 'target');
        fprintf(['* Completed Iteration ' num2str(j) '; AUC = ' num2str(AUC(j)) '\n']);
        fprintf(['* Hard vs. Lure AUC = ' num2str(AUC_hard(j)) '\n']);
        fprintf(['* Easy vs. Lure AUC = ' num2str(AUC_easy(j)) '\n']);
    
    % trainPats: is in nTRs, nVoxels
        % trainTargs: is in nCat, nTRs
        trainTargs2 = trainTargs(1,:);
        trainTargs2 = trainTargs2';
        svmStruct = svmtrain(trainPats,trainTargs2);
        group = svmclassify(svmStruct,testPats);
        correct = testTargs(1,:)';
        
        accuracySVM = group == correct;
        avg_accuracySVM(s,j) = mean(accuracySVM);
        avg_hard_accuracySVM(s,j) = mean(accuracySVM(testTH));
        avg_easy_accuracySVM(s,j) = mean(accuracySVM(testTE));
    
     % look at accuracy
        [m,guess] = max(acts,[],1);
        accuracy = zeros(1,length(acts));
        %calculate AUC for JUST TARGET vs. LURE
        for i = 1:length(acts)
            condition = find(testTargs(:,i));
            if condition == 1
                labels{i} = 'target';
            elseif condition == 2
                labels{i} = 'lure';
            end
            if condition == guess(i)
                accuracy(i) = 1;
            end
        end
     avg_accuracy(s,j) = mean(accuracy);

        hard_accuracy = accuracy(testTH);
        easy_accuracy = accuracy(testTE);
        avg_hard_accuracy(s,j) = mean(hard_accuracy);
        avg_easy_accuracy(s,j) = mean(easy_accuracy);
    end
%     TH_meanTC(s,:) = mean(TH_timecourse);
%     TE_meanTC(s,:) = mean(TE_timecourse);
%     LH_meanTC(s,:) = mean(LH_timecourse);
%     LE_meanTC(s,:) = mean(LE_timecourse);
%     took out cross val with AUC just to get average timecourses
    average_AUC(s,ktr) = mean(AUC);
%     std_AUC(s,ktr) = std(AUC)/sqrt(nIter-1);
    average_hardAUC(s,ktr) = mean(AUC_hard);
    average_easyAUC(s,ktr) = mean(AUC_easy);
%     std_hardAUC = std(AUC_hard)/sqrt(nIter-1);
%     std_easyAUC = std(AUC_easy)/sqrt(nIter-1);
    xvaltime = toc(startXVAL); %end timing
    %print cross-validation results
    fprintf('\n*********************************************\n');
    fprintf('finished cross-validation...\n');
    fprintf(['* Average AUC over Iterations: ' num2str(average_AUC(s,ktr)) '\n' ]);
    fprintf(['* Average Hard vs. Lure AUC over Iterations: ' num2str(average_hardAUC(s,ktr)) '\n']);
    fprintf(['* Average Easy vs. Lure AUC over Iterations: ' num2str(average_easyAUC(s,ktr))  '\n']);
    fprintf('Cross-validation model training time: \t%.3f\n',xvaltime); 
    
    
    end
end
%%


%first average over iterations
avg_subj = mean(avg_accuracy,2);
total_avg = mean(avg_subj);
error = std(avg_subj)/sqrt(nsub -1);

hard_avg_subj = mean(avg_hard_accuracy,2);
total_hard = mean(hard_avg_subj);
error_hard = std(hard_avg_subj)/sqrt(nsub - 1);

easy_avg_subj = mean(avg_easy_accuracy,2);
total_easy = mean(easy_avg_subj);
error_easy = std(easy_avg_subj)/sqrt(nsub - 1);
fprintf(sprintf('* Average Using Ridge Regression Classifier = %.2f +/- %.2f \n',total_avg,error));
fprintf(sprintf('* Hard vs. Lure Accuracy = %.2f +/- %.2f \n',total_hard,error_hard));
fprintf(sprintf('* Easy vs. Lure Accuracy = %.2f +/- %.2f \n',total_easy,error_easy));


%% average accuracy for SVM
%first average over iterations
avg_subj = mean(avg_accuracySVM,2);
total_avg = mean(avg_subj);
error = std(avg_subj)/sqrt(nsub -1);

hard_avg_subj = mean(avg_hard_accuracySVM,2);
total_hard = mean(hard_avg_subj);
error_hard = std(hard_avg_subj)/sqrt(nsub - 1);

easy_avg_subj = mean(avg_easy_accuracySVM,2);
total_easy = mean(easy_avg_subj);
error_easy = std(easy_avg_subj)/sqrt(nsub - 1);

fprintf(sprintf('* Average Using Ridge SVM Classifier = %.2f +/- %.2f \n',total_avg,error));
fprintf(sprintf('* Hard vs. Lure Accuracy = %.2f +/- %.2f \n',total_hard,error_hard));
fprintf(sprintf('* Easy vs. Lure Accuracy = %.2f +/- %.2f \n',total_easy,error_easy));
%% now analyze over all subjects
allavg = [mean(average_AUC) ;mean(average_hardAUC); mean(average_easyAUC)];
eallavg = [std(average_AUC)/sqrt(nsub-1); std(average_hardAUC)/sqrt(nsub-1); std(average_easyAUC)/sqrt(nsub-1)];


h = figure;
barwitherr(eallavg,allavg)
set(gca,'XTickLabel' , ['All '; 'Hard'; 'Easy']);
title('Average Crossval AUC')
xlabel('Trial Type')
ylabel('AUC')
ylim([0.5 0.75])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20)
%legend('Pre MOT', 'Post MOT')
print(h, sprintf('%sxvalresults.pdf', allplotDir), '-dpdf')

%% now plot average timecourse for each category: mseb plots by subject

timeHigh = [ nanmean(TH_meanTC,1); nanmean(TE_meanTC,1) ; nanmean(LH_meanTC,1) ; nanmean(LE_meanTC,1)];
eHigh = [nanstd(TH_meanTC,[],1)/sqrt(nsub-1) ;nanstd(TE_meanTC,[],1)/sqrt(nsub-1); nanstd(LH_meanTC,[],1)/sqrt(nsub-1) ;nanstd(LE_meanTC,[],1)/sqrt(nsub-1)];
h = figure;
npts = size(TH_meanTC,2);
mseb(1:npts,timeHigh, eHigh);
title(sprintf('High Timecourse'))
xlim([3 npts-2])
set(gca, 'XTick', [3:npts-2])
set(gca,'XTickLabel',[' 1'; ' 2'; ' 3'; ' 4'; ' 5'; ' 6'; ' 7'; ' 8'; ' 9'; '10'; '11'; '12'; '13'; '14'; '15']);
ylabel('Retrieval - Control Evidence')
xlabel('Time Points')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend('Retrieve-Fast', 'Retrieve-Slow', 'Control-Fast', 'Control-Slow')
print(h, sprintf('%sevidenceFOURCAT.pdf', allplotDir), '-dpdf') 

