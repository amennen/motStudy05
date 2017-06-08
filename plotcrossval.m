
NITER = 8;
N_SUB = 3;
for s = 1:N_SUB
    for i = 1:NITER
        
        %now see if evidence for target is related to fast vs. slow dot motion
        % targetEvidence = results_FR.iterations(i).acts(1,:);
        % training = results_FR.iterations(i).train_idx;
        % testing = results_FR.iterations(i).test_idx;
        % acc = results_FR.iterations(i).perfmet.corrects;
        
        targetEvidence = results(s).iterations(i).acts(1,:);
        training = results(s).iterations(i).train_idx;
        testing = results(s).iterations(i).test_idx;
        acc = results(s).iterations(i).perfmet.corrects; %only for testing trials
        
        target = find(results(s).iterations(i).perfmet.desireds==1);
        lure = find(results(s).iterations(i).perfmet.desireds==2);
        
        acc_target(s,i) = mean(acc(target));
        acc_lure(s,i) = mean(acc(lure));
        
         
        
    end
end
mIter_T = mean(acc_target,2);
mIter_L = mean(acc_lure,2);

meanA(1) = mean(mIter_T);
meanA(2) = mean(mIter_L);
stdA(1) = std(mIter_T)/sqrt(N_SUB - 1);
stdA(2) = std(mIter_L)/sqrt(N_SUB - 1);
figure;
barwitherr(stdA,meanA)
set(gca,'XTickLabel' , ['Targ'; 'Lure']);
title('Average Accuracy')
xlabel('Trial Type')
ylabel('Accuracy (%)')
ylim([0 0.9])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',18)

% meanA(1) = mean(evidenceTH);
% meanA(2) = mean(evidenceTE);
% meanA(3) = mean(evidenceLH);
% meanA(4) = mean(evidenceLE);
% stdA(1) = std(evidenceTH)/sqrt(length(evidenceTH)-1);
% stdA(2) = std(evidenceTE)/sqrt(length(evidenceTE)-1);
% stdA(3) = std(evidenceLH)/sqrt(length(evidenceLH)-1);
% stdA(4) = std(evidenceLE)/sqrt(length(evidenceLE)-1);
% plot(evidenceTH)
% hold on;
% plot(evidenceTE, 'r')
% plot(evidenceLE, 'g')
% plo
if crossval
meanA(1) = mean(mean(evidenceTH));
meanA(2) = mean(mean(evidenceTE));
meanA(3) = mean(mean(evidenceLH));
meanA(4) = mean(mean(evidenceLE));
stdA(1) = std(evidenceTH)/sqrt(numel(evidenceTH)-1);
stdA(2) = std(evidenceTE)/sqrt(numel(evidenceTE)-1);
stdA(3) = std(evidenceLH)/sqrt(numel(evidenceLH)-1);
stdA(4) = std(evidenceLE)/sqrt(numel(evidenceLE)-1);


    
figure;
barwitherr(stdA,meanA)
set(gca,'XTickLabel' , ['TH'; 'TE'; 'LH'; 'LE']);
title('Average Recall Evidence')
xlabel('Trial Type')
ylabel('Average Recall Evidence')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',18)

else
meanA(1) = mean(mean(evidenceTH));
meanA(2) = mean(mean(evidenceTE));
meanA(3) = mean(mean(evidenceLH));
stdA(1) = std(evidenceTH)/sqrt(numel(evidenceTH)-1);
stdA(2) = std(evidenceTE)/sqrt(numel(evidenceTE)-1);
stdA(3) = std(evidenceLH)/sqrt(numel(evidenceLH)-1);
    
figure;
barwitherr(stdA,meanA)
set(gca,'XTickLabel' , ['TH'; 'TE'; 'OM']);
title('Average Recall Evidence')
xlabel('Trial Type')
ylabel('Average Recall Evidence')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',18)

meanA(1) = mean(mean(evidenceTH));
meanA(2) = mean(mean(evidenceTE));
meanA(3) = mean(mean(evidenceLH));
stdA(1) = std(evidenceTH)/sqrt(numel(evidenceTH)-1);
stdA(2) = std(evidenceTE)/sqrt(numel(evidenceTE)-1);
stdA(3) = std(evidenceLH)/sqrt(numel(evidenceLH)-1);
    
figure;
barwitherr(stdA,meanA)
set(gca,'XTickLabel' , ['TH'; 'TE'; 'OM']);
title('Average Recall Evidence')
xlabel('Trial Type')
ylabel('Average Recall Evidence')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',18)

end
%organize into trials and plot average accuracy for each trial
%% now look at accuracy

if crossval
meanacc(1) = mean(aTH);
meanacc(2) = mean(aTE);
meanacc(3) = mean(aLH);
meanacc(4) = mean(aLE);

stdacc(1) = std(aTH)/sqrt(length(aTH)-1);
stdacc(2) = std(aTE)/sqrt(length(aTE)-1);
stdacc(3) = std(aLH)/sqrt(length(aLH)-1);
stdacc(4) = std(aLE)/sqrt(length(aLE)-1);

figure;
barwitherr(stdacc,meanacc)
set(gca,'XTickLabel' , ['TH'; 'TE'; 'LH'; 'LE']);
title('Average Accuracy')
xlabel('Trial Type')
ylabel('Average Recall Evidence')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',18)

else
     meanacc(1) = mean(aTH);
meanacc(2) = mean(aTE);
meanacc(3) = mean(aLH);

stdacc(1) = std(aTH)/sqrt(length(aTH)-1);
stdacc(2) = std(aTE)/sqrt(length(aTE)-1);
stdacc(3) = std(aLH)/sqrt(length(aLH)-1);

figure;
barwitherr(stdacc,meanacc)
set(gca,'XTickLabel' , ['TH'; 'TE'; 'LH']);
title('Average Accuracy')
xlabel('Trial Type')
ylabel('Average Recall Evidence')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',18)
end
  
    
%% plot by trial

nTR = trialDurTR+1;
tTH = reshape(evidenceTH, nTR,length(evidenceTE)/nTR );
tTE = reshape(evidenceTE, nTR ,length(evidenceTE)/nTR );
tLH = reshape(evidenceLH, nTR ,length(evidenceLH)/nTR );
byTrialTH = tTH';
byTrialTE = tTE';
byTrialLH = tLH';
if ~isempty(LE)
    tLE = reshape(evidenceLE, nTR ,length(evidenceLH)/nTR );
    %evLE = targetEvidence(tLE');
end
% evTH = targetEvidence(tTH');
% evTE = targetEvidence(tTE');
% evLH = targetEvidence(tLH');

figure;
plot(mean(tTH'));
hold on;
plot(mean(tTE'), 'g');
plot(mean(tLH'), 'r');
if crossval
    plot(mean(tLE'), 'r')
end

% figure;
% plot( mean(byTrialTH(11:20,:)) );
% hold on;
% plot(mean(byTrialTE(11:20,:)) , 'm');
% plot(mean(byTrialLH(11:20,:)) , 'r');
% legend(['HARD'; 'EASY'; 'OMIT'])
% 
% figure;
% plot(mean(byTrialTH(11:20,:)) - mean(byTrialTH(1:10,:)) );
% hold on;
% plot(mean(byTrialTE(11:20,:)) - mean(byTrialTE(1:10,:)) , 'm');
% plot(mean(byTrialLH(11:20,:)) - mean(byTrialLH(1:10,:)) , 'r');
% if crossval
%     plot(mean(tLE'), 'k')
% end
% 
