%cd /Volumes/norman/amennen/behav_test_anne/Participant' Data'/1
%number of participants here
cd /Volumes/norman/amennen/behav_test_anne/Participant' Data'/
num_subjects = 1:21;
exclude_subj = [2 14];
subvec = setdiff(num_subjects,exclude_subj);
N_mTurk = 9;
%subvec = [1 3 4 5 6 7 8];
subvec = setdiff(num_subjects,exclude_subj);
Ndraws = 1;
trialcolumns= [2:25]; %this maps onto trials 1-24
skip = 0;
for Ntake = 1:N_mTurk
    for s = 1:length(subvec)
        
        cd(num2str(subvec(s)))
        setup = load(['behav_subj_' num2str(subvec(s)) '_stimAssignment.mat']);
        r1F = dir('EK19*mat');
        r1 = load(r1F.name);
        trials = table2cell(r1.datastruct.trials);
        stimID = cell2mat(trials(:,9));
        cond = cell2mat(trials(:,10));
        sorted = sort(setup.pics);
        resp = xlsread(['S' num2str(subvec(s)) 'R1.xlsx']);
        resp = resp(:,trialcolumns);
        
        r2F = dir('EK23*mat');
        r2 = load(r2F.name);
        trials2 = table2cell(r2.datastruct.trials);
        stim2 = cell2mat(trials2(:,9));
        cond2 = cell2mat(trials2(:,10));
        resp2 = xlsread(['S' num2str(subvec(s)) 'R2.xlsx']);
        resp2 = resp2(:,trialcolumns);
        if skip
            skipPics = findrepeats(setup.pics);
            [c TC1] = setdiff(stimID,skipPics, 'stable');
            [c TC2] = setdiff(stim2,skipPics, 'stable');
        else
            TC1 = trialcolumns;
            TC2 = trialcolumns;
        end
        
        for d = 1:Ndraws
            
            peep = randperm(N_mTurk,Ntake);
            subResp1 = resp(peep,:);
            
            peep = randperm(N_mTurk,Ntake);
            subResp2 = resp2(peep,:);
            for i=1:length(TC1)
                if skip
                    col1 = TC1(i);
                    col2 = TC2(i);
                else
                    col1 = i;
                    col2 = i;
                end
                allresp1 = subResp1(:,col1);
                [mode1(s,i,d,Ntake) agree1(s,i,d,Ntake)] = mode(allresp1);
                rightchoiceIndex1 = stimID(col1);
                rightfile1 = setup.pics{rightchoiceIndex1};
                rightchoice1 = find(strcmp(sorted,rightfile1));
                %if length(find(allresp1==rightchoice1))==Ntake
                %if length(find(allresp1==rightchoice1))/Ntake >=.5
                if mode(allresp1) == rightchoice1
                    acc1(1,rightchoiceIndex1) = 1;
                else
                    acc1(1,rightchoiceIndex1) = 0;
                end
                acc1(2,rightchoiceIndex1) = cond(col1);
             
                allresp2 = subResp2(:,col2);
                [mode2(s,i,d,Ntake) agree2(s,i,d,Ntake)] = mode(allresp2);
                rightchoiceIndex2 = stim2(col2);
                rightfile2 = setup.pics{rightchoiceIndex2};
                rightchoice2 = find(strcmp(sorted,rightfile2));
                %if length(find(allresp2==rightchoice2))== Ntake
                %if length(find(allresp2==rightchoice2))/Ntake >=.5 
                if mode(allresp2) == rightchoice2
                    acc2(1,rightchoiceIndex2) = 1;
                else
                    acc2(1,rightchoiceIndex2) = 0;
                end
                acc2(2,rightchoiceIndex2) = cond2(col2);
            end
            
            %now index accuracy by stimulus index
            
            %omit 3 easy 2 hard 1
            %accuracydiff = acc(1,:) - acc2(1,:); %1-2
            acc1 = acc1(:,find(acc1(2,:)~=0));
            acc2 = acc2(:,find(acc2(2,:)~=0));
            omit = find(acc2(2,:)==3);
            easy = find(acc2(2,:)==2);
            hard = find(acc2(2,:)==1);
            
            hAvg = [mean(acc1(1,hard)) mean(acc2(1,hard))];
            eAvg = [mean(acc1(1,easy)) mean(acc2(1,easy))];
            oAvg = [mean(acc1(1,omit)) mean(acc2(1,omit))];
            ALLDATA(s,:,d,Ntake) = [hAvg eAvg oAvg];
            
%             EhAvg = [std(acc1(1,hard,d)) std(acc2(1,hard,d))]/sqrt(Ntake-1);
%             EeAvg = [std(acc1(1,easy,d)) std(acc2(1,easy,d))]/sqrt(Ntake-1);
%             EoAvg = [std(acc1(1,omit,d)) std(acc2(1,omit,d))]/sqrt(Ntake-1);
%             EALLDATA(s,:,d,Ntake) = [EhAvg EeAvg EoAvg];
          clear acc1
          clear acc2
        end
        cd ..
    end
end

%%IN ALLDATA(A,B,C,D)
% A = behavioral subject 1-19
% B = HARD1 HARD2 EASY1 EASY2 OMIT1 OMIT2 (1-6)
% C = draw number from 1- 100
% D = number of mTurkers that were being randomly selected (ranging from
% 1-9)

avgperm = mean(ALLDATA,3);
avgperm = reshape(avgperm, size(ALLDATA,1), size(ALLDATA,2),N_mTurk);

% now avgperm is 19 X 6 X 9 
%now take average across subjects
avgsubj = mean(avgperm,1);
avgsubj = reshape(avgsubj,size(ALLDATA,2),N_mTurk);
avgsubj = permute(avgsubj, [2 1]);
Eavgsubj = std(avgperm, [],1)/sqrt(length(subvec) - 1);
Eavgsubj = reshape(Eavgsubj,size(ALLDATA,2),N_mTurk);
Eavgsubj = permute(Eavgsubj, [2 1]);

figure;
plot(avgsubj(:,1)*100, '-r.')
hold on;
plot(avgsubj(:,2)*100, '--r.')
plot(avgsubj(:,3)*100, '-b.')
plot(avgsubj(:,4)*100, '--b.')
plot(avgsubj(:,5)*100,  '-k.') 
plot(avgsubj(:,6)*100,  '--k.')
title('Results Changing with N mTurk Raters Sampled')
xlabel('N mTurk Raters Sampled')
ylabel('Response Accuracy (%)')
ylim([75 100])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)
legend('R1 Hard', 'R2 Hard', 'R1 Easy', 'R2 Easy', 'R1 Omit', 'R2 Omit')
%%
figure;
errorbar(avgsubj(:,1), Eavgsubj(:,1),'-r.')
hold on;
errorbar(avgsubj(:,2), Eavgsubj(:,2), '--r.')
errorbar(avgsubj(:,3),  Eavgsubj(:,3),'-b.')
errorbar(avgsubj(:,4),  Eavgsubj(:,4),'--b.')
errorbar(avgsubj(:,5), Eavgsubj(:,5), '-k.') 
errorbar(avgsubj(:,6), Eavgsubj(:,6), '--k.')

%see how variability changes with draws
%agree1(s,trial,d,ntake) so 19 x 24 x 100 x 9
avgdraws1 = mean(agree1,3);
avgdraws1 = reshape(avgdraws1, length(subvec), length(trialcolumns),N_mTurk);
figure;
hist(avgdraws1(:,:,2))

%%
%now take different groups of the hard difficulty levels
HS1 = 1:9; %subjects 1-10 is really 1-9 number subjects
HS2 = 10:12;
HS4 = 13:14;
HS3 = 15:19; %should be 17-19 but minus 2 so it's 15-17

H1 = avgperm(HS1,:,:);
H2 = avgperm(HS2,:,:);
H3 = avgperm(HS3,:,:);
H4 = avgperm(HS4,:,:);

subH1 = mean(H1,1);
subH1 = reshape(subH1,size(ALLDATA,2),N_mTurk);
subH1 = permute(subH1, [2 1]);

subH2 = mean(H2,1);
subH2 = reshape(subH2,size(ALLDATA,2),N_mTurk);
subH2 = permute(subH2, [2 1]);

subH3 = mean(H3,1);
subH3 = reshape(subH3,size(ALLDATA,2),N_mTurk);
subH3 = permute(subH3, [2 1]);

subH4 = mean(H4,1);
subH4 = reshape(subH4,size(ALLDATA,2),N_mTurk);
subH4 = permute(subH4, [2 1]);

figure;
plot(subH1(:,1)*100, '-r.')
hold on;
plot(subH1(:,2)*100, '--r.')
plot(subH2(:,1)*100, '-b.')
plot(subH2(:,2)*100, '--b.')
plot(subH3(:,1)*100, '-k.') 
plot(subH3(:,4)*100, '--k.')
plot(subH4(:,1)*100, '-m.') 
plot(subH4(:,2)*100, '--m.')
title('Data broken down by MOT hard condition - HARD results only')
xlabel('N mTurk Raters Sampled')
ylabel('Response Accuracy (%)')
ylim([75 100])
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)
legend('R1 Speed1', 'R2 Speed1', 'R1 Speed2', 'R2 Speed2', 'R1 Speed3', 'R2 Speed3', 'R1 Speed4', 'R2 Speed4')


EH1 = std(ALLDATA(HS1,1:2),1)/sqrt(length(HS1) - 1);
EH2 = std(ALLDATA(HS2,1:2),1)/sqrt(length(HS2) - 1);
EH3 = std(ALLDATA(HS3,1:2),1)/sqrt(length(HS3) - 1);
EH4 = std(ALLDATA(HS4,1:2),1)/sqrt(length(HS4) - 1);
EALLHS = [EH1; EH2; EH3 ; EH4];

figure;
barwitherr(EALLHS*100,ALLHS*100)
ylim([30 105])
set(gca,'XTickLabel' , ['HS=17'; 'HS=20'; 'HS=25';'HS=32']);
title('Average Recall Matching Rate, Hard Pairs ONLY')
xlabel('Hard MOT Speed')
ylabel('Average Recall Rate (%)')
legend('Pre-mot', 'Post-mot')


%%
figure;
allsubavg = mean(ALLDATA(:,:,1,9),1);
allsubavg = reshape(allsubavg,2,3);
Eallsubavg = std(ALLDATA(:,:,1,9),[],1)/sqrt(length(subvec) - 1);
Eallsubavg = reshape(Eallsubavg,2,3);
barwitherr(Eallsubavg'*100,allsubavg'*100)
ylim([30 105])
set(gca,'XTickLabel' , ['Hard'; 'Easy'; 'Omit']);
title('Average Recall Matching Rate')
xlabel('MOT Category')
ylabel('Average Recall Rate (%)')
legend('Pre-mot', 'Post-mot')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)
%%
%plotting for each hard speed
ALLH1 = mean(ALLDATA(HS1,:,1,9),1);
ALLH1 = reshape(ALLH1,2,3);
EALLH1 = std(ALLDATA(HS1,:,1,9),[],1)/sqrt(length(HS1) - 1);
EALLH1= reshape(EALLH1,2,3);

ALLH2 = mean(ALLDATA(HS2,:,1,9),1);
ALLH2 = reshape(ALLH2,2,3);
EALLH2 = std(ALLDATA(HS2,:,1,9),[],1)/sqrt(length(HS2) - 1);
EALLH2= reshape(EALLH2,2,3);

ALLH3 = mean(ALLDATA(HS3,:,1,9),1);
ALLH3 = reshape(ALLH3,2,3);
EALLH3 = std(ALLDATA(HS3,:,1,9),[],1)/sqrt(length(HS3) - 1);
EALLH3= reshape(EALLH3,2,3);

ALLH4 = mean(ALLDATA(HS4,:,1,9),1);
ALLH4 = reshape(ALLH4,2,3);
EALLH4 = std(ALLDATA(HS4,:,1,9),[],1)/sqrt(length(HS4) - 1);
EALLH4= reshape(EALLH4,2,3);

%first plot hard cases only
ALLHARD = [ALLH1(:,1) ALLH2(:,1) ALLH3(:,1) ALLH4(:,1)];
EALLHARD = [EALLH1(:,1) EALLH2(:,1) EALLH3(:,1) EALLH4(:,1)];
figure;
barwitherr(EALLHARD'*100,ALLHARD'*100)
ylim([30 105])
set(gca,'XTickLabel' , ['SPEED 1'; 'SPEED 2'; 'SPEED 3'; 'SPEED 4']);
title('Average Recall Matching Rate: HARD TRIALS ONLY')
xlabel('HARD MOT DOT SPEED')
ylabel('Average Recall Rate (%)')
legend('Pre-mot', 'Post-mot')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)

%now easy
ALLEASY = [ALLH1(:,2) ALLH2(:,2) ALLH3(:,2) ALLH4(:,2)];
EALLEASY = [EALLH1(:,2) EALLH2(:,2) EALLH3(:,2) EALLH4(:,2)];
figure;
barwitherr(EALLEASY'*100,ALLEASY'*100)
ylim([30 105])
set(gca,'XTickLabel' , ['SPEED 1'; 'SPEED 2'; 'SPEED 3'; 'SPEED 4']);
title('Average Recall Matching Rate: EASY TRIALS ONLY')
xlabel('HARD MOT DOT SPEED')
ylabel('Average Recall Rate (%)')
legend('Pre-mot', 'Post-mot')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)

%now OMIT
ALLOMIT = [ALLH1(:,3) ALLH2(:,3) ALLH3(:,3) ALLH4(:,3)];
EALLOMIT = [EALLH1(:,3) EALLH2(:,3) EALLH3(:,3) EALLH4(:,3)];
figure;
barwitherr(EALLOMIT'*100,ALLOMIT'*100)
ylim([30 105])
set(gca,'XTickLabel' , ['SPEED 1'; 'SPEED 2'; 'SPEED 3'; 'SPEED 4']);
title('Average Recall Matching Rate: OMIT TRIALS ONLY')
xlabel('HARD MOT DOT SPEED')
ylabel('Average Recall Rate (%)')
legend('Pre-mot', 'Post-mot')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)

