% check that the data is filerred properly
%addpath('/Data1/code/motStudy05/code/')
close all;
%clear all;
subject = 12;
run = 2;

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

fprintf('now checking subject number: %i, run number: %i\n\n\n', subject,run);
%cd(num2str(subject));
r = load(['../data/' num2str(subject) '/reg/retrieval_mask.mat']);
roiInds = find(r.mask_brain);
[gX gY gZ] = ind2sub(size(r.mask_brain),roiInds);


folder = ['../data/' num2str(subject) '/motRun' num2str(run)];
fname = findOldestFile(folder,fullfile(folder,['motpatterns'  '*.mat']));
z0 = load(fname);

folder = ['../data/' num2str(subject) '/motRun' num2str(run)];
fname = findNewestFile(folder,fullfile(folder,['motpatterns'  '*.mat']));
z = load(fname);
nTRs = size(z.patterns.raw,1);
% first check that the ones that it's setting to zero are faulty
badsigvox = [];
guess_orig = [];
guess_new = [];
if ~isempty(z.patterns.allLow)
    fprintf('found %i bad voxels: :( \n', length(z.patterns.allLow));
    for i = 1:length(z.patterns.allLow)
        fprintf('voxel # %i: (%i,%i,%i)\n', z.patterns.allLow(i),gX(z.patterns.allLow(i)),gY(z.patterns.allLow(i)),gZ(z.patterns.allLow(i)));
        issig = ismember(z.patterns.allLow(i),sigVox);
        fprintf('voxel # %i in sigVox: %i\n' , z.patterns.allLow(i),issig);
        if issig
            badsigvox = [badsigvox z.patterns.allLow(i)];
            TRrel = find(any(z.patterns.regressor.twoCond,1)) + 2;
            % now check how much weight the voxel/activity had compared to
            % others
            for t =1:length(TRrel);
                ti = TRrel(t);
                sigInd = find(sigVox==z.patterns.allLow(i));
                w1 = abs(z0.patterns.raw_sm_filt_z(ti,sigVox)) .* abs(weights(:,1)');
                % ratio1 = mean(w1/d1);
                % this will give percentage of values under current influence-maybe as long as it's < 70%?
                q1 = tiedrank(w1)/length(w1);
                quart1(t) = q1(sigInd);
                % maybe get above .75 influence of all voxels into it??
                d1 = abs(z0.patterns.raw_sm_filt_z(ti,sigVox)) * abs(weights(:,1));
                score1(t) = w1(sigInd)/mean(w1);
                w2 = abs(z0.patterns.raw_sm_filt_z(ti,sigVox)) .* abs(weights(:,2)');
                q2 = tiedrank(w1)/length(w2);
                quart2(t) = q2(sigInd);
                d2 = abs(z0.patterns.raw_sm_filt_z(ti,sigVox)) * abs(weights(:,2));
                score2(t) = w2(sigInd)/mean(w2);
                
                
                % calculate original activation
                
            end
            
            %             figure;
            %             plot(quart1);
            %             hold on;
            %             plot(quart2, 'r');
            %             ylim([0 1]);
            %             ylabel('Quartile of Influence');
            %             xlabel('TR number');
            %             title(sprintf('voxel # %i: (%i,%i,%i)\n', z.patterns.allLow(i),gX(z.patterns.allLow(i)),gY(z.patterns.allLow(i)),gZ(z.patterns.allLow(i))));
            fprintf('voxel # %i median quart %.2f\n', z.patterns.allLow(i), median(quart1));
            figure;
            plot(score1);
            hold on;
            plot(score2, 'r');
            %ylim([0 20])
            xlabel('TR number');
            ylabel('Ratio to Mean Value')
            title(sprintf('voxel # %i: (%i,%i,%i)\n', z.patterns.allLow(i),gX(z.patterns.allLow(i)),gY(z.patterns.allLow(i)),gZ(z.patterns.allLow(i))));
        end
    end
else
    fprintf('no bad voxels, yay!\n');
end

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
% check low var

% for i = 1:length(highLowVar)
%     figure;
%     plot(z.patterns.raw(:,highLowVar(i)))
%     hold on;
%     plot(z.patterns.raw_sm_filt_z(:,highLowVar(i)), 'r')
%     plot(z.patterns.raw_sm_filt(:,highLowVar(i)), 'm')
%     legend('raw', 'z', 'filtered')
%     title(['Voxel number: ' num2str(highLowVar(i)) ' or ' '(' num2str(gX(highLowVar(i))) ',' num2str(gY(highLowVar(i))) ',' num2str(gZ(highLowVar(i))) ')'])
% end

% check that it didn't miss any
allStd = std(z.patterns.raw,[],1);

% figure;
% [sorted, ind] = sort(allStd);
% plot(ind,sorted, '.')
% title('Standard Deviation by Voxel')
% xlabel('Voxel #')
%
% vox = 1;
% figure;
% plot(z.patterns.raw(:,vox))
% hold on;
% plot(z.patterns.raw_sm_filt(:,vox), 'r')
% title(['Voxel number: ' num2str(vox) ' or ' '(' num2str(gX(vox)) ',' num2str(gY(vox)) ',' num2str(gZ(vox)) ')'])
% now make sure these voxels for these trials are not in brain
%
% for i = 1:length(z.patterns.allLow)
%     figure;
%     plot(z.patterns.raw(:,z.patterns.allLow(i)))
%     hold on;
%     plot(z0.patterns.raw_sm_filt_z(:,z.patterns.allLow(i)), 'r')
%     plot(z0.patterns.raw_sm_filt(:,z.patterns.allLow(i)), 'm')
%     legend('raw', 'z', 'filtered')
%     title(['Voxel number: ' num2str(z.patterns.allLow(i)) ' or ' '(' num2str(gX(z.patterns.allLow(i))) ',' num2str(gY(z.patterns.allLow(i))) ',' num2str(gZ(z.patterns.allLow(i))) ')'])
% end



% plot brain
% figure;
% for i = 1:36
%     subplot(6,6,i)
%     imagesc(r.mask_brain(:,:,i));
% end

tr = 1:10;
vox = 3;
r2 = z.patterns.raw(tr,vox);
t2 = z.patterns.raw_sm_filt(tr,vox);
a = z.patterns.raw_sm_filt_z(tr,vox);

m = mean(z.patterns.raw_sm_filt(tr,vox),1);
s = std(z.patterns.raw_sm_filt(tr,vox),1,1);
%%
% figure;
% [n,x] = hist(changeds);
% bar(x,n/150)
% title(sprintf('Subject %i mean = %4.2f', subject, mean(changeds(changeds>0))));
% xlabel('Speed Difference')
% ylabel('Fraction of RT TRs')
% set(findall(gcf,'-property','FontSize'),'FontSize',16)