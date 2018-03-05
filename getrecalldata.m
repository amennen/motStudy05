% purpose: get real time recall classification during TR's
laptop = 0;
if laptop
    motpath = '/Users/amennen/motStudy05/';
else
    motpath = '/Users/amennen/Documents/Norman/MOT/motStudy05/';
end
behavioral_data = [motpath  'BehavioralData/'];
allSub = [1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40];
recallSession = [20 24];
nsub = length(allSub);
nstim = 10;
nTRs=4;
nruns = 2;
% easy activation, hard activation x 4 x 2 runs x n sub

%%
for s = 1:nsub
    subjectNum = allSub(s);
    easy_activation = zeros(nstim,nTRs,nruns);
    hard_activation = zeros(nstim,nTRs,nruns);
    subj_dir =['/Volumes/norman/amennen/motStudy05_transferred/datafromintelrt/data/' num2str(subjectNum)];
    behavioral_dir = [behavioral_data num2str(subjectNum)];
    for i = 1:nruns
        SESSION = recallSession(i);
        r = dir(fullfile(subj_dir, ['recallpatternsdata_' num2str(recallSession(i)) '_' '*.mat']));
        r = load(fullfile(subj_dir,r(end).name));
        reg = r.patterns.regressor.twoCond;
        reg_rel = find(reg(1,:)) + 2;
        categsep = r.patterns.categSep(reg_rel);
        sepbytrial = reshape(categsep,nTRs,20); %right now TR x Trial
        sepbytrial = sepbytrial'; %results by trial number x TR number so 20 x 4
        
        r = dir(fullfile(behavioral_dir, ['EK' num2str(recallSession(i)) '_' 'SUB'  '*.mat']));
        r = load(fullfile(behavioral_dir,r(end).name));
        trials = table2cell(r.datastruct.trials);
        stimID = cell2mat(trials(:,8));
        cond = cell2mat(trials(:,9));
        rating = cell2mat(trials(:,12));
        [~,indSort] = sort(stimID);
        sepOrdered = sepbytrial(indSort,:);
        condOrdered = cond(indSort);
        easy = find(condOrdered==2);
        hard = find(condOrdered==1);
        
        easy_activation(:,:,i) = sepOrdered(easy,:);
        hard_activation(:,:,i) = sepOrdered(hard,:);
    end
    save(fullfile(subj_dir, 'recallactivations.mat'), 'easy_activation', 'hard_activation');
end