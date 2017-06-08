% analyza data from behavioral experiment
% for each subject, calculate dot tracking accuracy, do for each condition

cd /Volumes/norman/amennen/behav_test_anne/Participant' Data'/
num_subjects = 1:21;
exclude_subj = [2 15];
subvec = setdiff(num_subjects,exclude_subj);
NUM_TASK_RUNS = 3;
all_h = [];
all_e = [];
for s = 1:length(subvec)
    cd(num2str(subvec(s)))
    for MOT_trial = 1:NUM_TASK_RUNS
        flist = dir('*DOT*mat');
        allfiles = {flist.name};
        load(allfiles{MOT_trial})
        cond = datastruct.trials.cond;
        hard_trials = find(cond==1);
        easy_trials = find(cond==2);
        acc = datastruct.trials.acc;
        average_acc = mean(acc);
        hard_acc = mean(acc(hard_trials));
        easy_acc = mean(acc(easy_trials));
        hard_all(MOT_trial) = hard_acc;
        easy_all(MOT_trial) = easy_acc;
        average_all(MOT_trial) = average_acc;
        
    end
    hard_avg = mean(hard_all);
    easy_avg = mean(easy_all);
    all_avg = mean(average_all);
    all_h = [all_h hard_avg];
    all_e = [all_e easy_avg];
    cd ..
end
hard = mean(all_h);
easy=mean(all_e);

[sorted,rank] = sort(all_h);
halfdata = floor(length(rank)/2);
lowacc = rank(1:halfdata);
highacc = rank(halfdata+1:end);