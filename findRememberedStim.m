function [allRemID,easyRemID,hardRemID] = findRememberedStim(subjectNum)
projectName = 'motStudy03';
base_path = [fileparts(which('mot_realtime02.m')) filesep];
behavioral_dir = [base_path 'BehavioralData/' num2str(subjectNum) '/'];
SESSION = 20;
r = dir(fullfile(behavioral_dir, ['EK' num2str(SESSION) '_' 'SUB'  '*.mat']));
r = load(fullfile(behavioral_dir,r(end).name));
trials = table2cell(r.datastruct.trials);
stimID = cell2mat(trials(:,8));
cond = cell2mat(trials(:,9));
rating = cell2mat(trials(:,12));
easy = find(cond==2);
hard = find(cond==1);
rating_easy = rating(easy);
rating_hard = rating(hard);
[~, horder] = sort(stimID(find(cond==1)));
[~, eorder] = sort(stimID(find(cond==2)));
easy_ordered = rating_easy(eorder);
hard_ordered = rating_hard(horder);

easy_rem = find(easy_ordered>1);
hard_rem = find(hard_ordered>1);

IDeasy = stimID(easy);
IDhard = stimID(hard);

easyRemID = IDeasy(easy_rem);
hardRemID = IDhard(hard_rem);

allRemID = [easyRemID' hardRemID'];
end