function [INDEXREP] = findrepeats(picsused)

base_path = [fileparts(which('behav_test_anne.m')) filesep];
PICLISTFILE = [base_path 'stimuli/FIGRIM/image_data_targets.mat'];

load(PICLISTFILE);
Data = struct2table(image_data_targets);
Data = table2cell(Data);

for i=1:length(picsused)
    IMIndex(i) = find(not(cellfun('isempty', (strfind(Data(:,1),picsused{i})))));
    IMcategory{i} = Data{IMIndex(i),3};
end
[a b c] = unique(IMcategory);
d = hist(c, length(a));
repI = find(d == 2);
repcat = a(repI);
[lia,locb] = ismember(IMcategory,repcat);
INDEXREP = find(lia==1);
end