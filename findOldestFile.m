% find the most recent file in the specified location

function [filetoload] = findOldestFile(filepath,filetodir)
temp = dir(filetodir);
if ~isempty(temp)
    dates = [temp.datenum];
    names = {temp.name};
    [~,newest] = min(dates);
    filetoload = fullfile(filepath, names{newest});
else
    filetoload = [];
end