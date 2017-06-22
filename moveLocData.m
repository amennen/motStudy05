function moveLocData(subjectNum)
folder= '/jukebox/norman/amennen/PythonMot5/Loc/';
projectName = 'motStudy05';
behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
loc_dir = ['/Data1/code/' projectName '/' 'data' '/' num2str(subjectNum) '/Localizer/'];
fname = findNewestFile(loc_dir, fullfile(loc_dir, ['locpatterns' '*.mat']));
newname = ['locpat' num2str(subjectNum) '.mat'];
unix(['scp ' fname ' amennen@apps.pni.princeton.edu:' folder '/' newname])
end