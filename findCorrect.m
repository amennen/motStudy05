% find correct images

clear all;
projectName = 'motStudy02';
base_path = fileparts(which('mot_realtime01.m'));
% don't put in 22 until have subject
svec = [8 12 14 15 16 18 20:21 24 26 27 28 29 30 31 32];
RT = [8 12 14 15 15 18 21 ];
YC = [16 20 26 27 28 29];
iRT = find(ismember(svec,RT));
iYC = find(ismember(svec,YC));

NSUB = length(svec);
recallSession = [19 23];
nstim = 10;
s = 8;
i = 1;
picsPath = [base_path '/stimuli/STIM/ALLIMAGES/'];
%for s = 1:NSUB
    behavioral_dir = [base_path '/BehavioralData/' num2str(svec(s)) '/'];
    allstim = dir(fullfile(behavioral_dir, [ filesep 'mot_realtime01_subj_' num2str(svec(s)) '_' 'stim' '*.mat']));
    z = load(fullfile(behavioral_dir, allstim(end).name));
    picsOrder = z.pics;
   % for i = 1:length(recallSession)
        r = dir(fullfile(behavioral_dir, ['EK' num2str(recallSession(i)) '_' 'SUB'  '*.mat'])); 
        r = load(fullfile(behavioral_dir,r(end).name)); 
        
        trials = table2cell(r.datastruct.trials);
        stimID = cell2mat(trials(:,8));
        
        pictures = picsOrder(stimID);
        z = makeMap(pictures,1:length(pictures))
        
    %end
    
    
%end
