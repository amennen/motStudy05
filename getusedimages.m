%find any stimuli that weren't used only for training or not used at all
%(any of the other tasks or recog lures)

%if used in training it would have an id number between 21-28
%for words you have prepared cues and lurewords used for localizer and mot
%practice and rsvp!
%load stim assignment, pics and recogLures to get all images
projectName = 'motStudy02';
svec = 8:14;
NSUB = length(svec);
base_path = [fileparts(which('mot_realtime01.m')) filesep];
%load all images
ALLMATERIALS = -1;
base_path = [fileparts(which('mot_realtime01.m')) filesep];
PICLISTFILE = [base_path 'stimuli/SCREENNAMES.txt'];
candidates = readStimulusFile_evenIO(PICLISTFILE,ALLMATERIALS);
candidates_pics = transpose(candidates);

%load all words
CUELISTFILE_TARGETS = [base_path 'stimuli/text/wordpool_targets_anne.txt'];
candidates_words = readStimulusFile(CUELISTFILE_TARGETS,[]);
for s = 1:NSUB
    behavioral_dir = [base_path '/BehavioralData/' num2str(svec(s)) '/'];
    q = load([behavioral_dir 'mot_realtime01_subj_' num2str(svec(s)) '_stimAssignment.mat']);
    allPics = cat(2,q.pics,q.recogLures);
    allWords = cat(2,q.preparedCues,q.lureWords);
    fname = findNewestFile(behavioral_dir,fullfile(behavioral_dir, ['mot_realtime01_' num2str(svec(s)) '_' num2str(3)  '*.mat']));
    z = load(fname);
    training = z.stim.associate;
    trainCues = z.stim.stim;
    usedBad = setdiff(allPics,training);
    %usedWords = setdiff(allWords,trainCues);
    usedWords = allWords;
    %then open familiarization to get training pics
  
    %good pool are the training pics and unused images
    unused{s} = setdiff(candidates_pics,usedBad);
    unused_words{s} = setdiff(candidates_words,usedWords);
end
%see if any overlap between subjects
U12 = intersect(unused{1},unused{2});
U34 = intersect(unused{3},unused{4});
U56 = intersect(unused{5},unused{6});
UB1 = intersect(U12,U34);
UB2 = intersect(U56,unused{7});
anyunused = intersect(UB1,UB2);

W12 = intersect(unused_words{1},unused_words{2});
W34 = intersect(unused_words{3},unused_words{4});
W56 = intersect(unused_words{5},unused_words{6});
WB1 = intersect(W12,W34);
WB2 = intersect(W56,unused_words{7});
anyunusedwords = intersect(WB1,WB2);
