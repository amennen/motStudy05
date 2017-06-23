% look for missed triggers

subject = 3;
allreceived = [];
for s  = 21:23
    behavioral_dir = ['BehavioralData/' num2str(subject) '/'];
    z = dir(fullfile(behavioral_dir, ['mot_realtime05_*' num2str(s) '_' '*.mat']));
    z = load(fullfile(behavioral_dir,z(end).name));
    allreceived = [allreceived z.timing.trig.motion_Success];
end
Total(subject) = sum(sum(allreceived));
