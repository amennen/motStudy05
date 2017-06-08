% get dot speeds for every subject

%subvec = [12 14:15 18 20 21];
MOT_PREP = 5;
MAX_SPEED = 30;
for s = 1:length(svec)
    behavioral_dir = [fileparts(which('mot_realtime01.m')) '/BehavioralData/' num2str(svec(s)) '/'];
    fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime01_' num2str(svec(s)) '_' num2str(MOT_PREP)  '*.mat']));
    if ~isempty(fileSpeed)
        matlabOpenFile = [behavioral_dir '/' fileSpeed(end).name];
        lastRun = load(matlabOpenFile);
        hardSpeed(s) = MAX_SPEED - lastRun.stim.tGuess(end);
    end
end

% other idea: too difficult to do beforehand without setting everythign
% exactly the same so could calibrate things if they're different?