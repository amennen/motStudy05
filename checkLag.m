% check lag in file receiving
z = load('BehavioralData/1/mot_realtime05_1_21_15Jun17_1056.mat')
newFile = z.rtData.newestFile;
delay = [];
for f =1:length(newFile)
    name = newFile{f};
    if ~isempty(name)
    startI = 5;
    fileNumber = str2double(name(startI:(length(name)-4)));
    delay(end+1) = f-fileNumber;
    
    end
end