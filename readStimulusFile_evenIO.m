% read in stimulus file (string fileName) and return each line as cell in array
function stimuli = readStimulusFile_evenIO(fileName,numberNeeded,serialRead)

%% read file
stimuli = '';
fileHandle = fopen(fileName,'r+');
buffer = 0;
i=0;
while buffer ~= -1
    buffer = fgetl(fileHandle);
    i=i+1;
    if buffer ~= -1
        input{i} = buffer;
    end
end


%% did we grab enough materials?
Ntotal(1) = floor(length(input)/2);
Ntotal(2) = ceil(length(input)/2);
if numberNeeded < 1
    numberNeeded = length(input);
    numeach = Ntotal;
else
    numeach = numberNeeded/2;
    numeach = repmat(numeach,1,2);
end
takefromoutdoor = randperm(Ntotal(1),numeach(1));
takefromindoor = randperm(Ntotal(2),numeach(2)) + Ntotal(1);
if numeach == Ntotal
    allindices = Shuffle([takefromoutdoor takefromindoor]);
else
    allindices = [takefromoutdoor ; takefromindoor];
    allindices = reshape(allindices, 1,numberNeeded);
end

for i=1:numberNeeded
    stimuli{i} = input{allindices(i)};
end
fclose(fileHandle);

return

