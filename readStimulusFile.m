% read in stimulus file (string fileName) and return each line as cell in array
function stimuli = readStimulusFile(fileName,numberNeeded,serialRead)

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
    
    %% input handling
    if ~exist('numberNeeded','var') || isempty(numberNeeded)
        numberNeeded = length(input);
    end
    if ~exist('serialRead','var') || isempty(serialRead) || ~serialRead
        input = input(randperm(length(input)));
    end
    
    %% did we grab enough materials?
    if length(input) < numberNeeded
        error('There are not enough stimuli to provide that many.');
    elseif numberNeeded < 0
        stimuli = input;
    else
        for i=1:numberNeeded
            stimuli{i} = input{i};
        end
    end
    fclose(fileHandle);
    
return