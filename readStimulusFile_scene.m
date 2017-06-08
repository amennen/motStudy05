% read in stimulus file (string fileName) and return each line as cell in array
function stimuli = readStimulusFile_scene2(fileName,numberNeeded,picfolder)

    %% read file

    stimuli = cell(1,numberNeeded);
    %this reads in every line in the file  (i = 2 starts the first image)
    T = readtable(fileName);
    Data = table2cell(T); %Data{i,1} = filename, Data{i,3} = category (where to find it too), Data{i,6} = ID number
    input = size(Data,1);
    %% input handling
    if ~exist('numberNeeded','var') || isempty(numberNeeded) %if don't specify values needed, take all stimuli
        numberNeeded = length(input);
    end
    input = randperm(input); %permutate the file numbers that you're taking here
    
    
    %% did we grab enough materials?
 
    ALLCAT = { 'airport_terminal', 'amusement_park', 'badlands' ,'bathroom','bedroom','bridge','castle','cockpit','conference_room','dining_room',...
        'golf_course','highway','house','kitchen','lighthouse','living_room','mountain','pasture','playground','skyscraper','tower'};
    N_categories = size(ALLCAT,2);
    CATNUM = zeros(1,N_categories);
   
    if mod(N_categories,numberNeeded) == 0
        max_stim = (numberNeeded/N_categories);
        min_stim = max_stim;
    else
        max_stim = ceil(numberNeeded/N_categories);
        min_stim = max_stim - 1;
    end
    if length(input) < numberNeeded
        error('There are not enough stimuli to provide that many.');
    elseif numberNeeded < 0
        stimuli = Data(:,1);
    else
        j = 1;
        jused=[];
        for i=1:numberNeeded
            redraw =1;
            while redraw
                
                row = input(j);
                category = Data{row,3};
                IndexC = strfind(ALLCAT, category);
                catIndex = find(not(cellfun('isempty', IndexC)));
               
                if length(catIndex) >1
                    catIndex = catIndex(1);
                end
                if CATNUM(catIndex) < min_stim && exist(char(strcat(picfolder, Data(row,1)))) %check that the picture file exists
                    stimuli{i} = Data{row,1};  %take the file name
                    CATNUM(catIndex) = CATNUM(catIndex) + 1;
                    redraw = 0;
                    jused = [jused j];
                    
                elseif length(find(CATNUM>=min_stim)) == length(CATNUM)
                    if CATNUM(catIndex) < max_stim && ~ismember(j,jused) && exist(char(strcat(picfolder, Data(row,1)))) %then use this 
                        stimuli{i} = Data{row,1}; %take the file name
                        CATNUM(catIndex) = CATNUM(catIndex) + 1;
                        redraw = 0;
                        
                    else
                    end
                end
                if j == length(input) %restart
                    j=1;
                else
                    j = j+1;
                end
            end
        end
    end
    
end
