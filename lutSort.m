function sortedAssociates = lutSort(sortedStim,unsortedStim,unsortedAssociates)
    
    % initialize
    num_to_find = length(sortedStim);
    if iscell(unsortedAssociates)
        sortedAssociates{num_to_find} = [];
    else sortedAssociates = zeros(num_to_find,1);
    end
    
    % input checking
    if length(unsortedStim) ~= length(unsortedAssociates)
        error('unsorted stimulus vectors must be of equal length');
    end
    
    % organize
    for i = 1:num_to_find
        % sorted stimuli are strings in cells
        if iscell(sortedStim)
            this_stim = sortedStim{i};
            pos = find(strcmp(unsortedStim,this_stim),1,'first');
        % sorted stimuli are values in a numeric array
        else
            this_stim = sortedStim(i);
            pos = find(unsortedStim == this_stim,1,'first');
        end
        
        % associates are strings in cells
        if iscell(unsortedAssociates)
            if isempty(pos)
                error(['sorted stimulus ' this_stim ' not found in list of unsorted stimuli!']);
            else
                sortedAssociates{i} = unsortedAssociates{pos};
            end
        % associates are values in a numeric array
        else
            if isempty(pos)
                error(['sorted stimulus ' num2str(this_stim) ' not found in list of unsorted stimuli!']);
            else
                sortedAssociates(i) = unsortedAssociates(pos);
            end
        end
    end
    

return