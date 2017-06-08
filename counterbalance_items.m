% each arg reflects # of items in a given condition
function [itemOrder itemLabels sortedStim] = counterbalance_items(conditionItems,condStrings,dont_repeat)

    % get info about conditions
    num_conds = length(conditionItems);
    if iscell(conditionItems) 
        items_per_cond = cellfun('length',conditionItems);
    else items_per_cond = conditionItems;
    end
    max_items = max(items_per_cond);
    
    % establish general order
    itemOrder = [];
    for i=1:max_items
        itemOrder = [itemOrder randperm(num_conds)];
    end
    
    % prune out excess items in each condition
    for cond = 1:num_conds
        cond_items = Shuffle(find(itemOrder == cond));
        itemOrder(cond_items(1:max_items-items_per_cond(cond))) = [];
    end
    
    % assemble condition strings
    if exist('condStrings','var') && ~isempty(condStrings)
        for i=1:length(itemOrder)
            itemLabels{i} = condStrings{itemOrder(i)};
        end
    end
    
    % randomly assign items based on condition assignment
    counter = zeros(num_conds,1);
    if iscell(conditionItems)
        % scramble sequence
        for cond = 1:length(conditionItems)
            conditionItems{cond} = Shuffle(conditionItems{cond});
        end
        % assign items
        for i = 1:length(itemOrder)
            this_cond = itemOrder(i);
            counter(this_cond) = counter(this_cond)+1;
            sortedStim{i} = conditionItems{this_cond}{counter(this_cond)};
        end
    end
    
    % swap within condition if there are any repeats
    if exist('dont_repeat','var') && dont_repeat
        doubles = 1;
        % find doubles
        while doubles
            if ischar(sortedStim{1})
                matches = strcmp(sortedStim(1:end-1),sortedStim(2:end));
            else
                matches = zeros(length(sortedStim),1);
                for i = 2:length(sortedStim)
                    if sortedStim{i} == sortedStim{i-1}
                        matches(i) = 1;
                    end
                end
            end
            
            % fix doubles
            pos = find(matches);
            if isempty(pos), doubles = 0;
            else
                for repeated = pos
                    this_item = sortedStim{repeated};
                    this_itemOrder = itemOrder(repeated);
                    alternatives = (itemLabels(repeated) == itemLabels);
                    alternatives(repeated) = 0;
                    alt_pos = find(alternatives);
                    alt_pos = alt_pos(randi(length(alt_pos)));
                    itemOrder(repeated) = itemOrder(alt_pos);
                    itemOrder(alt_pos) = this_itemOrder;
                    sortedStim{repeated} = sortedStim{alt_pos};
                    sortedStim{alt_pos} = this_item;
                end
            end
        end
    end
    
return