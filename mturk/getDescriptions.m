nstim = 20;
Q1 = 36;
cols = 1:nstim;
cols = Q1 + 9*cols; % these are all the rows in the data table where the responses are
nsub = 54;
addpath(genpath('/Users/amennen/Documents/Norman/WordVec'))
HashVocab;

% delete rows where there aren't subjects and just hard code that part
t = readtable('NEWDESC_sp.csv'); %spellchecked first!
allquestions = table2cell(t(3:end,3:end));
firstrow = 5;
for s =1:nstim
    description{s} = allquestions{firstrow,s};
end
stopWords = {'a', 'about', 'above', 'above', 'across', 'after', 'afterwards', 'again', 'against', 'all', 'almost', 'alone', 'along', 'already', 'also','although','always','am','among', 'amongst', 'amoungst', 'amount',  'an', 'and', 'another', 'any','anyhow','anyone','anything','anyway', 'anywhere', 'are', 'around', 'as',  'at', 'back','be','became', 'because','become','becomes', 'becoming', 'been', 'before', 'beforehand', 'behind', 'being', 'below', 'beside', 'besides', 'between', 'beyond', 'bill', 'both', 'bottom','but', 'by', 'call', 'can', 'cannot', 'cant', 'co', 'con', 'could', 'couldnt', 'cry', 'de', 'describe', 'detail', 'do', 'done', 'down', 'due', 'during', 'each', 'eg', 'eight', 'either', 'eleven','else', 'elsewhere', 'empty', 'enough', 'etc', 'even', 'ever', 'every', 'everyone', 'everything', 'everywhere', 'except', 'few', 'fifteen', 'fify', 'fill', 'find', 'fire', 'first', 'five', 'for', 'former', 'formerly', 'forty', 'found', 'four', 'from', 'front', 'full', 'further', 'get', 'give', 'go', 'had', 'has', 'hasnt', 'have', 'he', 'hence', 'her', 'here', 'hereafter', 'hereby', 'herein', 'hereupon', 'hers', 'herself', 'him', 'himself', 'his', 'how', 'however', 'hundred', 'ie', 'if', 'in', 'inc', 'indeed', 'interest', 'into', 'is', 'it', 'its', 'itself', 'keep', 'last', 'latter', 'latterly', 'least', 'less', 'ltd', 'made', 'many', 'may', 'me', 'meanwhile', 'might', 'mill', 'mine', 'more', 'moreover', 'most', 'mostly', 'move', 'much', 'must', 'my', 'myself', 'name', 'namely', 'neither', 'never', 'nevertheless', 'next', 'nine', 'no', 'nobody', 'none', 'noone', 'nor', 'not', 'nothing', 'now', 'nowhere', 'of', 'off', 'often', 'on', 'once', 'one', 'only', 'onto', 'or', 'other', 'others', 'otherwise', 'our', 'ours', 'ourselves', 'out', 'over', 'own','part', 'per', 'perhaps', 'please', 'put', 'rather', 're', 'same', 'see', 'seem', 'seemed', 'seeming', 'seems', 'serious', 'several', 'she', 'should', 'show', 'side', 'since', 'sincere', 'six', 'sixty', 'so', 'some', 'somehow', 'someone', 'something', 'sometime', 'sometimes', 'somewhere', 'still', 'such', 'system', 'take', 'ten', 'than', 'that', 'the', 'their', 'them', 'themselves', 'then', 'thence', 'there', 'thereafter', 'thereby', 'therefore', 'therein', 'thereupon', 'these', 'they', 'thickv', 'thin', 'third', 'this', 'those', 'though', 'three', 'through', 'throughout', 'thru', 'thus', 'to', 'together', 'too', 'top', 'toward', 'towards', 'twelve', 'twenty', 'two', 'un', 'under', 'until', 'up', 'upon', 'us', 'very', 'via', 'was', 'we', 'well', 'were', 'what', 'whatever', 'when', 'whence', 'whenever', 'where', 'whereafter', 'whereas', 'whereby', 'wherein', 'whereupon', 'wherever', 'whether', 'which', 'while', 'whither', 'who', 'whoever', 'whole', 'whom', 'whose', 'why', 'will', 'with', 'within', 'without', 'would', 'yet', 'you', 'your', 'yours', 'yourself', 'yourselves', 'the'};
str_space='\s';
str_caps='[A-Z]';
str_ch='[a-z]';
str_nums='[0-9]';
str_hyph='-';
str_com = ',';
str_sl = '/';
% for each recall: take each word and average recall vectors
% over all recalls, average those vectors?
clear STIMULUS
stim=1;
for stim = 1:nstim
    clear stimavg
    fprintf('STIMULUS %i \n',stim)
    for s = 1:nsub
        recall = cell2mat(allquestions(s,stim));
        % take out punctuation
        Lstr1=length(recall);
        ind_space=regexp(recall,str_space);
        ind_caps=regexp(recall,str_caps);
        ind_chrs=regexp(recall,str_ch);
        ind_nums=regexp(recall,str_nums);
        ind_hyph=regexp(recall,str_hyph);
        ind_com = regexp(recall,str_com);
        ind_sl = regexp(recall,str_sl);
        mask=[ind_space ind_caps ind_chrs ind_nums ind_hyph ind_com ind_sl];
        num_str2=1:1:Lstr1;
        num_str2(mask)=[];
        str3=recall;
        str3([ind_hyph ind_com ind_sl]) = ' ';
        str3(num_str2)=[];
        
        
        % now take lowercase, split into separate words, take out stop words
        c = lower(strsplit(str3));
        out_str = c(~ismember(c,stopWords));
        nwords = length(out_str);
        if nwords > 1 && strcmp(c(end),'')
            c = c(1:end-1);
            out_str = c(~ismember(c,stopWords));
            nwords = length(out_str);
        end
        clear vec;
        if all(~strcmp(out_str, '')) % if there is a response, go through all words
            for w = 1:nwords
                word = out_str{w};
                test = VecLookup(word,vocabhash);
                found = 0;
                if ~isnan(test)
                    found = 1;
                    vec(w,:) = test;
                end
                
                while ~found
                    fprintf('\n\n\n>>>>The word %s is not in the dictionary. \n',upper(word));
                    fprintf('The sentence is: %s.\n', recall);
                    prompt = 'Do you want to enter a word if typo? type new word or y for two words\n';
                    new = input(prompt, 's');
                    if isempty(new) % if you don't want to type another word
                        vec(w,:) = nan(1,300); % make nan into vector
                        found = 1; % then get out of loop
                        %then keep vector as nan
                    elseif strcmp(new,'y')
                        w1 = input('What is the first word?\n','s');
                        w2 = input('What is the second word?\n','s');
                        % just average these two together first
                        if ~ismember(w1,stopWords)
                            w1vec = VecLookup(w1,vocabhash);
                        else
                            w1vec = nan(1,300);
                        end
                        if ~ismember(w2,stopWords)
                            w2vec = VecLookup(w2,vocabhash);
                        else
                            w2vec = nan(1,300);
                        end
                        vec(w,:) = nanmean([w1vec; w2vec]);
                        found = 1;
                    else
                        word = new;
                        vec(w,:) = VecLookup(new,vocabhash);
                        % but don't include if it's a stop word!!
                        if ~isnan(vec(w,:))
                            found = 1;
                        else
                            fprintf('The word %s is still not in the dictionary. \n',upper(new));
                        end
                        if ismember(new,stopWords)
                            vec(w,:) = nan(1,300);
                            fprintf('Stop word! Going to remove meow.\n')
                        end
                        if found
                            fprintf('We found %s thank you!\n', new);
                        end
                    end
                end
            end
            stimavg(s,:) = nanmean(vec);
        else % if there is no response
            stimavg(s,:) = nan(1,300);
        end
    end
    % now take the average for this stimulus
    STIMULUS(stim,:) = nanmean(stimavg);
end

save('stimulusVectors', 'STIMULUS')