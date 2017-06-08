clear all; close all; clc; 

dirs.working = pwd;
%% folder with paintings (with subfolders of filler, lure, target)
%dirs.pics = '//Users/yshin/Dropbox/_database/abstract_64_64_400px_memory103/';
dirs.pics = pwd;
pictypes = {'filler','ALLSIMILARIMAGES','ALLIMAGES'};
% filler: padding images
% lure: half of these guys are actually study images (1_012a), while their
%       pairs are going to be lures at test (1_012b)
% target: study items that are going to be targets at test

%how may randomized sequence do you want?
n_seq = 2;

for cseq = 1:n_seq
    my_seed = sum(100*clock);
    s = RandStream.create('mt19937ar','seed',my_seed);
    RandStream.setGlobalStream(s);

    this_dir = sprintf('memory_%.2d',cseq-1);
    mkdir(this_dir);

   cd(dirs.pics);
    % loading file names
    for i = 1:length(pictypes)   
        tmp = dir(pictypes{i});
        tmp = struct2cell(tmp);
        tmp = tmp(1,:)';
        tmp = tmp(~strncmp(tmp,'.',1));


        picnames.(pictypes{i}) = tmp;
        if i == 2
            picnames.lure = picnames.lure(2:2:end);
        end
        pictypelength.(pictypes{i}) = length(picnames.(pictypes{i}));
    end

    cd(dirs.working);

    % you can adjust them if you want
    nlearning_only = 20;
    nblock = 8;
    nprimacy = 4;
    nafter = 4;
    ntest = 8;
    nstudy_to_test = ntest*2;
    nfiller = nprimacy + nafter;

    % randomizing picture orders
    for i = 1:length(pictypes)  
        randidx.(pictypes{i}) = randperm(pictypelength.(pictypes{i}))';
    end

    randidx.filler = reshape(randidx.filler,[nfiller nblock]);
    randidx.lure = reshape(randidx.lure,[ntest nblock]);
    randidx.target = reshape(randidx.target,[ntest nblock]);

    % write txt files
    for cb = 1:nblock
        % study sequence
        fid_study(1) = fopen(fullfile(dirs.working,this_dir,sprintf('stim_study_%d_0.txt',cb-1+nlearning_only)),'w');
    %     fid_study(2) = fopen(fullfile(dirs.working,this_dir,sprintf('stim_study_%d_1.txt',cb-1+nlearning_only)),'w');
        % test sequence for left and right
        fid_test(1) = fopen(fullfile(dirs.working,this_dir,sprintf('stim_test_l_%d.txt',cb-1+nlearning_only)),'w');
        fid_test(2) = fopen(fullfile(dirs.working,this_dir,sprintf('stim_test_r_%d.txt',cb-1+nlearning_only)),'w');
        % answer key
        fid_test_answer = fopen(fullfile(dirs.working,this_dir,sprintf('answer_key_%d.txt',cb-1+nlearning_only)),'w');

        % write study sequence 
        % primacy filler
        for n = 1:nprimacy
            fprintf(fid_study(1),'%s\n',picnames.filler{randidx.filler(n,cb)}(1:end-4));
        end
        % study items that are going to be tested. 
        % to-be-targets (t) and to-be-lures (l) are alternating one another. 
        % the order of them are random draw.
        % e.g., t l t l l t t l l t l t ...
        for n = 1:ntest

            % answer key (left and right)
            oldloc = randi(2);
            newloc = 3-oldloc;

            % random draw of order
            if randi(2) == 2
                % target first
                fprintf(fid_study(1),'%s\n',picnames.target{randidx.target(n,cb)}(1:end-4));
                fprintf(fid_study(1),'%sa\n',picnames.lure{randidx.lure(n,cb)}(1:end-5));
                order_target(n,cb) = 1;
            else
                % lure first
                fprintf(fid_study(1),'%sa\n',picnames.lure{randidx.lure(n,cb)}(1:end-5));
                fprintf(fid_study(1),'%s\n',picnames.target{randidx.target(n,cb)}(1:end-4));
                order_target(n,cb) = 2;
            end

            % prepare answer keys
            fprintf(fid_test(oldloc),'%s\n',picnames.target{randidx.target(n,cb)}(1:end-4));
            fprintf(fid_test(newloc),'%sb\n',picnames.lure{randidx.lure(n,cb)}(1:end-5));
            fprintf(fid_test_answer,'%d\n',oldloc);
            target_location(n,cb) = oldloc;

        end

        % recency filler
        for n = nprimacy+1:nfiller
            fprintf(fid_study(1),'%s\n',picnames.filler{randidx.filler(n,cb)}(1:end-4));
        end
        fclose all;
    end

    save(fullfile(dirs.working,this_dir,'memory_seq.mat'),'my_seed','picnames','pictypelength','pictypes','randidx','order_target','target_location');
    clearvars -except cseq dirs pictypes
end