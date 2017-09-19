for i = 1:3
    folder = ['motRun' num2str(i)];
    cd(folder);
    fn = dir(['motpatternsdata_' '*']);
    n = load(fullfile(fn(end).name));
    
    test_sd = std(n.patterns.raw(1:iTrial,:),[],1);
    bad = find(test_sd==0 & n.patterns.raw(iTrial,:)<20);
    cd ../
        %z(iTrial) = length(bad);
        
        %if ~isempty(bad)
        %    notMoving = [notMoving bad];
        %    patterns.constantVoxels = unique(notMoving); % MKAE THIS VALUES NONT INDECIES!!!
        %   patterns.raw_sm_filt_z(iTrial,patterns.constantVoxels) = 0;
        %end
end