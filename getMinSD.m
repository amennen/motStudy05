nrun = 3;
for i = 1:nrun
  
    folder = ['motRun' num2str(i)];
    cd(folder);
    fn = dir(['motpatternsdata_' '*']);
    n = load(fullfile(fn(end).name));
    sd = n.patterns.runStd;
    min(sd)
    cd ../
end