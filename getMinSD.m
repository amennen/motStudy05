nrun = 3;
for i = 1:nrun
  
    folder = ['motRun' num2str(i)];
    cd(i);
    fn = dir([lastRunHeader 'motpatternsdata_' '*']);
    n = load(fullfile(lastRunHeader,fn(end).name));
    sd = n.patterns.runStd;
    min(sd)
    cd ../
end