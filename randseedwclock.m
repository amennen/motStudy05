function seed = randseedwclock()
% function seed = randseedwclock()

seed = sum(100*clock);
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));

end