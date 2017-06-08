readStimulusFile([base_path 'stimuli/text/wordpool_ONLYTARGETS.txt'],stim.num_learn);
Order = randperm(59);
locI = Order(1:16);
RTI = Order(17:36);
lureI = Order(37:end);
input(locI)

for i=1:length(locI)
fprintf('%s\n',input{locI(i)})
end

for i=1:length(RTI)
fprintf('%s\n',input{RTI(i)})
end

for i=1:length(lureI)
fprintf('%s\n',input{lureI(i)})
end
