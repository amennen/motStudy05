function [p] = run_mathworks_anova(pat,regs)

% -- building of vector conds and vector groups

nVox    = size(pat,1);
nConds  = size(regs,1); 
groups  = [];
dataIdx = [];

for c=1:nConds
  theseIdx = find(regs(c,:)==1);
  dataIdx  =[dataIdx,theseIdx];
  groups   =[groups,repmat(c,1,length(theseIdx))];
end

% run the anova and save the p's
p = zeros(nVox,1);

for j=1:nVox
  if mod(j,10000) == 0
    % disp( sprintf('anova on %i of %i',j,nVox) );
    fprintf('.');
  end
  p(j) = anova1(pat(j,dataIdx'),groups,'off');
end   