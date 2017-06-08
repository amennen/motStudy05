act = -1:1;
OptimalForget = 0.05;
InitialSpeed = 5;
Relation = 'linear';
MeanEv = 0;
StdEv = 0;
maxIncrement = 2;
Scale = 150;
IntegrateTime = 3; %total number of time points INCLUDING current one
ControlFunc = 'cubic';
dotSpeed = 0:30;
recallEvidence = -(.05)*dotSpeed + 0.5 + StdEv.*randn + MeanEv;

thisfig = figure;
plot(dotSpeed,recallEvidence,'color', [207 127 102]/255 ,'LineWidth', 2);
xlabel('Dot Speed')
ylabel('Retrieval - Control Evidence')
print(thisfig, sprintf('%sfb_goodspeedsdist.pdf', '/Users/amennen/Documents/Norman/generalsMISC/'), '-dpdf')
