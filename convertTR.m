%% convert time to TR number

function [TRnumber] = convertTR(onsetTime,offsetTime,TRlength)
%offsetTime = time of actual stimulus
%onsetTime is begginging of your counter
dt = offsetTime - onsetTime;
TRnumber = floor(dt/TRlength) + 1;


end