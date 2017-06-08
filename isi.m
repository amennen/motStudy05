%% inter-stimulus interval, exit on trigger
function [lag onset_time offset_time] = isi(window,duration,fontColor,trigger,device)

    if ~exist('trigger','var') || ~exist('device','var')
        trigger = 'null'; device = -1;
    end
    startTime = GetSecs();
    DrawFormattedText(window,'+','center','center',fontColor,2);
    onset_time = Screen('Flip',window);
    elapsedTime = 0;
    while (elapsedTime < duration)
        pause(0.005)
        elapsedTime = GetSecs()-startTime;
    end
    if ~strcmp(trigger,'null')
        waitForKeyboard(trigger,device)
    end
    offset_time = GetSecs();
    lag = (offset_time-startTime)-duration;
        
return