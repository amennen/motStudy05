%% display some text
function startTime = displayText_specific(mainWindow,text,horiz,COLOR,WRAPCHARS,timespec)

    DrawFormattedText(mainWindow,text,'center',horiz,COLOR,WRAPCHARS);
    if ~exist('timespec', 'var')
        timespec = GetSecs;
    end
      
    startTime = Screen('Flip',mainWindow, timespec);
           
return
