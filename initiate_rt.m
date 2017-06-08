function [mainWindow windowSize colors device trigger workingDir logName matlabSaveFile ivx fmri slack] = initiate_rt(subject,session,sessionStrings,short_string,long_string,speed,colors)

    %% declarations
    AssertOpenGL
    DEFAULT_MONITOR_SIZE = [1024 768];
    DEBUG_MONITOR_SIZE = [700 500];
    MAINFONT = 'Arial';
    MAINFONTSIZE = 30;
    KEYBOARD_TRIGGER = 'Return';
    SCAN_TRIGGER = '='; % skyra PST
    ivx = [];
    subjectDir = [];
    fmri = 0;
    %% colors
    colors.WHITE = [255 255 255];
    colors.BLACK = [0 0 0];
    colors.GREY = (colors.BLACK + colors.WHITE)./2;
    colors.GREEN = [0 255 0];
    colors.RED = [255 0 0];
    colors.BLUE = [0 0 255];

     % default keypress handling
    device = -1;
    trigger = KEYBOARD_TRIGGER; % keyboard
    %% helpful input reminder for incorrect syntax
    if subject <= 0 || session <= 0 || session > length(sessionStrings)
        disp(' '); disp('********SESSION MAP********')
        for i=1:length(sessionStrings)
            disp([num2str(i) ': ' sessionStrings{i}]);
        end
        disp('***************************'); disp(' '); 
        disp('use valid session and syntax: my_experiment(subject_number,session,[debug_SPEED/optional])')
        mainWindow=[]; windowSize=[]; device=[]; trigger=[]; workingDir=[]; subjectDir=[]; logName=[]; matlabSaveFile=[];
        return
    end

    %% disable input / catch errors
    dbstop if error
    ListenChar(2); % disables keyboard input to Matlab command window
    

    %% find working dir (figure out which computer we're using)
    try %work computer
        ls('/Volumes/Macintosh HD/Users/amennen/Documents/Norman/MOT/motStudy03/');
        workingDir = '/Volumes/Macintosh HD/Users/amennen/Documents/Norman/MOT/motStudy04/';
        windowSize.degrees = [35 30];
        [keyboardIndices, productNames, allInfos] = GetKeyboardIndices;
        device = keyboardIndices(find(strcmp(productNames, '')));
        addpath(genpath('/Users/amennen/mychanges_sptb/'));
    catch
        try % my laptop
            ls('/Users/amennen/motStudy04/')
            workingDir = '/Users/amennen/motStudy04/';
            windowSize.degrees = [35 30];
        catch
            try  %computer testing room
                ls('/Users/normanlab/mot_study/')
                workingDir = '/Users/normanlab/motStudy04/';
                windowSize.degrees = [35 30];
                addpath(genpath('/Users/normanlab/mychanges_sptb/'));
                catch
                    try %Skyra
                        ls('/Data1/code/motStudy03/')
                        workingDir = '/Data1/code/motStudy04/code/';
                        fmri  = 1;
                        % special scanner keypress input
%                         if debug_mode
%                             device = -1;
%                         else
%                             device = PniKbGetDeviceNamed('Xkeys');
%                             %                         device = PniKbStartEverything;
%                         end
                        addpath(genpath('/Data1/code/SPTBanne'))
                        trigger = SCAN_TRIGGER;
                        [keyboardIndices, productNames, allInfos] = GetKeyboardIndices;
                        z = strfind(productNames, 'Xkeys');
                        %z = strfind(productNames, 'Dell Dell');
                        deviceIND = find(~cellfun(@isempty,z));
                        device = keyboardIndices(deviceIND);

                        windowSize.degrees = [51 30];
                    catch
                        error('Can''t find working directory');
                    end
               % end
            end
        end
    end

        %% initiate graphics
    sca            
    pause on
    Screen('Preference', 'SkipSyncTests', 2); % 0 for screentest
    if speed > 0
        debug_mode = 1;
        screenNumber = 0;
        windowSize.pixels = DEBUG_MONITOR_SIZE;
    else
%         HideCursor;
        debug_mode = 0;
        Screen('Preference', 'SkipSyncTests', 2);
        screens = Screen('Screens');
        screenNumber = max(screens);
        resolution = Screen('Resolution',screenNumber);
        if fmri
            windowSize.pixels = [resolution.width/2 resolution.height];
        else
            windowSize.pixels = [resolution.width resolution.height];
        end
    end
    [mainWindow, null] = Screen('OpenWindow',screenNumber,colors.BGCOLOR,[0 0 windowSize.pixels]);
    ifi = Screen('GetFlipInterval', mainWindow);
    slack  = ifi/2;
    if windowSize.pixels(2) > windowSize.pixels(1)
        SIZE_AXIS = 1;
    else SIZE_AXIS = 2;
    end
    MAINFONTSIZE = round(30 * (windowSize.pixels(SIZE_AXIS) / DEFAULT_MONITOR_SIZE(SIZE_AXIS))); % scales font by window size

    Screen('TextFont',mainWindow,MAINFONT);
    Screen('TextSize',mainWindow,MAINFONTSIZE);
    Screen('Flip',mainWindow);
    
   
    windowSize.degrees_per_pixel = windowSize.degrees ./ windowSize.pixels;


    %% initiate file handling

    % log file
    logName = ['subj' int2str(subject) '.txt'];


    % matlab save file
    matlabSaveFile = ['mot_realtime04_' num2str(subject) '_' num2str(session) '_' datestr(now,'ddmmmyy_HHMM') '.mat'];

return