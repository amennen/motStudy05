SET_SPEED = 0;
SUBJECT = 100;
SESSION = 1;
exp_string_short = ['textlog_sess' int2str(SESSION)];
exp_string_long = ['EK' int2str(SESSION)];
COLORS.MAINFONTCOLOR = [200 200 200];
COLORS.BGCOLOR = [50 50 50];
SETUP = 1;
SESSIONSTRINGS{1} = 'GENERATE PAIRINGS'; % set up rsvp study learn associates

[mainWindow WINDOWSIZE COLORS DEVICE TRIGGER WORKING_DIR LOG_NAME MATLAB_SAVE_FILE ivx fmri SLACK] = ...
    initiate_rt(SUBJECT,SESSION,SESSIONSTRINGS,exp_string_short,exp_string_long,SET_SPEED,COLORS);
CENTER = WINDOWSIZE.pixels/2;

% prepare time-gating
% define colors for drawing
square_col = COLORS.BLACK;
stim.square_dims = round([20 20] ./ mean(WINDOWSIZE.degrees_per_pixel)); % 20ï¿½ visual angle in pixels

square_bounds = [CENTER-(stim.square_dims/2) CENTER+(stim.square_dims/2)-1];


% determine if we want a picture or black background

% draw square
offset = square_bounds(1:2);
Screen('FillRect', mainWindow, square_col, square_bounds)

% show the result

timeon = Screen('Flip', mainWindow);
