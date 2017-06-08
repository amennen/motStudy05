% wait_for_files
%  waits for target files in a specified directory and returns control to
%  calling function once they appear. A timer object reports which of the
%  files has been found and which the function is still waiting for.
%
% needed_files: cell array containng target file names
% search_dir: directory to search
% relative_path: set to 1 if needed_files do not include full file path
% timeout: when to give up? in seconds
%
% Written by J Poppenk --- 20-09-2011

function timeout_error = wait_for_files(needed_files,search_dir,relative_path,timeout)
    tic
    if ~iscell('needed_files')
        needed_files = {needed_files};
    end
    if ~exist('relative_path','var')
        relative_path = 1;
    end
    if ~exist('timeout','var')
        timeout = inf;
    end
    timeout_error = false;
	got_file(length(needed_files)) = 0;
    search_fn = [pwd '/file_search_status.mat'];
    save(search_fn,'needed_files','got_file');
	disp(['Waiting for files to resume processing...']);
	t=timer('period',10);
	set(t,'ExecutionMode','fixedrate','StartDelay',0.5);
	set(t,'TimerFcn',@(x,y)report_status(search_fn));
    start(t);
    if relative_path
        for i=1:length(needed_files)
            needed_files{i} = [pwd '/' needed_files{i}];
        end
    end
	while ~all(got_file) && ~timeout_error
		available = dirrec(search_dir);
		for i = 1:length(available)
			matches = strcmp(needed_files,available{i});
			got_file = got_file + matches;
        end
        if ~all(got_file)
            pause(0.05);
        end
        save(search_fn,'needed_files','got_file');
        if toc > timeout
            timeout_error = true;
        end
	end
	stop(t); delete(t);
    delete(search_fn)
	pause(0.05);
    delay = toc;
    if ~timeout_error
    	disp(['All files detected. Total wait time was ' num2str(delay) ' seconds. Processing will now continue.'])
    else disp(['Files were not detected within the ' num2str(timeout) ' second timeout window, so the search was terminated.'])
    end
return

function report_status(search_fn)
    load(search_fn)
    disp(['******File search update ******'])
    for i = 1:length(needed_files)
        if got_file(i), report = ' --> found';
        else report = ' --> searching';
        end
        disp([needed_files{i} report]);
    end
return
