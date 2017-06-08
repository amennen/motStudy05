function printlog(fn,msg,varargin)

	fid = fopen(fn,'a');
    if exist('varargin','var') && ~isempty(varargin)
        msg2 = [];
        for i=1:length(varargin)
            msg2=[msg2 ', varargin{' num2str(i) '}'];
        end
        eval(['fprintf(fid,msg' msg2 ');']);
        eval(['fprintf(msg' msg2 ');']);
    else  
        fprintf(fid,msg);
        fprintf(msg);
    end
	fclose(fid);

return