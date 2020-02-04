function log_string(logfid,logstr,logging,displaying)

% Write a string to log file and (optionally) display

if (nargin < 1); logfid = []; end;
if (nargin < 2); logstr = ''; end;
if (nargin < 3); logging = []; end;
if (nargin < 3); displaying = []; end;

if isempty(logfid)
    logging = false;
elseif isempty(logging)
    logging = true;
end

if isempty(displaying)
    displaying = true;
end

if logging; fprintf(logfid,'%s\n',logstr); end

if displaying; disp(logstr); end

return