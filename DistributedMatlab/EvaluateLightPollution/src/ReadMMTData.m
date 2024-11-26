function [tabl,lines,comments] = ReadMMTData(file)
% ReadMMTData - Read a data file of satellite brigtness downloaded from the 
%               MMT website.
% Syntax: [tabl,lines,comments] = ReadMMTData(file);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Read data into a table
opts = detectImportOptions(file);
opts = setvartype(opts,opts.VariableNames{1},'char');
tabl = readtable(file,opts);

% Fix the table columns
tabl.Date     = strcat(tabl.x_, repmat({' '},size(tabl.Date)), char(tabl.Date));
tabl.x_       = [];
tabl.Track    = tabl.Channel;
tabl.Channel  = tabl.Phase;
tabl.Phase    = tabl.Distance;
tabl.Distance = tabl.Penumbra;
tabl.Penumbra = tabl.Filter;
tabl.Filter   = tabl.Mag;
tabl.Mag      = tabl.StdMag;
tabl.StdMag   = tabl.Time;
tabl.Time     = [];

% Ensure that the date strings have full yyyy-mm-dd HH:MM:SS.FFFFFF format

N = numel(tabl.Date);
for n=1:N
    l = length(tabl.Date{n});
    if l < 26
        if l < 19
            error('Incomplete time');
        elseif l == 19
            tabl.Date{n} = [tabl.Date{n} '.000000'];
        elseif l < 26
            d = tabl.Date{n};
            if strcmpi(d(20),'.')
                a = l - 20;
                tabl.Date{n} = [d repmat('0',[1 a])];
            else
                error('Invalid time');
            end
        end
    elseif l > 26
        error('Invalid time');
    end
end

tabl.DateNum = datenum(tabl.Date,'yyyy-mm-dd HH:MM:SS.FFF');

if nargout > 0
    
    % Read the file into memory
    lines = fileread(file);
    lines = regexp(lines, '\r\n|\r|\n', 'split')'; 

    % Separate comments from non-comments
    ndx = startsWith(lines,'#');
    comments = lines(ndx);
    lines = lines(~ndx);
    
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================