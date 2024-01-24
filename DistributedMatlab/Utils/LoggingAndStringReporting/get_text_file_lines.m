function [lines,Nlines] = get_text_file_lines(text_file,no_empty_lines,errmode)

% Open a text file and put all of the lines into a cell column vector

% Initializations

Nargin = nargin;

if Nargin < 2; no_empty_lines = false; end
if Nargin < 3; errmode = 1; end

% Check if file exists

lines = [];
Nlines = 0;

if ~exist(text_file,'file')
    if errmode == 1
        error(['Cannot find <' text_file '>']);
    elseif errmode == 2
        warning(['Cannot find <' text_file '>']);
    end
    return
end

% Open file

fid = fopen(text_file,'r');

% Allocate initial cell column vector buffer

Nlinesbuf = 100;
lines = cell(Nlinesbuf,1);

% Read and 

Nlines = 0;
while ~feof(fid)
    
    % Get next line
    
    tline = fgetl(fid);
    if ~ischar(tline)
        error(['Invalid file: ' text_file]);
    end
    
    tline = strtrim(tline);
    
    if no_empty_lines && isempty(tline)
        store_line = false;
    else
        store_line = true;
    end
    
    if store_line
        % Increment line counter
        Nlines = Nlines+1;
        % Expand line buffer if required
        if (Nlines > Nlinesbuf)
            lines = cat(1,lines,cell(Nlinesbuf,1));
            Nlinesbuf = 2*Nlinesbuf;
        end
        % Store line
        lines{Nlines} = tline;
    end
    
end

% Close file

fclose(fid);

% Trim cell column vector buffer

if (Nlines < Nlinesbuf)
    lines = lines(1:Nlines);
end

return
end
