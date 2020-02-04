function [strparts, Nstrparts] = string_parts(string, delimiter)
%
% Calculate the parts of a string separated by delimiters.
%

% Default delimiter is space.

if (nargin == 1)
    delimiter = ' ';
end

% Eliminate adjacent occurances of delimiter (turns out not to be
% required).

% still_eliminating = true;
% delimiter_twice   = [delimiter delimiter];
% oldstring         = string;
% 
% while (still_eliminating)    
%     newstring = strrep(oldstring, delimiter_twice, delimiter);
%     if (length(newstring) == length(oldstring))
%         still_eliminating = false;
%     else
%         oldstring         = newstring;
%     end
% end
%
% string = newstring;

% Separate the delimited parts.

rem  = string;
Ntok = 0;
Nrem = 1;

while (Nrem > 0)
   [tok, rem] = strtok(rem, delimiter); 
   if (Ntok == 0) 
       strparts = cellstr(tok);
   else
       strparts = [strparts cellstr(tok)];
   end
   Ntok = Ntok+1;
   Nrem = numel(rem);
end

Nstrparts = Ntok;

% Check for null string.

if (Nstrparts == 1)
    if (length(char(strparts(1))) == 0) Nstrparts=0; end
end