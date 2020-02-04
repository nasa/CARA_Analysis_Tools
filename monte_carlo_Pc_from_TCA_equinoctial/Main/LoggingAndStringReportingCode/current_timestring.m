function [ts] = current_timestring(delimiter)

% Generate the current time string

if (~exist('delimiter','var'))
    delimiter = '';
end

if (numel(delimiter) == 0)
    d1 = '-';
    d2 = ' ';
    d3 = ':';
else
    d1 = delimiter;
    d2 = delimiter;
    d3 = delimiter;
end

cl = fix(clock);

ts = [num2str(cl(1),'%04d') d1 ...
      num2str(cl(2),'%02d') d1 ...
      num2str(cl(3),'%02d') d2 ...
      num2str(cl(4),'%02d') d3 ...
      num2str(cl(5),'%02d') d3 ...
      num2str(cl(6),'%02d')];

return