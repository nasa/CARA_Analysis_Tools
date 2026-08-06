function [VCMFile] = getVCMFile(SatID,inputPath)
% getVCMFile - Get a satellite's VCM file name, checking for 5 or 9
%              digit satellite ID number variants. Return empty set if
%              no file found.
%
% Syntax: VCMFile = getVCMFile(SatID, inputPath);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    SatID      -   Satellite ID string (5 or 9 digit)
%
%    inputPath  -   Path to directory containing VCM files
%
% =========================================================================
%
% Output:
%
%   VCMFile    -   VCM file name (empty if not found)
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initialize output

VCMFile  = [];

% Search using input SatID

template = fullfile(inputPath,[SatID '*.vcm']);
dl = dir(template);
if ~isempty(dl)
    if numel(dl) > 1
        warning('More than one candidate VCM file found; using first');
        VCMFile = [dl(1).name ' (among ' numel(dl) ' candidates ']; 
    else
        VCMFile = dl.name;
    end
    return;
end

% If default file template doesn't work, then try the 5 or 9 digit name
% variants, if possible.

NumID = str2double(SatID);
if isempty(NumID) || isnan(NumID) || (NumID < 0)
    % SatID does not define a usable number, so no variant file name can be
    % constructed
    return;
end

LenID = length(SatID);

if LenID == 5
    % Try 9 digit variant
    SatIDAlt = ['0000' SatID];
    template = fullfile(inputPath,[SatIDAlt '*.vcm']);
    dl = dir(template);
elseif (LenID == 9) && strcmpi(SatID(1:4),'0000')
    % Try 5 digit variant
    SatIDAlt = SatID(5:9);
    template = fullfile(inputPath,[SatIDAlt '*.vcm']);
    dl = dir(template);
end

if ~isempty(dl)
    if numel(dl) > 1
        warning('More than one candidate VCM file found; using first');
        VCMFile = [dl(1).name ' (among ' numel(dl) ' candidates ']; 
    else
        VCMFile = dl.name;
    end
end

return;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================