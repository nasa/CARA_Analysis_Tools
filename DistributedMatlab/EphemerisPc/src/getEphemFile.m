function [EphemSatID,EphemFile] = getEphemFile(SatID,inputPath)
% getEphemFile - Get a satellite's ephemeris file name, checking for 5
%                or 9 digit satellite ID number variants. Return empty
%                set if no file found.
%
% Syntax: [EphemSatID, EphemFile] = getEphemFile(SatID, inputPath);
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
%    inputPath  -   Path to directory containing ephemeris files
%
% =========================================================================
%
% Output:
%
%   EphemSatID -   Satellite ID string of found ephemeris file
%                  (empty if not found)
%
%   EphemFile  -   Full path to ephemeris file
%                  (empty if not found)
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Try the default ephemeris file name

EphemSatID = SatID;
EphemFile = fullfile(inputPath,[SatID '.eci']);
EphemFileExists = exist(EphemFile,'file');

if EphemFileExists
    return;
end

% If default file doesn't exist, then try the 5 or 9 digit name variants,
% if possible

EphemSatID = [];
EphemFile  = [];

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
    EphemFileAlt = fullfile(inputPath,[SatIDAlt '.eci']);
    if exist(EphemFileAlt,'file')
        EphemSatID = SatIDAlt;
        EphemFile  = EphemFileAlt;
    end
elseif (LenID == 9) && strcmpi(SatID(1:4),'0000')
    % Try 5 digit variant
    SatIDAlt = SatID(5:9);
    EphemFileAlt = fullfile(inputPath,[SatIDAlt '.eci']);
    if exist(EphemFileAlt,'file')
        EphemSatID = SatIDAlt;
        EphemFile  = EphemFileAlt;
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