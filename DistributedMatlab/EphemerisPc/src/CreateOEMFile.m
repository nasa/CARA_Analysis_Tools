function Success = CreateOEMFile(obj,eph,EphBeg,OEMPath,ForceCreation,verbose)
% CreateOEMFile - Create a CCSDS OEM file from a BFMC ephemeris structure.
%                 (For CARA analysis team internal use)
%
% Syntax: Success = CreateOEMFile(obj, eph, EphBeg);
%         Success = CreateOEMFile(obj, eph, EphBeg, OEMPath);
%         Success = CreateOEMFile(obj, eph, EphBeg, OEMPath, ForceCreation);
%         Success = CreateOEMFile(obj, eph, EphBeg, OEMPath, ForceCreation, verbose);
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
%    obj            -   Object number (1 = Primary, 2 = Secondary)
%
%    eph            -   Ephemeris structure with fields:
%                         .T - Ephemeris times              [1xN]
%                         .X - State vectors                [6xN]
%                         .P - Covariance matrices          [6x6xN]
%                         .N - Number of points
%
%    EphBeg         -   Ephemeris begin time (days after
%                       1969-12-31 00:00:00.0)
%
%    OEMPath        -   Output directory path
%                       (optional, default = pwd)
%
%    ForceCreation  -   Flag to force file creation even if
%                       file already exists
%                       (optional, default = false)
%
%    verbose        -   Flag to display command window output
%                       (optional, default = true)
%
% =========================================================================
%
% Output:
%
%   Success        -   File creation success flag
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

%% Initializations

Nargin = nargin;
if Nargin < 6; verbose = []; end
if isempty(verbose); verbose = true; end
if Nargin < 5; ForceCreation = []; end
if isempty(ForceCreation); ForceCreation = false; end
if Nargin < 4; OEMPath = ''; end
if isempty(OEMPath); OEMPath = pwd; end

%% Define persistent variables

persistent TimeFormat DateFormat dnReference
if isempty(TimeFormat)
    % Time format
    TimeFormat = 'yyyy-mm-dd HH:MM:SS.FFF';
    % Date format
    DateFormat = 'yyyymmddHHMMSS';
    % Days after the reference time of 1969-12-31 00:00:00.0
    dnReference = datenum([1969 12 31 0 0 0]);
end

%% Define CCSDS header structure

Head.OBJECT_ID = num2str(obj);
if obj == 1
    Head.OBJECT_NAME = 'Primary';
elseif obj == 2
    Head.OBJECT_NAME = 'Secondary';
else
    Head.OBJECT_NAME = ['Obj' num2str(obj)];
end

%% Check if OEM file needs to be created

% Number of eph. times
Neph = numel(eph.T);

% Bounding dates
dnEphBeg = EphBeg+dnReference;
UT1 = datestr(dnEphBeg+eph.T(1)  ,DateFormat);
UT2 = datestr(dnEphBeg+eph.T(end),DateFormat);
% Trim seconds part of date strings if both zero
if strcmpi(UT1(13:14),'00') && strcmpi(UT2(13:14),'00')
    UT1 = UT1(1:12); UT2 = UT2(1:12);
end

OEMFile = [Head.OBJECT_NAME '_' UT1 '_to_' UT2 '_' num2str(Neph) 'pts.oem'];
OEMFull = fullfile(OEMPath,OEMFile);
if ForceCreation || ~exist(OEMFull,'file')
    if verbose > 0
        disp(['Creating CCSDS OEM file ' OEMFile]);
    end
else
    if verbose > 1
        disp(['Found previous CCSDS OEM file ' OEMFile]);
    end
    Success = true;
    return;
end

%% Create CCSDS OEM data structure

% The small angle approximation is the appropriate model when converting
% between ASW TEME and J2K
NUTModel = 'SmallAngleApprox';

% Epoch in Matlab date number format
EphBegRef = EphBeg+dnReference;
Data.Epoch = eph.T+EphBegRef;

% Allocate state and covariance arrays
Data.State = NaN(Neph,6);
Data.Cov = NaN(Neph,36); DimCovVec = [1 36];

% Convert TEME states and covariances into EME2000 (ie., J2K)
for n=1:Neph

    % OEM Epoch string in yyyy-mm-dd HH:MM:SS.FFF format.
    % This string will be written into the OEM file.
    EpochUTC = datestr(Data.Epoch(n),TimeFormat);
    
    % Calculate the Matlab date number corresonding to the EpochUTC, 
    % which can be slightly different than the original epoch (by < 1 ms)
    EpochDN = datenum(EpochUTC,TimeFormat);
    
    % Interpolate the ephemeris to the OEM Epoch
    T = EpochDN-EphBegRef; % Exact eph. time for OEM Epoch
    T = max(eph.T(1),min(eph.T(eph.N),T)); % Prevents extrapolation
    [rTEME,vTEME,cTEME] = interpStateCov(T,eph,1); % Unweighted interp.
    
    [rJ2K,vJ2K,cJ2K] = TEME2J2K_PosVelCov(rTEME,vTEME,cTEME, ...
                                          EpochUTC);
    Data.State(n,:) = [rJ2K vJ2K];
    Data.Cov(n,:) = reshape(cJ2K,DimCovVec);
    
end


% Write the OEM file
Success = CCSDSWriter(OEMFull,Data,Head);

return
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