function [Success] = CCSDSWriter(EphemFile,EphemData,EphemHead)
% CCSDSWriter - Write a CCSDS Orbital Ephemeris Message (OEM) file.
%
% Syntax: [Success] = CCSDSWriter(EphemFile, EphemData);
%         [Success] = CCSDSWriter(EphemFile, EphemData, EphemHead);
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
%    EphemFile  -   File name of output CCSDS OEM format ephemeris file
%
%    EphemData  -   Structure holding ephemeris data, using the same
%                   format produced by CCSDSParser.m:
%                     .Epoch - Matlab datenum              [Nx1]
%                     .State - Position and velocity vector [Nx6]
%                     .Cov   - Covariance 6x6              [Nx36]
%
%    EphemHead  -   Structure holding header information (optional)
%                   Notes: See below for default header values.
%                          The function assumes EME2000 states and 
%                          covariances.
%
% =========================================================================
%
% Output:
%
%   Success     -   File creation success status flag
%
% =========================================================================
%
% References:
%
%   The Consultive Committee for Space Data Systems, "Orbit Data
%   Messages," CCSDS 502.0-B-3, April 2023.
%
% =========================================================================
%
% Dependencies:
%
%   set_default_param.m
%   isempty_cell.m
%
% =========================================================================
%
% Initial version: Jan 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

%% Define persistent variables
persistent pathsAdded TimeFormat
if isempty(pathsAdded)
    pathsAdded = true;
    % Path for set_default_param function
    addpath('..\Utils\General');    
    % CCSDS OEM time format
    TimeFormat = 'yyyy-mm-ddTHH:MM:SS.FFF';
    % Add the required Matlab paths
    addpath('..\Utils\General');
end

%% Initializations and defaults

% Initialize the output
Success = false;

% Initialize the optional input
if nargin < 3; EphemHead = []; end

% Initialize the header parameters
EphemHead = set_default_param(EphemHead, 'CCSDS_OEM_VERS' , '2.0');
EphemHead = set_default_param(EphemHead, 'CREATION_DATE'  , ''); % Current time later used, if empty
EphemHead = set_default_param(EphemHead, 'ORIGINATOR'     , 'CARA_ANALYSIS');
EphemHead = set_default_param(EphemHead, 'OBJECT_ID'      , '0');
EphemHead = set_default_param(EphemHead, 'OBJECT_NAME'    , ['Obj' EphemHead.OBJECT_ID]);
EphemHead = set_default_param(EphemHead, 'CENTER_NAME'    , 'Earth');
EphemHead = set_default_param(EphemHead, 'REF_FRAME'      , 'EME2000');
EphemHead = set_default_param(EphemHead, 'TIME_SYSTEM'    , 'UTC');

% Populate empty header parameters
if isempty(EphemHead.CREATION_DATE)
    EphemHead.CREATION_DATE = datestr(clock,TimeFormat);
end

% Ensure epoch and state arrays have same number of entries
Ne = numel(EphemData.Epoch);
sizeState = size(EphemData.State);
if sizeState(1) ~= Ne
    error('Dimension mismatch between Epoch and State arrays');
end

% Ensure state arrays have valid dimensions for either
%   6-D position/velocity states
% or
%   9-D position/velocity/acceleration states
Ns = sizeState(2);
if Ns ~= 6 && Ns ~= 9
    error('State vector dimension must be 6 or 9');
end

% Generate start and stop times from the eph. data
START_TIME = datestr(EphemData.Epoch(1) ,TimeFormat);
STOP_TIME  = datestr(EphemData.Epoch(Ne),TimeFormat);

%% Write the CCSDS OEM file

% Open the output file
[FID, msg] = fopen(EphemFile,'wt');
if ~isempty(msg)
    fprintf('Error opening CCSDS OEM file: %s\n     %s\n', EphemFile, msg);
    return;
end

% Write the header
fprintf(FID,'CCSDS_OEM_VERS = %s\n', EphemHead.CCSDS_OEM_VERS);
fprintf(FID,'CREATION_DATE = %s\n',  EphemHead.CREATION_DATE);
fprintf(FID,'ORIGINATOR = %s\n',     EphemHead.ORIGINATOR);

fprintf(FID,' \n');
fprintf(FID,'META_START\n');

fprintf(FID,'OBJECT_NAME = %s\n',    EphemHead.OBJECT_NAME);
fprintf(FID,'OBJECT_ID = %s\n',      EphemHead.OBJECT_ID);
fprintf(FID,'CENTER_NAME = %s\n',    EphemHead.CENTER_NAME);
fprintf(FID,'REF_FRAME = %s\n',      EphemHead.REF_FRAME);
fprintf(FID,'TIME_SYSTEM = %s\n',    EphemHead.TIME_SYSTEM);
fprintf(FID,'START_TIME = %s\n',     START_TIME);
fprintf(FID,'STOP_TIME = %s\n',      STOP_TIME);

fprintf(FID,'META_STOP\n');
fprintf(FID,' \n');

% Generate the epoch date strings
dsEpoch = cell(Ne,1);
for n=1:Ne
    ds = datestr(EphemData.Epoch(n),TimeFormat);
    dsEpoch{n} = ds;
end    

% Format for State vector lines
FMT = '%s';
for n=1:Ns
    FMT = cat(2,FMT,' %+0.15e');
end
FMT = cat(2,FMT,' \n');

% Write the State vector lines
if Ns == 6
    % Position/velocity states
    for n=1:Ne
        fprintf(FID,FMT, dsEpoch{n}, ...
                EphemData.State(n,1),EphemData.State(n,2),EphemData.State(n,3), ...
                EphemData.State(n,4),EphemData.State(n,5),EphemData.State(n,6));
    end    
else
    for n=1:Ne
        % Position/velocity/acceleration states
        fprintf(FID,FMT, dsEpoch{n}, ...
                EphemData.State(n,1),EphemData.State(n,2),EphemData.State(n,3), ...
                EphemData.State(n,4),EphemData.State(n,5),EphemData.State(n,6), ...
                EphemData.State(n,7),EphemData.State(n,8),EphemData.State(n,9));
    end    
end

% Write the Cov vector lines

if ~isempty(EphemData.Cov)
    
    % Size of covariance data
    sizeCov = size(EphemData.Cov);
    if sizeCov(1) ~= Ne
        error('Dimension mismatch between Epoch and Cov arrays');
    end
    
    % Get the dimension of the cov. matrix (Dc) from the number of elements
    % in the lower-triangular cov. vector
    Nc = sizeCov(2);
    if Nc == 36
        Dc = [6 6];
    else
        error('Invalid/unrecogniz3ed covariance dimension');
    end
    
    % Format strings
    FMTc = cell(6,1);
    FMTc{1} = '%+0.15e \n';
    FMTc{2} = '%+0.15e %+0.15e \n';
    FMTc{3} = '%+0.15e %+0.15e %+0.15e \n';
    FMTc{4} = '%+0.15e %+0.15e %+0.15e %+0.15e \n';
    FMTc{5} = '%+0.15e %+0.15e %+0.15e %+0.15e %+0.15e \n';
    FMTc{6} = '%+0.15e %+0.15e %+0.15e %+0.15e %+0.15e %+0.15e \n';
    
    % Write the covariance data
    fprintf(FID,' \n');
    fprintf(FID,'COVARIANCE_START\n');
    
    for n=1:Ne
        % ds = datestr(EphemData.Epoch(n),TimeFormat); ds(11) = 'T';
        fprintf(FID,'EPOCH = %s\n', dsEpoch{n});
        fprintf(FID,'COV_REF_FRAME = %s\n', EphemHead.REF_FRAME);
        cov = reshape(EphemData.Cov(n,:),Dc);
        % Lower triangle
        i=1; fprintf(FID,FMTc{i},cov(i,1));
        i=2; fprintf(FID,FMTc{i},cov(i,1),cov(i,2));
        i=3; fprintf(FID,FMTc{i},cov(i,1),cov(i,2),cov(i,3));
        i=4; fprintf(FID,FMTc{i},cov(i,1),cov(i,2),cov(i,3),cov(i,4));
        i=5; fprintf(FID,FMTc{i},cov(i,1),cov(i,2),cov(i,3),cov(i,4),cov(i,5));
        i=6; fprintf(FID,FMTc{i},cov(i,1),cov(i,2),cov(i,3),cov(i,4),cov(i,5),cov(i,6));
    end    

    fprintf(FID,'COVARIANCE_STOP\n');
    fprintf(FID,' \n');

end

% Close the output file
status = fclose(FID);
if status ~= 0
    fprintf('Error closing CCSDS OEM file: %s\n     %s\n', EphemFile);
else
    Success = true;
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2024-Jan-31 | Initial development, based on
%                |             | compatibility with the pre-existing
%                |             | function CCSDSParser.m, and focused on
%                |             | writing EME2000 frame OEM files.
% J. Halpin      | 2026-Apr-23 | Reformatted header to match code
%                |             | standards.