function [PcTable,DB] = BuildPcTable(OCMDBfile,PcTableFile,params)
% BuildPcTable - Build an OCMDB file augmented with a Pc table for use by
% EventRate
%
% Syntax: [PcTable,DB] = BuildPcTable(OCMDBfile,PcTableFile,params)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Build an OCMDB file augmented with a Pc table, and with 
% sensitive quantities replaced with NaN values
%
% =========================================================================
%
% Input:
%
%   OCMDBfile   - 
%
%   PcTableFile - Optional - Filename under which to save the generated
%                  PcTable
%
%   Params      - Optional - Struct. See 'EventRate_default_params.m' for
%                  more documentation
%
% =========================================================================
%
% Output:
%
%   DB         - OCM DB table with filtered entries removed. See the file
%                'EventRate Truncated OCMDB Field Definitions.xlsx' in the 
%                'doc' directory for more details
%
%   PcTable    - Output structure containing parameters relevant to the
%                generation of reference Pc values in the output DB table:
%
%                 OCMDBfile - The name of the OCMDB file used to generate
%                             the PcTable
%
%                 HBRmin_meters - The minimum secondary HBR used for
%                                 precomputing reference Pc values
%
%                 HBRmax_meters - The maximum secondary HBR used for
%                                 precomputing reference Pc values
%
%                 HBRnum - The number of HBR values for which reference Pc 
%                          is precomputed
%
%                 HBR - Array of HBR values for precomputing reference Pc
%
%                 log10HBR - The log10 of values in HBR
%
%                 OriginalDBColumns - The number of columns in the output 
%                                     DB table
%
%                 ExpandedDBColumns - The total number of columns in the 
%                                     original DB table
%
%                 DeletedDBColumns  - The number of columns removed from 
%                                     the DB table
%
% =========================================================================
%
% Dependencies:
%
%   ../../ProbabilityOfCollision
%   ../../Utils/AugmentedMath
%   ../../Utils/CovarianceTransformations
%   ../../Utils/General
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------

% Defaults and initializations
Nargin = nargin;
if Nargin < 2
    PcTableFile = [];
end
if Nargin < 3
    params = [];
end

% Set Pc table parameters, if not already specified
params = set_default_param(params, 'PcTable_HBRmin_meters',10^(-2.0)); % 1 cm
params = set_default_param(params, 'PcTable_HBRmax_meters',10^( 2.0)); % 100 m
params = set_default_param(params, 'PcTable_HBRnum',41); % five points per logarithmic decade

% Ensure sensible HBR limits in Pc table
if params.PcTable_HBRmin_meters >= params.PcTable_HBRmax_meters
    error('Invalid Pc table HBR bounds');
end

% Ensure min number of points in Pc table
if params.PcTable_HBRnum < 7
    warning('Increasing PcTable to have at least seven points');
    params.PcTable_HBRnum = 7;
end

% Set up Pc table 
NHBR = params.PcTable_HBRnum;
log10HBRTable = linspace(log10(params.PcTable_HBRmin_meters), ...
                         log10(params.PcTable_HBRmax_meters), ...
                         NHBR); 
HBRTable = 10.^log10HBRTable;

% Check for the existence of the OCMDB file
[~,f,e] = fileparts(OCMDBfile);
if ~exist(OCMDBfile,'file')
    warning(['Could not find OCMDB data file: ' [f e]]);
    PcTable = []; DB = [];
    return;
end

% Load the OCMDB file creating the DB array
disp([current_timestring() ' Building OCMDB PcTable file']);
disp([' Loading OCMDB file ' f e]);
DB = load(OCMDBfile); DB = DB.DB;

% Number of columns in DB
Ncol = size(DB,2);
if Ncol < 218
    error('Insufficient number of OCMDB columns (< 218)');
elseif Ncol > 234
    error('Unexpectedly large number of OCMDB columns (> 254)');
end

% Eliminate bad DB entries
DB = EliminateBadDBEntries(DB,1,1,1,0);
Nconj = size(DB,1); Nconjstr = num2str(Nconj);

% Build the covariance matrix blocks
CovBuild6 = triu(ones(6,6))';
CovBuild6(CovBuild6~=0)  = 1:21;
CovBuild6 = CovBuild6 + triu(CovBuild6',1);

% Radius of the Earth in meters
Rearth = 6378.137e3;

% Degrees per radian
DegPerRad = 180/pi;

% Allocate Pc table array and the altitude/latitude table array
PcArray = NaN(Nconj,NHBR);
ABArray = NaN(Nconj,2);

% Expand HBR table for parfor operations
HBR = repmat(HBRTable,[Nconj 1]);

% Extract states and covariances (m units)
DBX1 = DB(:,172:177)*1e3;
DBC1 = DB(:,73:93);
DBX2 = DB(:,178:183)*1e3;
DBC2 = DB(:,133:153);

% Set up PcMultiStep parameters
PcMultiStep_params = [];
PcMultiStep_params.OnlyPc2DCalculation = true;

% Set up progress bar
hbar = parfor_progressbar(Nconj,'','Name','Pc Table Calculation'); 
verbose = false;

% Loop to create Pc table
parfor n=1:Nconj % Parfor enabled
    
    if verbose
        disp(['  Processing conjunction ' num2str(n) ' of ' Nconjstr]);
    end
    
    % Primary object ECI state
    Xp = DBX1(n,:);
    
    % Primary object covariance 
    TriEls = DBC1(n,:);
    C1UVW = TriEls(CovBuild6); % [m^2 m^2/s m^2/s^2]
    
    % Convert UVW (same as the RIC) to ECI
    Pp = RIC2ECI(C1UVW,Xp(1:3),Xp(4:6));  % [m^2 m^2/s m^2/s^2]
    Pp = cov_make_symmetric(Pp);
    
    % Secondary object ECI state
    Xs = DBX2(n,:);

    % Secondary object covariance 
    TriEls = DBC2(n,:);
    C2UVW = TriEls(CovBuild6);  % [m^2 m^2/s m^2/s^2]
    
    % Convert UVW (same as the RIC) to ECI
    Ps = RIC2ECI(C2UVW,Xs(1:3),Xs(4:6));  % [m^2 m^2/s m^2/s^2]
    Ps = cov_make_symmetric(Ps);
    
    % Allocate Pc buffer
    PcBuffer = NaN(1,NHBR);

    % Calculate Pc values
    for i=1:NHBR
        PcBuffer(i) = PcMultiStep( ...
            Xp(1:3),Xp(4:6),Pp, ...
            Xs(1:3),Xs(4:6),Ps, ...
            HBR(n,i),PcMultiStep_params);
    end

    % Store buffered Pc values into array
    PcArray(n,:) = PcBuffer';
    
    % Calculate the altitude and latitude of the primary at TCA
    
    % ECI coordinates for the conjunctions in this event
    X = Xp(1);
    Y = Xp(2);
    Z = Xp(3);

    % Radius (R), altitude (A), and latitude (B).
    R = sqrt(X.^2+Y.^2+Z.^2);
    A = R-Rearth;
    B = asin(Z./R)*DegPerRad;
    B(B < -90) = -90; B(B > +90) = +90;
    
    % Store altitude and latitude
    ABArray(n,:) = [A B];
    
    hbar.iterate(1);

end

% Close progress bar
close(hbar);

% Clear temporary variables
clear DBX1 DBC1 DBX2 DBC2 HBR;

% Add tables to DB array, and calculate expanded number of columns
DB = [DB PcArray ABArray];
NcolExpanded = size(DB,2);

% OCMDB file columns to delete when building the PcTable, which hold
% unnecessary information (or unreleasable information)
DelColumns = [10 19:36 37:43 44:50 54:96 97:103 104:110 114:164 167:168 170:171 172:234];

% Delete DB array columns holding unnecessary information
ndx = DelColumns <= Ncol;
DelColumns = DelColumns(ndx);
keep = true([1 NcolExpanded]);
keep(DelColumns) = false;
DBTable = DB(:,keep); % Contracted DB to be saved into a matfile later

% Define the PcTable structure
% PcTable.OCMDBfile = OCMDBfile;
PcTable.HBRmin_meters = params.PcTable_HBRmin_meters;
PcTable.HBRmax_meters = params.PcTable_HBRmax_meters;
PcTable.HBRnum = params.PcTable_HBRnum;
PcTable.HBR = HBRTable;
PcTable.log10HBR = log10HBRTable;
PcTable.OriginalDBColumns = Ncol;
PcTable.ExpandedDBColumns = NcolExpanded;
PcTable.DeletedDBColumns  = DelColumns;

% Save if required
if ~isempty(PcTableFile)
    [~,f,e] = fileparts(PcTableFile);
    disp([current_timestring() ' Saving OCMDB PcTable file ' f e]);    
    save(PcTableFile,'PcTable','DBTable');
end

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
% D. Hall   | 2022-Apr-11 | Initial version.
% N. Ravago | 2025-Oct-29 | Integrated PcMultiStep
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================