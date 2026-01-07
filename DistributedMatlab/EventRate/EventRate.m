function out = EventRate(paramsInfo)
%%
% EventRate - Estimate the rate conjunctions exceeding a threshold Pc value 
% occur for a model mission occupying an orbit/environment similar to that
% of a primary satellite with data archived in the CARA database.
%
% Syntax: out = EventRate(paramsInfo);
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: 
%
% Estimate the rate conjunctions exceeding a threshold Pc value occur for
% a model mission occupying an orbit/environment similar to that of a
% primary satellite or satellites with data archived in the CARA database.
%
% The analysis uses the method for semi-empirical conjunction risk analysis
% as described in
%
%  Doyle Hall (2019) "Determining Appropriate Risk Remediation Thresholds
%  from Empirical Conjunction Data using Survival Probability Methods" 
%  AAS 19-631.
%
% Specifically, the "last-update" Pc value for an event must exceed a
% user-specified Pc threshold to contribute to the estimated event rate.
% For a model mission that conducts remediation mitigation maneuvers
% (RMMs), last-update values are those estimated just before the
% mission's RMM "commit" time (but also after the mission's RMM "consider" 
% time).  The last-update Pc then can be regarded as the "perceived"
% collision risk of the event. The algorithm uses the archived record of
% last-update Pc values to estimate event rates, as well as the perceived
% mission-cumulative collision risk. 
%
% =========================================================================
%
% Input:
%   
%   paramsInfo - Input parameter structure containing information regarding
%                the input conjunction data, analysis configuration and the 
%                prospective mission to be analyzed. 
% 
%                See EventRate_default_params.m and 
%                EventRate_ConjDist_default_params.m for full 
%                documentation.
%                
%                Required fields:
%                   OCMDBFile OR PcTableFile - Name of an OCMDB data file 
%                   for processing. If OCMDBFile is used, EventRate will 
%                   precompute reference Pc values and save them to an
%                   augmented PcTableFile. This file may be provided as a
%                   PcTableFile parameter instead to skip preprocessing on
%                   future runs.
%
%                   priset - A set of one or more mission surrogate 
%                   satellite NORAD IDs. These IDs must be present in the 
%                   OCMDBFile or PcTableFile data
%
%                All other fields will be populated with default values if
%                not set
%
% =========================================================================
%
% Output:
%
%   out        -   output structure
%
%       CumulativePcNoRMMs   - [1x3] array of the [median, lower bound, 
%                              upper bound] of the 95% confidence interval 
%                              of estimated cumulative Pc for the 
%                              prospective mission duration without RMMs 
%                              considered 
%
%       CumulativePcWithRMMs - [1x3] array of the [median, lower bound, 
%                              upper bound] of the 95% confidence interval 
%                              of estimated cumulative Pc for the 
%                              prospective mission duration with RMMs
%                              considered
%   
%       MissionEventRate     - [1x3] array of the [median, lower bound,
%                              upper bound] of the 95% confidencer interval
%                              of estimated mission red event rate for the
%                              prospective mission red threshold Pc.
%
%       ConjDist             - Structure containing various intermediate
%                              values calculated within EventRate_ConjDist
%
% =========================================================================
%
% References:
%
%   Doyle Hall (2019) "Determining Appropriate Risk Remediation Thresholds
%   from Empirical Conjunction Data using Survival Probability Methods" 
%   AAS 19-631.
%
% =========================================================================
%
% Dependencies:
%
%   ProcessPcTable.m
%   EventRate_ConjDist.m
%   EventRate_default_params.m
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Disclaimer:
%
%    No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
%    WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
%    INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE
%    WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
%    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
%    INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
%    FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
%    THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
%    CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT
%    OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY
%    OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.
%    FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES
%    REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE,
%    AND DISTRIBUTES IT "AS IS."
%
%    Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
%    AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
%    SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF
%    THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES,
%    EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM
%    PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT
%    SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED
%    STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
%    PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE
%    REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL
%    TERMINATION OF THIS AGREEMENT.
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------
%% Initializations
persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'/params')); addpath(s.path);
    s = what(fullfile(p,'/src')); addpath(s.path);
    s = what(fullfile(p,'/src/matlab_lib')); addpath(s.path);
    s = what(fullfile(p,'/src/ConjDist')); addpath(s.path);
    s = what(fullfile(p,'../ProbabilityOfCollision')); addpath(s.path);
    s = what(fullfile(p,'../Utils/AugmentedMath')); addpath(s.path);
    s = what(fullfile(p,'../Utils/CovarianceTransformations')); addpath(s.path);
    s = what(fullfile(p,'../Utils/General')); addpath(s.path);
    s = what(fullfile(p,'../Utils/Plotting/Utils')); addpath(s.path);
    s = what(fullfile(p,'../Utils/TimeTransformations')); addpath(s.path);

    pathsAdded = true;
end

% Initialize the program execution parameters
if nargin < 1
    params = [];
else
    params = initialize_params(paramsInfo);
end

%% Augment params structure with the default EventRate parameters
params = EventRate_default_params(params);

%% Check for OCMDB file and associated Pc table file
if params.PcTable_mode == 0
    
    % Check for the existence of the OCMDB file
    if ~exist(params.OCMDBfile,'file')
        [~,f,e] = fileparts(params.OCMDBfile);
        warning(['Could not find OCMDB data file: ' [f e]]);
        out = [];
        return;
    end
    
    % Set PcTable to null to indicate Pc table mode is not being used
    params.PcTable = [];
    
else
    
    % Do the processesing required to use Pc table mode
    params = ProcessPcTable(params);
    
    % Ensure Pc table and DB exist
    if isempty(params.PcTable) || isempty(params.DB)
        warning('Invalid Pc table or DB array');
        out = [];
        return;
    end
    
    % Return now if only building the PcTable, and not doing any other
    % EventRate processing
    if params.PcTable_mode == 3
        out = [];
        return;
    end

    % Return with warning if using exclude_noncatastrophic parameter, which
    % is incompatible with PcTable mode
    if params.exclude_noncatastrophic
        warning('Pc table mode incompatible with exclude_noncatastrophic parameter');
        out = [];
        return;
    end
    
end

%% Perform the EventRate processing
% (using the EventRate version of the ConjDist function)
out = EventRate_ConjDist(params);

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
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================