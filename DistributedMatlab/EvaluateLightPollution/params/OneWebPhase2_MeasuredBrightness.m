% Function EvaluateLightPollution parameter specification file
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

%% Output path
params.Output.outputpath = fullfile('output','OneWebPhase2_MeasuredBrightness');

%% Parameters of the new or proposed constellation to evaluate

% Number and altitudes of satellites in each shell of new constellation
params.New.Nc              = [1764 2304 2304];
params.New.Altitude_km     = [1200 1200 1200];
params.New.Inclination_deg = [87.9 40.0 55.0];

%% Evaluation ruleset
params.Evaluation.RuleSet = 1;
params.Plotting.NbrightShadingLevel = 2;

%% Parameters for the analog satellite(s)
params.Analog.Type        = 'SameSatelliteDesign';
params.Analog.Altitude_km = 1200;
params.Analog.datapath    = 'data\MMT\OneWeb'; % MMT data path
params.Analog.UTbegin     = '2021-01-01 00:00:00.000'; % MMT data start
params.Analog.UTend       = '2021-03-01 00:00:00.000'; % MMT data end

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