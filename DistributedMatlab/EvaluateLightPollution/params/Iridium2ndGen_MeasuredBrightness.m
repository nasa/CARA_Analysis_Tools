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
params.Output.outputpath = fullfile('output','Iridium2ndGen_MeasuredBrightness');

%% Parameters of the new or proposed constellation to evaluate

% Number of satellites in constellation
params.New.Nc = 75;

% Altitude of new constellation
params.New.Altitude_km = 780;

% Inclination of new constellation (degrees)
params.New.Inclination_deg = 86.4;

%% Parameters for the analog satellite(s)

params.Analog.Type        = 'SameSatelliteDesign';
params.Analog.Altitude_km = 780;
params.Analog.datapath    = 'data/MMT/Iridium_2ndGen'; % MMT data path
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