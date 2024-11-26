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
params.Output.outputpath = fullfile('output','OneWebPhase2_ReducedBrightness');

%% Parameters of the new or proposed constellation to evaluate

% Number and altitudes of satellites in each shell of new constellation
params.New.Nc              = [1764 2304 2304];
params.New.Altitude_km     = [1200 1200 1200];
params.New.Inclination_deg = [87.9 40.0 55.0];

%% Parameters for the analog satellite(s)
params.Analog.Type         = 'DifferentSatelliteDesign';
params.Analog.Altitude_km  = 1200; % OneWeb altitude
params.Analog.datapath     = 'data/MMT/OneWeb'; % MMT data path
params.Analog.UTbegin      = '2021-01-01 00:00:00.000'; % MMT data start
params.Analog.UTend        = '2021-03-01 00:00:00.000'; % MMT data end

% OneWeb 1.0m x 1.3m nadir face with two 0.6m diameter solar dishes
params.Analog.ReflBoxHWL_m = 1.9; % = 1.0*1.3 + 2*(pi*0.3^2);

% Brightness reduction of constellation satellite design relative to the
% analog satellite design, in stellar magnitudes
Mreduce = 1.7;

% Photometric flux reduction factor
Freduce = 10^(-0.4*Mreduce);

% Scaled nadir-facing area (assuming constant albedo)
params.New.ReflBoxHWL_m = Freduce * params.Analog.ReflBoxHWL_m;

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