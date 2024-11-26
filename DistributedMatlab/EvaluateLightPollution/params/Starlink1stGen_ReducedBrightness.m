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
params.Output.outputpath = fullfile('output','Starlink1stGen_ReducedBrightness');

%% Parameters of the new or proposed constellation to evaluate

% Number and altitudes of satellites in each shell of new constellation
params.New.Nc              = [1584 1584  720  348  172 2493 2478 2547];
params.New.Altitude_km     = [ 550  540  570  560  560  336  341  346];
params.New.Inclination_deg = [53.0 53.2 70.0 97.6 97.6 42.0 48.0 53.0];

%% Parameters for the analog satellite(s)
params.Analog.Type         = 'DifferentSatelliteDesign';
params.Analog.Altitude_km  = 550; % Starlink altitude
params.Analog.ReflBoxHWL_m = [0.1 1.6 3.2]; % Box enclosing reflective components
params.Analog.datapath     = 'data/MMT/Starlink_VisorSat'; % MMT data path
params.Analog.UTbegin      = '2021-01-01 00:00:00.000'; % MMT data start
params.Analog.UTend        = '2021-03-01 00:00:00.000'; % MMT data end

% Photometric flux scaling factor
FluxScaleFactor            = 0.1;
SizeScaleFactor            = sqrt(FluxScaleFactor);
params.New.ReflBoxHWL_m    = SizeScaleFactor * params.Analog.ReflBoxHWL_m;

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