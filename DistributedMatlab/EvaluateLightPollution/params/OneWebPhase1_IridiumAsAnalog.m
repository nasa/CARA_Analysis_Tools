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
params.Output.outputpath = fullfile('output','OneWebPhase1_IridiumAsAnalog');

%% Parameters of the new or proposed constellation to evaluate

% Number of satellites in constellation
params.New.Nc = 1980;          
% Altitude of new constellation
params.New.Altitude_km = 1200;
% Inclination of new constellation (degrees)
params.New.Inclination_deg = 87.9;

% OneWeb 1.0m x 1.3m nadir face area with two 0.6m diameter solar dishes
% (single nadir-area in m^2 supplied as a scalar parameter)
params.New.ReflBoxHWL_m = 1.9; % = 1.0*1.3 + 2*(pi*0.3^2);

%% Parameters for the analog satellite(s)

params.Analog.Type         = 'DifferentSatelliteDesign';
params.Analog.Altitude_km  = 780; % Iridium altitude

params.Analog.datapath     = 'data/MMT/Iridium_2ndGen'; % MMT data path
params.Analog.UTbegin      = '2021-01-01 00:00:00.000'; % MMT data start
params.Analog.UTend        = '2021-03-01 00:00:00.000'; % MMT data end

% Iridium 2nd Gen bus encl. box
% (bounding box specified as 3-D array)
% params.Analog.ReflBoxHWL_m = [1.5 2.4 3.1]; 

% Iridium 2nd Gen bus 3.1m x 2.4m nadir face
% (single nadir-area in m^2 supplied as a scalar parameter)
params.Analog.ReflBoxHWL_m = 7.4; % = 3.1*2.4;

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