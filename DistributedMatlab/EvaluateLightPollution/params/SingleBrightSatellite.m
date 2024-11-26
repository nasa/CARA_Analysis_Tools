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
params.Output.outputpath = fullfile('output','SingleBrightSatellite');

%% Parameters of the new or proposed constellation to evaluate

% Number of satellites in constellation
params.New.Nc = 1;          

% Altitude of new constellation
params.New.Altitude_km = 540;

% Inclination of new constellation (degrees)
params.New.Inclination_deg = 28.5;

% Brightness parameters
params.New.FracBrighterThanSatCon1Limit = 1; % NaN to calculate from analog satellites

%% Zenith range magnitude statistics (i.e., normalized to altitude range)
% Make these magnitudes so small, so that the satellite registers as 
% brighter than recommended essentially under all circumstances
params.New.Mzen50    = 6.12 - 50;
params.New.Mzen05    = 6.67 - 50;
params.New.Mzen95    = 4.68 - 50;
params.New.MzenAltkm = 1000;

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