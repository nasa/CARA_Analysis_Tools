% Function EvaluateLightPollution parameter specification file for the all
% eight orbital shells of the of the Starlink 1st generation constellation
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
params.Output.outputpath = fullfile('output','Starlink1stGen_MeasuredBrightness');

%% Parameters of the new or proposed constellation to evaluate

% Number, altitudes and inclinations of satellites in each shell of
% the new or proposed constellation

params.New.Nc              = [1584 1584  720  348  172 2493 2478 2547];
params.New.Altitude_km     = [ 550  540  570  560  560  336  341  346];
params.New.Inclination_deg = [53.0 53.2 70.0 97.6 97.6 42.0 48.0 53.0];

%% Parameters for the analog satellite(s)

params.Analog.Type        = 'SameSatelliteDesign';
params.Analog.Altitude_km = 550;
params.Analog.datapath    = 'data/MMT/Starlink_VisorSat'; % MMT data path
params.Analog.UTbegin     = '2021-01-01 00:00:00.000';    % MMT data start
params.Analog.UTend       = '2021-03-01 00:00:00.000';    % MMT data end

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