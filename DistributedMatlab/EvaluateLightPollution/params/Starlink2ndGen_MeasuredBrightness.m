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
params.Output.outputpath = fullfile('output','Starlink2ndGen_MeasuredBrightness');

%% Number and altitudes of satellites in each shell of new constellation
params.New.Nc              = [7178 7178 7178 2000 1998 4000   144   324];
params.New.Altitude_km     = [ 328  334  345  360  373  499   604   614];
params.New.Inclination_deg = [30.0 40.0 53.0 96.9 75.0 53.0 148.0 115.7];

%% Parameters for the analog satellite(s)
params.Analog.Type        = 'SameSatelliteDesign';
params.Analog.Altitude_km = 550;
params.Analog.datapath    = 'data/MMT/Starlink_VisorSat'; % MMT data path
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