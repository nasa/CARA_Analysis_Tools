function params = EventRate_default_params(params)
% EventRate_default_params - Set the default parameters for function 
%                            EventRate
%
% Syntax: params = EventRate_default_params(params)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Set the default parameters for function EventRate
%
% =========================================================================
%
% Input:
%
%   params - Struct containing EventRate parameters. Un-initialized
%            parameters will be set to default values. Preset parameters 
%            will not be overwritten
%
% =========================================================================
%
% Output:
%
%   params - Struct containing EventRate parameters
%
%       EventRate_estimation_mode - EventRate estimation mode:
%               1 => Estimates cumulative mission Pc and RMM rates, given
%                    a RMM Pc threshold (default)
%               2 => Estimates required RMM Pc threshold, given a lifetime 
%                    cumulative Pc goal (mostly obsolete analysis mode)
%
%       EventRate_Pc_value        - EventRate Pc value - specifies fixed 
%                                   RMM Pc threshold for mode 1; specifies 
%                                   cumulative mission Pc goal mode 2)
%                                   (default = 1.0E-4)
%
%       OCMDBfile                 - OCMDB file name (set to a test file by 
%                                   default)
%
%       PcTable_mode              - Integer indicating the mode to build 
%                                   and/or use PcTable mat file:
%               0 => Do not build or use a PcTable mat file
%                   (as programmed in first version)
%               1 => Build and use a PcTable mat file [default]
%                   (reuse old tables, if found)
%               2 => Build and use a PcTable mat file
%                   (do not reuse old tables)
%               3 => Build a PcTable mat file only
%                   (do not reuse old tables)
%
%       PcTable_HBRmin_meters     - Minimum HBR for secondary objects
%                                   (default = 10^-2 meters)
%
%       PcTable_HBRmax_meters     - Maximum HBR for secondary objects
%                                   (default = 10^2 meters)
%
%       PcTable_HBRnum            - Number of points to use when building
%                                   the PcTable
%                                   (default = 41, 10 points per 
%                                   logarithmic decade based on default 
%                                   min/max values)
%
%       RMM_execution_rates       - Bool indicating whether to estimate RMM
%                                   execution or planning rates. 'true' for
%                                   execution rates, 'false' for planning 
%                                   rates.
%                                   (default = true)
%
%       redyel_event_option       - String indicating whether to classify
%                                   red/yellow events based on maximum Pc 
%                                   update during the consider/commit time 
%                                   window ('mud') or the last Pc update 
%                                   ('lud'). Set automatically based on the
%                                   RMM_execution_rates mode.
%                                   (default = 'lud')
%
%       priset                    - Array of integers indicating the 
%                                   primary IDs to combine for 
%                                   semi-empirical conjunctions 
%                                   (default = _____)
%
%       prisetstr                 - Cell array of strings indicating the 
%                                   primary satellite names to combine for 
%                                   semi-empirical conjunctions. Can be set
%                                   to empty to initialize based on values 
%                                   in priset. 
%                                   (default = '')
%
%       mission_name              - New mission name 
%                                   (default = 'NewMission')
%
%       mission_HBR_meters        - New mission HBR in meters 
%                                   (default = 10)
%
%       mission_duration_years    - New mission duration in years 
%                                   (default = 15)
%
%       mission_on_orbit_mass_kg  - New mission on-orbit mass in kg
%                                   (default = 2000)
%
%       Tcommit_days              - RMM commit time in days before TCA
%                                   (default = 1.0)
%
%       Tconsider_days            - RMM consider time in days before TCA
%                                   (default = 3.5)
%
%       exclude_noncatastrophic   - Bool indicating whether to exclude 
%                                   likely non-catastrophic events from 
%                                   analysis
%                                   (default = false)
%
%       SecondaryCatalogGrowth    - Parameter used to model secondary
%                                   catalog growth. A factor of 1.3 
%                                   indicates that the prospecive mission 
%                                   encounters a future secondary satellite 
%                                   population 30% larger than in the 
%                                   archived data set.
%                                   (default = 1.0 (no growth))
%
%       report_and_plot_svi_rates - Mode indicating whether to report and
%                                   plot Screening Volume Incursion (SVI) 
%                                   rates 
%                                   (default = false)
%
%       make_time_bins            - Bool used in conjunction with SVI 
%                                   plotting mode 
%                                   (default = false)
%
%       make_time_plots           - Array of bools used in conjunction with 
%                                   SVI plotting mode 
%                                   (default = [false false false])
%
%       plot_spatial_distribution - Bool indicating whether to plot alt-lat
%                                   conjunction plots 
%                                   (default = true)
%
%       make_update_plots         - Bool indicating whether to make Pc
%                                   update sequence plots 
%                                   (default = false)
%
%       make_Pcum_plots           - Bool indicating whether to make Pcum
%                                   plots. Set automatically based on 
%                                   EventRate_estimation_mode. 
%                                   (default = RMM_execution_rates)
%
%       make_Pcum_images          - Bool indicating whether to make Pcum
%                                   images. Set automatically based on 
%                                   EventRate_estimation_mode. 
%                                   (default = false)
%
%       visible_plots             - Bool indicating whether plots should be
%                                   visible on screen. Plots are saved
%                                   regardless of visibility option.
%                                   (default = true)
% =========================================================================
%
% Dependencies:
%
%   set_default_param
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

% Set the default parameters for function EventRate

params = set_default_param(params, 'EventRate_estimation_mode', 1);
if params.EventRate_estimation_mode ~= 1 && ...
   params.EventRate_estimation_mode ~= 2
    error('Invalid EventRate_estimation_mode parameter');
end

% EventRate Pc value (specifies fixed RMM Pc threshold for mode 1;
%                     specifies cumulative mission Pc goal mode 2)
params = set_default_param(params, 'EventRate_Pc_value', 1e-4);

% Default OCMDB file
params = set_default_param(params, 'OCMDBfile', '');
params = set_default_param(params, 'PcTableFile', '');
params = set_default_param(params, 'PcTable_mode', 1);

% Pc table parameters
params = set_default_param(params, 'PcTable_HBRmin_meters',10^(-2.0)); % 1 cm
params = set_default_param(params, 'PcTable_HBRmax_meters',10^( 2.0)); % 100 m
params = set_default_param(params, 'PcTable_HBRnum',41); % five points per logarithmic decade

% EventRate mode to model RMM execution vs RMM planning rates
params = set_default_param(params,'RMM_execution_rates',true);
if ~isequal(params.RMM_execution_rates,true) && ...
   ~isequal(params.RMM_execution_rates,false)
    error('Invalid RMM_execution_rates parameter');
end
if ~params.RMM_execution_rates && params.EventRate_estimation_mode ~= 1
    error('RMM planning rates estimation mode requires EventRate_estimation_mode = 1');
end

% Set the red/yellow event option
if params.EventRate_estimation_mode == 1
    if params.RMM_execution_rates
        % Last update option within commit/consider time limits is 
        % used to model RMM execution rates
        params.redyel_event_option = 'lud';
    else
        % Max-Pc update option within commit/consider time limits is 
        % used to model RMM planning rates
        params.redyel_event_option = 'mud';
    end
end

% Set of CARA primaries to combine for semi-empirical conjunctions
params = set_default_param(params, 'priset', 20580);
params = set_default_param(params, 'prisetstr', '');
if isempty(params.prisetstr)
    if (numel(params.priset) <= 2)
        params.prisetstr = num2str(params.priset,'%05i ');
        params.prisetstr = strrep(params.prisetstr,' ','_');
    else
        params.prisetstr = [num2str(params.priset(1),'%05i ') ...
            '_' num2str(numel(params.priset)-1) 'more'];
    end
end

% Mission name
params = set_default_param(params, 'mission_name', 'NewMission');

% Mission HBR and duration
params = set_default_param(params, 'mission_HBR_meters', 10);
params = set_default_param(params, 'mission_duration_years', 15);
params = set_default_param(params, 'mission_on_orbit_mass_kg', 2000);

% Mission RMM commit and consider times (days from TCA)
params = set_default_param(params, 'Tcommit_days',   1.0);
params = set_default_param(params, 'Tconsider_days', 3.5);

% Exclude likely non-catastrophic events from analysis
params = set_default_param(params, 'exclude_noncatastrophic', false);

% Secondary catalog population growth factor
params = set_default_param(params, 'SecondaryCatalogGrowth', 1);

% Screening volume parameters
params = set_default_param(params, 'report_and_plot_svi_rates', false);
if params.report_and_plot_svi_rates
    warning('EventRate_ConjDist does not support SVI analysis (use ConjDist)');
    if ~isfield(params,'SVRICtype'  ) || isempty(params.SVRICtype  ) || ...
       ~isfield(params,'SVRICdims_m') || isempty(params.SVRICdims_m)
        error('If reporting SVI rates, {priorb,SVRICtype,SVRICdims_m} must be defined');
    end
    params.make_time_bins = true;
    params.make_time_plots = [true true true];
else
    params.make_time_bins = false;
    params.make_time_plots = [false false false];
end

% Plotting flags

% Make conjunction latitude vs altitude spatial distribution plots
params = set_default_param(params, 'plot_spatial_distribution', true);

% Make Pc update-sequence plots
params = set_default_param(params, 'make_update_plots', false);

% Flags to make Pcum plots and images
if params.EventRate_estimation_mode == 1
    % Default images for mode 1
    params = set_default_param(params, 'make_Pcum_plots',  params.RMM_execution_rates);
    params = set_default_param(params, 'make_Pcum_images', false);
else
    % Default images for mode 2
    params = set_default_param(params, 'make_Pcum_plots',  true);    
    params = set_default_param(params, 'make_Pcum_images', true);
end

% Make plots visible on the screen
params = set_default_param(params, 'visible_plots', true);

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
% N. Ravago | 2025-Oct-29 | Removed Pc function indicators
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================