function params = EventRate_ConjDist_default_params(params)
% EventRate_ConjDist_default_params - Set remaining default parameters for 
%                                     function EventRate
%
% Syntax: params = EventRate_ConjDist_default_params(params)
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Set remaining default parameters for function 
%              EventRate_ConjDist
%
% =========================================================================
%
% Input:
%
%   params - Struct containing EventRate parameters. Un-initialized
%            parameters will be set to default values. Preset parameters 
%            will not be overwritten. 
%
% =========================================================================
%
% Output:
%
%   params - Struct containing EventRate_ConjDist parameters. See
%            EventRate_default_params.m for parameters undocumented here.
%
%       start_time          - time at start of simulation
%
%       OCMDBpath           - path to OCMDB file
%
%       OCMDBroot           - OCMDB filename
%
%       OCMDBext            - OCMDB file extension
%
%       DB                  - Contents of a OCMDB file
%                             (default = empty)
%
%       secset              - Array of secondary IDs to consider
%                             (default = empty (consider all secondaries))
%
%       redyel_event_HBR_m  - Prospective mission HBR in meters (obsolete, 
%                             EventRate_ConjDist uses mission_HBR_meters 
%                             instead)
%
%       CommitConsiderHours - Array of the [commit, consider] time    [1x2]
%                             window limits in hours. Only updates
%                             within this window may be considered 
%                             in calculations.
%                             (default = [1.5, 7.0]*24)
%
%       postTCA_limit_hours - Maximum time in hours between last update and 
%                             OCMDB entry creation to plot in update 
%                             sequence plots
%                             (default = max(CommitConsiderHours,168))
%
%       Tmissionyrs         - Prospective mission duration in years
%
%       Tmissionyrsplot     - Mission durations to plot. Currently only
%                             supporting 1 duration
%
%       RemManeuver         - String indicating RMM type: 'Translational'
%                             or 'Rotational'
%                             (default = 'Translational')
%
%       RemReduction        - RMM reduction factor
%                             (default = 0.03 for 'Translational',
%                                        0.20 for 'Rotational)
%
%       runtag              - Output file tag (string) containing OCMDB
%                             filename and primary IDs
%
%       timetag             - Simulation start time timetag in string form; 
%                             derived from params.start_time
%
%       timestamp           - Bool indicating whether to include the the
%                             timetag parameter in the output filename
%
%       outtag -            - output file tag combining runtag and timetag
%                             if timestamp is 'true'
%
%       outputpath          - output folder
%
%       outputdir           - output directory within output folder
%                             (contains outtag to keep distinct runs 
%                             separated)
%
%       logging             - Bool indicating whether to maintain a log
%                             file
%                             (default = true)
%
%       displaying          - Bool indicating whether to display logged
%                             events in the console
%                             (default = true)
%
%       log_file_name       - filename for log file
%
%       logfid              - MATLAB file ID for log file
%
%       SVIRICtype_alt      - String indicating SV shape (alternative name, 
%                             derived from SVRICtype)
%
%       SVRICdims_m         - Array containing SV dimensions in  [Npri x 3]
%                             meters. If [1x3] input is provided 
%                             but Npri > 1, will be automatically
%                             replicated to [Npri x 3]
%                             (default = [2 44 51]*1000)
%
%       TCAtol_s            - Time window in seconds used to associate 
%                             updatesinvolving the same primary/secondary. 
%                             Updates with TCAs within this window will be 
%                             associated to the same event
%                             (default = 15*60 = 900)
%
%       TCRtol_s            - (obselete, currently unused)
%                             (default = 15)
%
%       TCAlimits_UT        - cell array of Minimum/maximum TCA times {1x2}
%                             to consider. Strings are dates in the
%                             format 'yyyy-mm-dd HH:MM:SS.FFF'
%
%       TCAbin_days                  - TCA bin width in days
%                                      (default = 7)
%
%       Pcred, log10_Pcred           - Pc threshold for red events
%                                      (default = 1.0e-4)
%
%       Pcyel, log10_Pcyel           - Pc threshold for yellow events
%                                      (default = 1.0e-7)
%        
%       Pcgre, log10_Pcgre           - Pc threshold for green events
%                                      (default = 1.0e-10)
%
%       acol                         - RGB triplet to use for all     [1x3]
%                                      events in ConjDist_time_plot
%
%       ycol                         - RGB triplet to use for yellow  [1x3]
%                                      events in ConjDist_time_plot
%
%       rcol                         - RGB triplet to use for red     [1x3]
%                                      events in ConjDist_time_plot
%
%       gcol                         - RGB triplet to use for green   [1x3]
%                                      events in ConjDist_time_plot
%
%       rng_mode                     - String indicating RNG mode: 
%                                      'default' or 'shuffle'. 
%                                      (default = 'default)
%
%       association_check            - Check for disagreements within OCMDB
%                                      file associations between  
%                                      conjunctions and events. Set to 
%                                      false automatically when using a
%                                      PcTable file.
%                                      (default = false)
%
%       make_event_file              - Bool indicating whether to create a 
%                                      .csv output file listing all 
%                                      associated events
%                                      (default = false)
%
%       make_repeating_event_file    - Bool indicating whether to create a 
%                                      repeating event file (not currently
%                                      supported)
%                                      (default = false)
%
%       repeating_event_span_days    - Repeating event logging currently 
%                                      not supported
%                                      (default = 7)
%
%       repeating_event_Pc_threshold - Repeating event logging currently 
%                                      not supported
%                                      (default = 1.0e-4)
%
%       plot_err_bars                - Array of bools indicating      [1x2]
%                                      whether to plot error bars. 
%                                      First element controls 
%                                      vertical error bars, second 
%                                      controls horizontal. If prior 
%                                      value is a single bool, will 
%                                      be replicated to [1x2]
%                                      (default = [false, true])
%
%       plot_conf_ranges             - Array of bools indicating      [1x2]
%                                      whether to plot confidence 
%                                      ranges in ConfDist_time_panel. 
%                                      First element controls unified 
%                                      bin-set conf range, second
%                                      controls bin-to-bin conf range
%                                      (default = [false, false])
%
%       Pc_cutoff_accum_risk         - Minimum risk for observed events to 
%                                      be used in calculating accumulated
%                                      risk
%                                      (default = params.Pcgre)
%
%       N_botstrap_accum_risk        - Number of particles to use for 
%                                      calculating accumulated risk
%                                      (default = 1000)
%
%       lnPc_boostrap_uncertainty    - Factor of uncertainty to assume in
%                                      observed event Pc values
%                                      (default = 0.3)
%
%       add_one_unrem_event          - Bool indicating whether to add one
%                                      unremediated event for conservatism
%                                      (default = true)
%
%       confidence_noncatastrophic   - Minimum non-catastrophic chance for
%                                      an event to be considered 
%                                      noncatastrophic and be disregarded 
%                                      in risk calculation
%                                      (default = 0.999)
%
%       FredVsPremImage              - Bool indicating whether to make a 
%                                      Fred vs. Prem image
%                                      (default = false)
%
%       Pc_for_update_plots          - If generating update sequence plots,
%                                      minimum Pc at last/maximum update 
%                                      (depending on mode) for which
%                                      to generate a plot.
%                                      (default = params.Pcred)
%
%       Nu_for_update_plots          - If generating update sequence plots,
%                                      minimum number of updates in the 
%                                      sequence for which to generate a 
%                                      plot
%                                      (default = 0)
%
%       make_altlat_plots            - Bool indicating whether or not to 
%                                      make altitude-latitude plots
%                                      (default = false)
%
%       plot_format                  - File extension for saved plot images
%                                      (default = '.png')
%
%       make_first_generation_plots  - Bool indicating whether to perform
%                                      first-generation remediation 
%                                      threshold analysis + generate plots
%                                      (default = false)
%
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
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

Nargin = nargin;
if (Nargin < 1); params =[]; end

% Note start time
params.start_time = current_timestring();

% Default OCMDB file
% if ~isfield(params,'OCMDBfile') || isempty(params.OCMDBfile)
%     error('The file name OCMDBfile must be specified');
% end
% [params.OCMDBpath, params.OCMDBroot, params.OCMDBext] = ...
%     fileparts(params.OCMDBfile);

% Check if DB array for the OCMDB file is provided as a parameter
if ~isfield(params,'DB'); params.DB = []; end

% Column vector of primary object numbers [Npri x 1]
if ~isfield(params,'priset') || isempty(params.priset)
    error('The set of primaries priset  must be specified');
end

if ~isfield(params,'secset'); params.secset = []; end

params.priset = sort(params.priset);
params.Npri = numel(params.priset);

if ~isfield(params,'prisetstr'); params.prisetstr = ''; end
    
if isempty(params.prisetstr)
    if (params.Npri <= 6)
        params.prisetstr = num2str(params.priset,'%05i ');
        params.prisetstr = strrep(params.prisetstr,' ','_');
    else
        params.prisetstr = [num2str(params.priset(1),'%05i ') ...
            '_among_' num2str(params.Npri)];
    end
end

% Options for registering red and yellow events
if ~isfield(params,'redyel_event_option') || isempty(params.redyel_event_option)
    if params.RMM_execution_rates
        % Modeling RMM execution rates requires the use of the 
        % last update event Pc within commit/consider/catastrophic limits
        params.redyel_event_option = 'lud';
    else
        % Modeling RMM planning rates requires the use of the 
        % maximum update Pc within commit/consider/catastrophic limits
        params.redyel_event_option = 'mud';
    end
end
params.redyel_event_option = lower(params.redyel_event_option(1:3));

switch params.redyel_event_option
    case {'lud', 'mud'}
        % Valid options for RMM execution and planning rate analyses
    case {'max', 'min' 'med', 'fud'}
        % Valid options, but not used for RMM rate analysis
        warning(['redyel_event_option = ' params.redyel_event_option ...
            ' is a non-standard calculation mode and not recommended']);
    otherwise
        error('Invalid redyel_event_option parameter');
end

if ~isfield(params,'redyel_event_HBR_m')
    params.redyel_event_HBR_m = [];
elseif isnan(params.redyel_event_HBR_m)
    params.redyel_event_HBR_m = [];
end

% Last update commit and consider times for maneuver execution or planning
if ~isfield(params,'CommitConsiderHours'); params.CommitConsiderHours = []; end
if isempty(params.CommitConsiderHours)
    params.CommitConsiderHours = [1.5 7.0]*24;
end
if any(params.CommitConsiderHours < 0)
    error('Last update parameter CommitConsiderHours must be nonnegative');
elseif params.CommitConsiderHours(1) == params.CommitConsiderHours(2)
    error('Last update parameter CommitConsiderHours span an interval, e.g., [72 168]');
end
params.CommitConsiderDays = params.CommitConsiderHours/24;

% Last update post-TCA limit
if ~isfield(params,'postTCA_limit_hours'); params.postTCA_limit_hours = []; end
if isempty(params.postTCA_limit_hours)
    params.postTCA_limit_hours = max(max(abs(params.CommitConsiderHours)),7*24);
end
if numel(params.postTCA_limit_hours) > 1 || ...
   ~isequal(class(params.postTCA_limit_hours),'double')
    error('Invalid postTCA_limit_hours parameter');
elseif params.postTCA_limit_hours < 0
    warning('The postTCA_limit_hours parameter cannot be negative; setting to zero');
    params.postTCA_limit_hours = 0;
end
params.postTCA_limit_days = params.postTCA_limit_hours/24;

% Mission Pcum limit duration for threshold Pc estimation
if ~isfield(params,'Pcmission'); params.Pcmission = []; end
if isempty(params.Pcmission)
    params.Pcmission = 1e-3;
end

% Mission duration
if ~isfield(params,'Tmissionyrs'); params.Tmissionyrs = []; end
if isempty(params.Tmissionyrs)
    params.Tmissionyrs = 10;
end

% Mission durations to plot
if ~isfield(params,'Tmissionyrsplot'); params.Tmissionyrsplot = []; end

% Mission secondary catalog population growth factor
if ~isfield(params,'SecondaryCatalogGrowth'); params.SecondaryCatalogGrowth = []; end
if isempty(params.SecondaryCatalogGrowth); params.SecondaryCatalogGrowth = 1; end
if (params.SecondaryCatalogGrowth <= 0) || ...
   isinf(params.SecondaryCatalogGrowth) || ...
   isnan(params.SecondaryCatalogGrowth)
    error('Invalid SecondaryCatalogGrowth parameter');
end

% RMM type and reduction factor
if ~isfield(params,'RemManeuver'); params.RemManeuver = []; end
if isempty(params.RemManeuver)
    params.RemManeuver = 'Translational';
end
if ~strcmpi(params.RemManeuver,'Translational') && ...
   ~strcmpi(params.RemManeuver,'Rotational')
    error('Invalid RemManeuver parameter');
end
if ~isfield(params,'RemReduction'); params.RemReduction = []; end
if isempty(params.RemReduction)
    if strcmpi(params.RemManeuver,'Translational')
        params.RemReduction = 0.03;
    else
        params.RemReduction = 0.20;
    end
end

% Set up output directory
if ~isfield(params,'runtag'); params.runtag = []; end
if isempty(params.runtag)
    params.runtag = [strrep(params.OCMDBroot,' ','_') ...
        '_' params.prisetstr];
    if ~strcmpi(params.redyel_event_option,'lud')
        params.runtag = [params.runtag '_Pc' params.redyel_event_option];
    end
end

params.timetag = make_timetag(params.start_time);

if ~isfield(params,'timestamp'); params.timestamp = []; end
if isempty(params.timestamp); params.timestamp = true; end

if params.timestamp
    params.outtag  = [params.runtag '_' params.timetag];
else
    params.outtag  = params.runtag;
end

if ~isfield(params,'outputpath'); params.outputpath = []; end
if isempty(params.outputpath); params.outputpath = 'output'; end

if ~isfield(params,'outputdir'); params.outputdir = []; end
if isempty(params.outputdir)
    params.outputdir = fullfile(params.outputpath,params.outtag);
end

if exist(params.outputdir,'dir')
    error('Output/log directory already exists.');
else
    mkdir(params.outputdir);
end

% Defaults for logging and displaying

if ~isfield(params,'logging'); params.logging = []; end
if isempty(params.logging); params.logging = true; end

if ~isfield(params,'displaying'); params.displaying = []; end
if isempty(params.displaying); params.displaying = true; end

% Set up log file

if params.logging
    params.log_file_name = fullfile(params.outputdir,[params.runtag '.log']);
    params.logfid = fopen(params.log_file_name,'wt');
else
    params.logfid = [];
end

% RIC screening volume type

if ~isfield(params,'SVRICtype') || isempty(params.SVRICtype)
    params.SVRICtype = 'box';
else
    params.SVRICtype = lower(params.SVRICtype);
end

switch params.SVRICtype
    case {'ellipsoid', 'ell'}
        params.SVRICtype_alt = 'RIC_ellipsoid';
    case 'box'
        params.SVRICtype_alt = 'RIC_box';
    otherwise
        error(['Invalid screening volume type: ' params.SVRICtype]);
end

% Array of RIC screening volumes for the primaries [1 x 3] or [Npri x 3]
if ~isfield(params,'SVRICdims_m') || isempty(params.SVRICdims_m)
    % params.SVRICdims_m = [0.5 17 20]*1000; % LEO 2-3 (A-train)
    params.SVRICdims_m = [2 44 51]*1000; % LEO 1-2 Grace-1 and -2
end

% Expand a [1 x 3] SVRICdims_m input array into an [Npri x 3] array, if
% required

szSVRICdims_m = size(params.SVRICdims_m);
if (szSVRICdims_m(1) == 1) && (szSVRICdims_m(2) == 3)
    params.SVRICdims_m = repmat(params.SVRICdims_m,[params.Npri,1]);
elseif (szSVRICdims_m(1) ~= params.Npri) || (szSVRICdims_m(2) ~= 3)
    error('Illegal dimension for SVRICdims_m parameter');
end

% Check if TCA tolerance for update association is provided. If not, use 15
% minutes
if ~isfield(params,'TCAtol_s') || isempty(params.TCAtol_s)
    params.TCAtol_s = 15*60;
end
if ~isfield(params,'TCRtol_s') || isempty(params.TCRtol_s)
    params.TCRtol_s = 15;
end

% Check if TCA limits are provided as a parameter

if ~isfield(params,'TCAlimits_UT'); params.TCAlimits_UT = []; end

% Check if a TCA bin width is provided as a parameter. If not, use 7 days.

if ~isfield(params,'TCAbin_days') || isempty(params.TCAbin_days)
    params.TCAbin_days = 7;
end

if ~isfield(params,'Pcred') || isempty(params.Pcred)
    params.Pcred = 1e-4;
end

if ~isfield(params,'Pcyel') || isempty(params.Pcyel)
    params.Pcyel = 1e-7;
end

if ~isfield(params,'Pcgre') || isempty(params.Pcgre)
    params.Pcgre = 1e-10;
end

params.log10_Pcred = log10(params.Pcred);
params.log10_Pcyel = log10(params.Pcyel);
params.log10_Pcgre = log10(params.Pcgre);

% Colors

params.acol = [0 0 0];
params.ycol = [255 (2/3)*255 0]/255;
params.rcol = [1 0 0];
params.gcol = [0 1 0];

% Random number generation mode either 'default' or 'shuffle'
if ~isfield(params,'rng_mode') || isempty(params.rng_mode)
    params.rng_mode = 'default';
end
rng(params.rng_mode); % Initialize random number generator

if strcmpi(params.rng_mode,'shuffle')
    warning('Random number mode set to shuffle - recommend default for reproducible results');
elseif ~strcmpi(params.rng_mode,'default')
    error('Illegal random number mode');
end

% Perform association checks

if ~isfield(params,'association_check') || isempty(params.association_check)
    params.association_check = false;
end

% Make event list CSV file

if ~isfield(params,'make_event_file') || isempty(params.make_event_file)
    params.make_event_file = false;
end
if ~isfield(params,'make_repeating_event_file') || isempty(params.make_repeating_event_file)
    params.make_repeating_event_file = false;
end
if ~isfield(params,'repeating_event_span_days') || isempty(params.repeating_event_span_days)
    params.repeating_event_span_days = 7;
end
if ~isfield(params,'repeating_event_Pc_threshold') || isempty(params.repeating_event_Pc_threshold)
    params.repeating_event_Pc_threshold = 1e-4;
end

% Make time bins and plots of rates, for all, yellow and red events

if ~isfield(params,'make_time_bins') || isempty(params.make_time_bins)
    params.make_time_bins = false;
end

if ~isfield(params,'make_time_plots') || isempty(params.make_time_plots)
    params.make_time_plots = false;
end

if numel(params.make_time_plots) == 1
    params.make_time_plots = repmat(params.make_time_plots,[3 1]);
end

if numel(params.make_time_plots) ~= 3
    error('Illegal dimension for for make_time_plots parameter');
end

if any(params.make_time_plots)
    
    params.make_time_bins = true;

    if ~isfield(params,'plot_err_bars') || isempty(params.plot_err_bars)
        params.plot_err_bars = [false true];
    end
    if numel(params.plot_err_bars) == 1
        params.plot_err_bars = repmat(params.plot_err_bars,[1 2]);
    end
    if numel(params.plot_err_bars) ~= 2
        error('Illegal dimension for for plot_err_bars parameter');
    end

    if ~isfield(params,'plot_conf_ranges') || isempty(params.plot_conf_ranges)
        params.plot_conf_ranges = [false false];
    end
    if numel(params.plot_conf_ranges) == 1
        params.plot_conf_ranges = repmat(params.plot_conf_ranges,[1 2]);
    end
    if numel(params.plot_conf_ranges) ~= 2
        error('Illegal dimension for for plot_conf_ranges parameter');
    end
    
end

if params.make_time_bins
    if ~isfield(params,'bin_samples') || isempty(params.bin_samples)
        params.bin_samples = 40; % Good for 95% ranges
    end
end

% Parameters for calculating accumulated risk and estimating Pc
% thresholds for remediation maneuvers

if ~isfield(params,'Pc_cutoff_accum_risk') || isempty(params.Pc_cutoff_accum_risk)
    params.Pc_cutoff_accum_risk = params.Pcgre;
end

if ~isfield(params,'N_bootstrap_accum_risk') || isempty(params.N_bootstrap_accum_risk)
    params.N_bootstrap_accum_risk = 1000;
end

if ~isfield(params,'lnPc_bootstrap_uncertainty') || isempty(params.lnPc_bootstrap_uncertainty)
    params.lnPc_bootstrap_uncertainty = 0.3;
end

if ~isfield(params,'add_one_unrem_event') || isempty(params.add_one_unrem_event)
    % TRUE  for modes 1 & 2 (conservative RMM thresholds)    
    params.add_one_unrem_event = true;
end

if ~isfield(params,'exclude_noncatastrophic') || isempty(params.exclude_noncatastrophic)
    % TRUE  for modes 2 & 4 (environmental protection)        
    params.exclude_noncatastrophic = false;
end

if ~isfield(params,'confidence_noncatastrophic') || isempty(params.confidence_noncatastrophic)
    params.confidence_noncatastrophic = 0.999;
end

if params.exclude_noncatastrophic
    if (params.confidence_noncatastrophic <= 0) || ...
       (params.confidence_noncatastrophic >= 1)
        error('Outside valid range: 0 < exclude_noncatastrophic < 1')
    else
        if (params.confidence_noncatastrophic < 0.95)
            warning('Recommend 0.95 <= exclude_noncatastrophic <= 0.999');
        end
    end
end

if ~isfield(params,'FredVsPremImage') || isempty(params.FredVsPremImage)
    params.FredVsPremImage = false;
end

% Parameters for plotting update sequences

if ~isfield(params,'make_update_plots') || isempty(params.make_update_plots)
    params.make_update_plots = false;
end

if ~isfield(params,'Pc_for_update_plots') || isempty(params.Pc_for_update_plots)
    params.Pc_for_update_plots = params.Pcred;
end

if ~isfield(params,'Nu_for_update_plots') || isempty(params.Nu_for_update_plots)
    params.Nu_for_update_plots = 0;
end

% Parameters for plotting radius-latitude distributions

if ~isfield(params,'make_altlat_plots') || isempty(params.make_altlat_plots)
    params.make_altlat_plots = false;
end

if numel(params.make_altlat_plots) == 1
    params.make_altlat_plots = repmat(params.make_altlat_plots,[3 1]);
end

if numel(params.make_altlat_plots) ~= 3
    error('Illegal dimension for make_altlat_plots parameter');
end

% Parameters for plots

if ~isfield(params,'visible_plots') || isempty(params.visible_plots)
    params.visible_plots = false;
end

if ~isfield(params,'plot_format') || isempty(params.plot_format)
    params.plot_format = '.png';
end
if ~strcmpi(params.plot_format(1),'.')
    params.plot_format = ['.' params.plot_format];
end

if ~isfield(params,'make_first_generation_plots') || isempty(params.make_first_generation_plots)
    params.make_first_generation_plots = false;
end

if ~isfield(params,'report_and_plot_svi_rates') || isempty(params.report_and_plot_svi_rates)
    params.report_and_plot_svi_rates = true;
end

return;
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