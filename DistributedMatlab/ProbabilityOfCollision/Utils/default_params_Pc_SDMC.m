function params = default_params_Pc_SDMC(params)
% default_params_Pc_SDMC - Add and/or set the defaults for the
%                          parameters used by the function Pc_SDMC.
%
% Syntax: params = default_params_Pc_SDMC(params);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% INPUT:
%
% params = (Optional) Empty or partially populated structure containing the
%          parameters used by the Pc_SDMC function, containing the
%          following fields:
%
%   Target95pctPcAccuracy = Defines the target accuracy for the 
%                           SDMC method Pc estimation (recommended = 0.1)
%
%   max_num_trials = Defines the maximum number of trials allowed for any
%                    SDMC method Pc estimation (recommended = 3.7e7)
%
%   default_num_trials = Defines the default number of trials to use when
%                        an analytical Pc cannot be calculated
%                        (recommended = 3.7e7)
%
%   trajectory_mode = Defines the type of motion to use when running the
%                     Monte Carlo trials
%                     0 => (default) 2-body motion
%                     1 => full rectilinear motion
%                     2 => rectilinear motion (position derivatives only)
%
%   sdmc_list_file = Defines where the outputs from the SDMC library are
%                    printed
%                    ''       => (default) no outputs are printed
%                    'STDOUT' => prints to STDOUT
%                    fileName => prints to file defined by fileName
%
%   use_parallel = Boolean defining whether or not the SDMC library should
%                  be run in parallel. Only available with Matlab's
%                  Parallel Computing Toolbox.
%                  Defaults to true
%
%   seed = Negative integer defining the seed for the SDMC Monte Carlo
%          random number generator. If use_parallel is true, the abs(seed)
%          is used to seed Matlab's random number generator, which is then
%          used to create random seeds for each parallel process. Valid
%          values are -1 to -2147483648.
%          Defaults to -1
%
%   conf_level = Confidence level for Pc uncertainty range representing a
%                percentage of certainty. Valid values are floating point
%                numbers between 0 and 1 (not inclusive).
%                Defaults to 0.95 (i.e. 95% certainty)
%
%   generate_ca_dist_plot = Boolean defining if the SDMC library should
%                           generate the data needed to plot close approach
%                           distributions.
%                           Defaults to false
%
%   num_trials = The number of Monte Carlo trials to run. If no value is
%                specified Pc_SDMC will calculate an estimated number of
%                trials using Pc3D_Hall in order to receive approximately
%                400 hits (i.e. 10% accuracy).
%                Defaults to []
%
%   num_workers = The number of separate parallel processes to use when
%                 running SDMC. If this number exceeds the number of
%                 processors available in the parallel pool, then SDMC will
%                 be split into num_workers jobs which will then be run
%                 within the parallel pool. If the value passed in is
%                 empty, then Pc_SDMC will use the number of processors in
%                 the parallel pool as num_workers. This option can be used
%                 to replicate results across machines with different
%                 number of processors. This option cannot be used without
%                 the parallel toolbox.
%                 Defaults to []
%
%   span_days = The span (in days) to propagate the orbits around the time
%               of closest approach (TCA). If no value is specified Pc_SDMC
%               will calculate an optimal span based on the Nc3D
%               conjunction limits (from Pc3D_Hall).
%               Defaults to []
%
%   Nc3D = The "Pc" output parameter from a call to Pc3D_Hall. Pc_SDMC will
%          call Pc_SDMC() itself if this parameter is not specified.
%          Defaults to []
%
%   Nc3DInfo = The "out" output parameter from a call to Pc3D_Hall. Pc_SDMC
%              will call Pc_SDMC() itself if this parameter is not
%              specified.
%              Defaults to []
%
%   apply_covXcorr_corrections = Determines if covariance cross-correlation
%                                effects should be applied to the
%                                conjunction. Depends on a properly filled
%                                out params.covXcorr sub-structure. See
%                                get_covXcorr_parameters documentation for
%                                a definition of that structure.
%                                Defaults to true
%
%   pri_objectid = Primary object ID
%                  Defaults to 1
%
%   sec_objectid = Secondary object ID
%                  Defaults to 2
%
%   TCA = Date/time field defining the time of closest approach for the
%         conjunction passed to Pc_SDMC
%         Defaults to 1/1/1990 00:00:00
%
%   pri_epoch = Date/time field defining the epoch of the primary object
%               state passed to Pc_SDMC. SDMC will propagate the state to
%               TCA before processing the conjunction.
%               Defaults to the TCA value
%
%   sec_epoch = Date/time field defining the epoch of the secondary object
%               state passed to Pc_SDMC. SDMC will propagate the state to
%               TCA before processing the conjunction.
%               Defaults to the TCA value
%
%   verbose = Verbosity (0 - none, 1 - some, 2 - more).
%
% =========================================================================
%
% OUTPUT:
%
%   params = Fully-populated structure containing parameters used by the 
%            function Pc_SDMC.
%
% =========================================================================
%
% Initial version: Feb 2023; Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------

    % Initializions and defaults
    Nargin = nargin;
    if Nargin < 1; params = []; end
    
    % Default trajectory mode
    %   0 = 2-body motion
    %   1 = full rectilinear motion
    %   2 = rectilinear motion (position deviations only)
    params = set_default_param(params,'trajectory_mode',0);
    if params.trajectory_mode ~= 0 && params.trajectory_mode ~= 1 && ...
            params.trajectory_mode ~= 2
        error('Supplied trajectory_mode must be 0, 1, or 2');
    end

    % Target 95% confidence Pc estimation accuracy
    % (0.1 corresponds to ~400 hits; 0.2 corresponds to ~100 hits)
    params = set_default_param(params,'Target95pctPcAccuracy',0.1);
    
    % Max number of trials allowed (to prevent needlessly long runs)
    params = set_default_param(params,'max_num_trials',3.7e7);
    
    % Default number of trials (when a Pc cannot be calculated
    % analytically). The value of 3.7e7 verifies that the 95% confidence
    % limit at 10% accuracy will be just below a Pc of 1e-7 if 0 hits are
    % encountered.
    params = set_default_param(params,'default_num_trials',3.7e7);
    
    % Outputs from SDMC library
    %  ''       = no outputs
    %  'STDOUT' = printed to standard out
    %  fileName = printed to a file name (fileName cannot exceed 255
    %             characters)
    params = set_default_param(params,'sdmc_list_file','');
    if ~ischar(params.sdmc_list_file) || length(params.sdmc_list_file) > 255
        error('Supplied sdmc_list_file must be a character vector less than or equal to 255 characters');
    end
    
    % Process SDMC with parallel processes
    params = set_default_param(params,'use_parallel',true);
    % Check if the parallel processing toolbox is available, and set the
    % use_parallel parameter to false if it isn't
    if ~license('test','distrib_computing_toolbox')
        params.use_parallel = false;
    end
    
    % Default seed
    params = set_default_param(params,'seed',-1);
    if params.seed >= 0 || isinf(params.seed) || isnan(params.seed) || ...
            floor(params.seed) ~= params.seed || params.seed < -2147483648
        error('Supplied seed must be a negative integer greater than -2,147,483,649');
    end
    
    % Default confidence level for Pc uncertainty
    params = set_default_param(params,'conf_level',0.95);
    if params.conf_level <= 0 || params.conf_level >= 1
        error('Supplied conf_level must be a number between 0 and 1 (not inclusive)');
    end
    
    % Conjunction close approach distribution control
    params = set_default_param(params,'generate_ca_dist_plot',false);
    
    % Check the number of trials parameter
    if ~isfield(params,'num_trials'); params.num_trials = []; end
    if ~isempty(params.num_trials)
        if isinf(params.num_trials) || isnan(params.num_trials) || ...
                params.num_trials <= 0 || floor(params.num_trials) ~= params.num_trials
            error('Supplied num_trials must be a positive integer value');
        elseif params.num_trials > 1e9
            warning(['Supplied num_trials (' num2str(params.num_trials) ') is above the recommended maximum of 1e9' newline ...
                'It is expected that this process will run for a long time']);
        end
    end
    % If the number of trials aren't passed in, then the necessary number
    % of trials will be calculated by the Pc_SDMC process
    
    % Set the number of workers
    if ~isfield(params,'num_workers'); params.num_workers = []; end
    if ~isempty(params.num_workers)
        if isinf(params.num_workers) || isnan(params.num_workers) || ...
                params.num_workers <= 0 || floor(params.num_workers) ~= params.num_workers
            error('Supplied num_workers must be a positive integer value');
        end
    end
    % If the number of workers aren't passed in, the Pc_SDMC will use the
    % parallel pool default size as the number of workers
    
    % Check the span days parameter
    if ~isfield(params,'span_days'); params.span_days = []; end
    if ~isfield(params,'tmid'); params.tmid = []; end
    % It is assumed that Pc_SDMC will fill in these parameters if they are
    % empty
    
    % InputPc is used to determine the number of trials to run
    if ~isfield(params,'InputPc'); params.InputPc = []; end

    % Initialize covariance cross correlation correction indicator flag
    params = set_default_param(params,'apply_covXcorr_corrections',true);
    
    % Check primary and secondary object IDs
    params = set_default_param(params,'pri_objectid',1);
    params = set_default_param(params,'sec_objectid',2);
    
    % Check TCA, pri epoch, and sec epoch. These values should all be the
    % same. If they are not, then provide a warning. Set a default value of
    % Jan 1, 1990 if no data is provided. It doesn't really matter what
    % date is used since SDMC only uses 2-body or rectilinear motion.
    params = set_default_param(params,'TCA',datetime(1990,1,1));
    params = set_default_param(params,'pri_epoch',params.TCA);
    params = set_default_param(params,'sec_epoch',params.TCA);
    if params.TCA ~= params.pri_epoch || params.TCA ~= params.sec_epoch
        warning(['Supplied TCA, pri_epoch, and sec_epoch do not match.' newline ...
            'This is not a recommended run configuration for SDMC.']);
    end
    
    % Set the default retrograde orbit reorientation mode
    %  0 => No retrograde orbit adjustment (causes 2D-Nc and 3D-Nc to fail for
    %       retrograde orbits, such as the Alfano 2009 test cases)
    %  1 => If either orbits is retrograde, try reorienting the reference 
    %       frame axes (recommended)
    %  2 => Always try reorienting the ref. frame axes (testing mode)
    %  3 => Reorient axes to force primary to be retrograde (testing mode)
    params = set_default_param(params,'RetrogradeReorientation',1);
    
    % Validate target accuracy parameter
    if numel(params.Target95pctPcAccuracy) == 1
        % Process single target accuracy
        if params.Target95pctPcAccuracy <= 0 || ...
           params.Target95pctPcAccuracy >= 1
            warning('Invalid Target95pctPcAccuracy parameter; setting to 10%');
            params.Target95pctPcAccuracy = 0.1;
        end
    else
        % Procsess target accuracy Nx2 table holding log10(Pc) and TgtAcc
        % values
        tbl = params.Target95pctPcAccuracy;
        siz = size(tbl);
        if siz(1) <= 1 || siz(2) ~= 2
            warning('Invalid Target95pctPcAccuracy table; setting to 10%');
            params.Target95pctPcAccuracy = 0.1;
        end
    end
    
    % Default warning level
    params = set_default_param(params,'warning_level',1);    

    % Default verbosity.
    params = set_default_param(params,'verbose',false);

end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-MAR-10 | Initial Development
% D. Hall        | 2023-07-27  | Added retrograde orbit processing
% D. Hall        | 2023-08-18  | Added SDMC-Pc target accuracy processing
% L. Baars       | 2023-08-31  | Added the default number of trials param
% L. Baars       | 2025-08-06  | Minor documentation updates necessary for
%                                public release
% L. Baars       | 2025-09-03  | Added num_workers parameter

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================