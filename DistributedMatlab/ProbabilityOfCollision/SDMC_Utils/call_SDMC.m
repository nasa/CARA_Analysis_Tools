function [SDMCPc,SDMCPcUnc,numHits,trialData] = call_SDMC(r1,v1,C1,r2,v2,C2,HBR,params,alwaysCompile)
% call_SDMC - Matlab wrapper to compile and call the mex code surrounding
%             the SDMC library.
%
% Syntax: [SDMCPc,SDMCPcUnc,numHits,trialData] = call_SDMC(r1,v1,C1,r2,v2,C2,HBR,params,alwaysCompile);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   This function uses a Monte Carlo method to calculate the probability of
%   collision (Pc) between two space objects for a single conjunction,
%   given input states and covariances at the nominal TCA.
%
%   IMPORTANT NOTE:
%     In order to run call_SDMC, a system environment variable must be set
%     using the absolute path to the following path:
%
%       SDK\DistributedMatlab\ProbabilityOfCollision\SDMC_Utils\lib
%
%     The variable used is different per operating system.
%
%       Linux uses: LD_LIBRARY_PATH
%       Windows uses: PATH
%
%     If the system environment variable isn't set correctly, the
%     call_SDMC.m routine will display an error message dictating the value
%     which should be stored within the appropriate system environment
%     variable.
%
%     This variable must be set outside of Matlab (e.g. via the "Edit
%     system environment variables" dialog in Windows or in the .bashrc or
%     .profile file on Linux) and Matlab must be restarted once the
%     variable is set correctly.
%
%     This has been updated so that on Linux for Matlab ver 9.10+ setting
%     the LD_LIBRARY_PATH is unnecessary. Also, note that if recompliation
%     of the mex files is required, then you must have the appropriate
%     Intel Fortran libraries in LD_LIBRARY_PATH.
%
% =========================================================================
%
% INPUT:
%
%    r1      - Primary object's ECI position vector (m)        [3x1 or 1x3]
%    v1      - Primary object's ECI velocity vector (m/s)      [3x1 or 1x3]
%    C1      - Primary object's ECI covariance matrix                 [6x6]
%    r2      - Secondary object's ECI position vector (m)      [3x1 or 1x3]
%    v2      - Secondary object's ECI velocity vector (m/s)    [3x1 or 1x3]
%    C2      - Secondary object's ECI covariance matrix               [6x6]  
%    HBR     - Combined primary+secondary hard-body radii (m)         [1x1]
%
%    params  - Auxilliary input parameter structure, described in detail
%              in function "default_params_Pc_SDMC".
%
%    alwaysCompile - (Optional) Boolean which indicates if the mex file
%                    should always be recompiled.
%                    Defaults to false
%
% =========================================================================
%
% OUTPUT:
%
%   SDMCPc    - The estimated conjunction Pc value.                   [1x1]
%   SDMCPcUnc - The uncertainty in the estimated Pc value based on    [1x2]
%               params.conf_level.
%   numHits   - The number of trials resulting in a "hit".            [1x1]
%   trialData - Detailed trial data information up to               [table]
%               params.max_output_trials number of rows. Includes
%               the following data:
%      hitIndicator - Indicates if the trial CA distance < HBR
%      hitTime - Time of impingmement (for hit) or global PCA (for miss)
%      hitTimeMiss - Time of hit miss distance or global PCA
%      hitRadius - Hit miss distance or global Pca distance (m)
%      pcaTime - Time of global PCA
%      pcaRadius - global PCA distance (m)
%      priPosX_km - primary inertial x position at hitTimeMiss (km)
%      priPosY_km - primary inertial y position at hitTimeMiss (km)
%      priPosZ_km - primary inertial z position at hitTimeMiss (km)
%      priVelX_kmps - primary inertial x velocity at hitTimeMiss (km/sec)
%      priVelY_kmps - primary inertial y velocity at hitTimeMiss (km/sec)
%      priVelZ_kmps - primary inertial z velocity at hitTimeMiss (km/sec)
%      secPosX_km - secondary inertial x position at hitTimeMiss (km)
%      secPosY_km - secondary inertial y position at hitTimeMiss (km)
%      secPosZ_km - secondary inertial z position at hitTimeMiss (km)
%      secVelX_kmps - secondary inertial x velocity at hitTimeMiss (km/sec)
%      secVelY_kmps - secondary inertial y velocity at hitTimeMiss (km/sec)
%      secVelZ_kmps - secondary inertial z velocity at hitTimeMiss (km/sec)
%
%      Note: All times are in days since the Jan 0, 1970 epoch.
%
% =========================================================================
%
% Initial version: Mar 2023; Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    %% Check input parameters
    if nargin == 8
        alwaysCompile = false;
    elseif nargin ~= 9
        error('Incorrect number of parameters passed in');
    end
    
    %% Setup proper pathing and make sure the mex interface has been compiled
    % Add paths
    mpath = fileparts(mfilename('fullpath'));
    persistent pathsAdded
    if isempty(pathsAdded)
        s = what(mpath); addpath(s.path);
        s = what(fullfile(mpath, 'lib')); addpath(s.path);
        pathsAdded = true;
    end
    
    % Make sure compile/runtime paths are properly set
    includepath = fullfile(mpath, 'include');
    libpath = fullfile(mpath, 'lib');
    ldpath = getenv('LD_LIBRARY_PATH');
    if ~endsWith(ldpath, [':' libpath])
        setenv('LD_LIBRARY_PATH',[ldpath ':' libpath]);
    end
    
    % We change directories to the path where the script resides so that
    % the mex file is created here and not where the script was run from.
    cwd = cd(mpath);
    % Always make sure to change directories back to the original directory
    % whenever this script exits
    cleanup_obj = onCleanup(@() cd(cwd));

    if isunix
        if verLessThan('matlab','9.4')
            error('Matlab must be version R2018a or greater to use SDMC on Linux');
        elseif verLessThan('matlab','9.10')
            matlab_sdmc_wrapper = 'matlab_sdmc_wrapper_R2018a';
        else
            matlab_sdmc_wrapper = 'matlab_sdmc_wrapper_R2021a';
        end
    elseif ispc
        if verLessThan('matlab','9.3')
            error('Matlab must be version R2017b or greater to use SDMC on Windows');
        else
            matlab_sdmc_wrapper = 'matlab_sdmc_wrapper.mexw64';
        end
    else
        error('SDMC does not support this OS');
    end
    
    % Check if the mex file needs to be compiled
    if ~isequal(fileparts(which(matlab_sdmc_wrapper)), mpath) || alwaysCompile
        % Compile the mex
        includeopt = ['-I' includepath];
        libopt = ['-L' libpath];
        if isunix
            if verLessThan('matlab','9.10')
                mex('-R2018a',includeopt,libopt,...
                    '-lsdmc','-lastroblkd','-lastrosupt','-lastrosuptts',...
                    '-ltimesupt','-lifcore','-lifport','-limf','-lintlc',...
                    'matlab_sdmc_wrapper.cpp','sdmcEntry.cpp','-output',matlab_sdmc_wrapper);
            else
                mex('-R2018a',includeopt,libopt,...
                    '-lsdmc_linux',...
                    '-lifcore','-lifport','-limf','-lintlc','-lsvml',...
                    'matlab_sdmc_wrapper.cpp','sdmcEntry.cpp', ...
                    'LDFLAGS=$LDFLAGS -Wl,-rpath=''$ORIGIN''./lib', ... % linker flag - avoid needing environment var
                    '-output', matlab_sdmc_wrapper);
            end
        elseif ispc
            mex('-R2018a',includeopt,libopt,...
                '-lsdmc_win','matlab_sdmc_wrapper.cpp','sdmcEntry.cpp');
        end
    end
    
    %% Setup data structure for calling the SDMC library
    % This variable allows the conversion to the SDMC epoch (12/31/1969)
    epoch1970 = datenum(1970, 0, 0, 0, 0, 0);
    
    % Primary satellite
    stateInfo.pri_satno = params.pri_objectid;
    stateInfo.pri_epoch = datenum(params.pri_epoch) - epoch1970; % days since 1970 epoch
    stateInfo.pri_pos_x = r1(1)./1e3;   % km
    stateInfo.pri_pos_y = r1(2)./1e3;   % km
    stateInfo.pri_pos_z = r1(3)./1e3;   % km
    stateInfo.pri_vel_x = v1(1)./1e3;   % km/sec
    stateInfo.pri_vel_y = v1(2)./1e3;   % km/sec
    stateInfo.pri_vel_z = v1(3)./1e3;   % km/sec
    stateInfo.pri_cov_11 = C1(1,1);     % m^2
    stateInfo.pri_cov_21 = C1(2,1);     % m^2
    stateInfo.pri_cov_22 = C1(2,2);     % m^2
    stateInfo.pri_cov_31 = C1(3,1);     % m^2
    stateInfo.pri_cov_32 = C1(3,2);     % m^2
    stateInfo.pri_cov_33 = C1(3,3);     % m^2
    stateInfo.pri_cov_41 = C1(4,1);     % m^2/s
    stateInfo.pri_cov_42 = C1(4,2);     % m^2/s
    stateInfo.pri_cov_43 = C1(4,3);     % m^2/s
    stateInfo.pri_cov_44 = C1(4,4);     % m^2/s^2
    stateInfo.pri_cov_51 = C1(5,1);     % m^2/s
    stateInfo.pri_cov_52 = C1(5,2);     % m^2/s
    stateInfo.pri_cov_53 = C1(5,3);     % m^2/s
    stateInfo.pri_cov_54 = C1(5,4);     % m^2/s^2
    stateInfo.pri_cov_55 = C1(5,5);     % m^2/s^2
    stateInfo.pri_cov_61 = C1(6,1);     % m^2/s
    stateInfo.pri_cov_62 = C1(6,2);     % m^2/s
    stateInfo.pri_cov_63 = C1(6,3);     % m^2/s
    stateInfo.pri_cov_64 = C1(6,4);     % m^2/s^2
    stateInfo.pri_cov_65 = C1(6,5);     % m^2/s^2
    stateInfo.pri_cov_66 = C1(6,6);     % m^2/s^2
    stateInfo.pri_pos_sensitivity_x = params.covXcorr.Gvecp(1); % m
    stateInfo.pri_pos_sensitivity_y = params.covXcorr.Gvecp(2); % m
    stateInfo.pri_pos_sensitivity_z = params.covXcorr.Gvecp(3); % m
    stateInfo.pri_vel_sensitivity_x = params.covXcorr.Gvecp(4); % m/s
    stateInfo.pri_vel_sensitivity_y = params.covXcorr.Gvecp(5); % m/s
    stateInfo.pri_vel_sensitivity_z = params.covXcorr.Gvecp(6); % m/s
    stateInfo.pri_forecast_uncertainty = params.covXcorr.sigp;  % unitless
    
    % Secondary satellite
    stateInfo.sec_satno = params.sec_objectid;
    stateInfo.sec_epoch = datenum(params.sec_epoch) - epoch1970; % days since 1970 epoch
    stateInfo.sec_pos_x = r2(1)./1e3;   % km
    stateInfo.sec_pos_y = r2(2)./1e3;   % km
    stateInfo.sec_pos_z = r2(3)./1e3;   % km
    stateInfo.sec_vel_x = v2(1)./1e3;   % km/sec
    stateInfo.sec_vel_y = v2(2)./1e3;   % km/sec
    stateInfo.sec_vel_z = v2(3)./1e3;   % km/sec
    stateInfo.sec_cov_11 = C2(1,1);     % m^2
    stateInfo.sec_cov_21 = C2(2,1);     % m^2
    stateInfo.sec_cov_22 = C2(2,2);     % m^2
    stateInfo.sec_cov_31 = C2(3,1);     % m^2
    stateInfo.sec_cov_32 = C2(3,2);     % m^2
    stateInfo.sec_cov_33 = C2(3,3);     % m^2
    stateInfo.sec_cov_41 = C2(4,1);     % m^2/s
    stateInfo.sec_cov_42 = C2(4,2);     % m^2/s
    stateInfo.sec_cov_43 = C2(4,3);     % m^2/s
    stateInfo.sec_cov_44 = C2(4,4);     % m^2/s^2
    stateInfo.sec_cov_51 = C2(5,1);     % m^2/s
    stateInfo.sec_cov_52 = C2(5,2);     % m^2/s
    stateInfo.sec_cov_53 = C2(5,3);     % m^2/s
    stateInfo.sec_cov_54 = C2(5,4);     % m^2/s^2
    stateInfo.sec_cov_55 = C2(5,5);     % m^2/s^2
    stateInfo.sec_cov_61 = C2(6,1);     % m^2/s
    stateInfo.sec_cov_62 = C2(6,2);     % m^2/s
    stateInfo.sec_cov_63 = C2(6,3);     % m^2/s
    stateInfo.sec_cov_64 = C2(6,4);     % m^2/s^2
    stateInfo.sec_cov_65 = C2(6,5);     % m^2/s^2
    stateInfo.sec_cov_66 = C2(6,6);     % m^2/s^2
    stateInfo.sec_pos_sensitivity_x = params.covXcorr.Gvecs(1); % m
    stateInfo.sec_pos_sensitivity_y = params.covXcorr.Gvecs(2); % m
    stateInfo.sec_pos_sensitivity_z = params.covXcorr.Gvecs(3); % m
    stateInfo.sec_vel_sensitivity_x = params.covXcorr.Gvecs(4); % m/s
    stateInfo.sec_vel_sensitivity_y = params.covXcorr.Gvecs(5); % m/s
    stateInfo.sec_vel_sensitivity_z = params.covXcorr.Gvecs(6); % m/s
    stateInfo.sec_forecast_uncertainty = params.covXcorr.sigs;  % unitless
    
    % General conjuction inputs
    stateInfo.tca                   = datenum(params.TCA) - epoch1970;     % days since 1970 epoch
    stateInfo.trajectory_mode       = params.trajectory_mode;              % 0 = 2-body, 1 = rectilinear, 2 = rectiliniear (position deviations only)
    stateInfo.span_days             = params.span_days;                    % days
    stateInfo.seed                  = params.seed;
    stateInfo.num_trials            = params.num_trials;
    stateInfo.hbr_m                 = HBR;                                 % m
    stateInfo.max_radius            = params.max_radius;                   % m
    stateInfo.max_num_output_trials = params.max_output_trials;
    stateInfo.sdmc_list_file        = params.sdmc_list_file;               % ' ' = None, 'STDOUT' = print to STDOUT, otherwise print to file name

    %% Run SDMC and retrieve outputs
    try
        try
            % str2fun allows us to use the full path, so the mex need not be on our path
            obj = mex_interface(str2fun([mpath '/' matlab_sdmc_wrapper]));
        catch ME
            if isunix
                envVarName = 'LD_LIBRARY_PATH';
                envVarType = 'Linux';
            elseif ispc
                envVarName = 'PATH';
                envVarType = 'Windows';
            else
                rethrow(ME);
            end
            error(['Could not run mex file, library path is likely set incorrectly.' newline ...
                'Set ' envVarName ' environment variable to include: ' libpath newline ...
                'Restart Matlab after setting the ' envVarType ' System environment variable.']);
        end

        % Initialize the state of the object based on data passed in
        retVal = obj.run_sdmc(stateInfo);
        
        % Check for invalid 6x6 covariance and a warning or error state was
        % returned
        if retVal ~= 0 && (isequal(zeros(3,3),C1(1:3,4:6)) || isequal(zeros(3,3),C2(1:3,4:6)))
            msgId = 'SDMC:invalid_6x6_cov';
            msgTxt = ['Invalid primary or secondary 6x6 covariance caused SDMC runtime warning or error: ' num2str(retVal)];
            error(msgId,msgTxt);
        end
        
        if retVal > 0
            % Greater than 0 are warnings
            if retVal == 1
                msgId  = 'SDMC:hyperbolic_primary';
                msgTxt = 'Hyperbolic primary epoch elements';
            elseif retVal == 2
                msgId  = 'SDMC:rectilinear_primary';
                msgTxt = 'Rectilinear primary epoch elements';
            elseif retVal == 3
                msgId  = 'SDMC:retrograde_primary';
                msgTxt = 'Retrograde equatorial primary epoch elements';
            elseif retVal == 4
                msgId  = 'SMDC:hyperbolic_secondary';
                msgTxt = 'Hyperbolic secondary epoch elements';
            elseif retVal == 5
                msgId  = 'SDMC:rectilinear_secondary';
                msgTxt = 'Rectilinear secondary epoch elements';
            elseif retVal == 6
                msgId  = 'SDMC:retrograde_secondary';
                msgTxt = 'Retrograde equatorial secondary epoch elements';
            elseif retVal == 7
                msgId  = 'SDMC:max_ephem_points';
                msgTxt = 'Max internal ephemeris points exceeded (should not occur)';
            elseif retVal == 8
                msgId  = 'SDMC:too_few_ephem_points';
                msgTxt = 'Too few internal ephemeris points (5)  (should not occur)';
            elseif retVal == 9
                msgId  = 'SDMC:too_many_bad_elems';
                msgTxt = 'Too many bad elements during MC trials';
            elseif retVal == 10
                msgId  = 'SDMC:cannot_open_list_file';
                msgTxt = 'Could not open output list file';
            else
                msgId  = 'SDMC:unknown';
                msgTxt = 'Unknown warning';
            end
            warning(msgId,msgTxt);
        elseif retVal < 0
            % Less than 0 are errors
            if retVal == -1
                msgId  = 'SDMC:primary_npd';
                msgTxt = 'NPD encountered for primary with remediation';
            elseif retVal == -2
                msgId  = 'SDMC:secondary_npd';
                msgTxt = 'NPD encountered for secondary with remediation';
            elseif retVal == -3
                msgId  = 'SDMC:small_rel_vel';
                msgTxt = 'Small TCA relative velocity condition';
            elseif retVal == -4
                msgId  = 'SDMC:small_cov';
                msgTxt = 'Small epoch covariance condition';
            elseif retVal == -5
                msgId  = 'SDMC:badElems';
                msgTxt = 'Bad elements encountered during MC trials';
            else
                msgId  = 'SDMC:unknown';
                msgTxt = 'Unknown error';
            end
            error(msgId,msgTxt);
        end

        % Get the SDMC outputs
        [numHits, trialData] = obj.get_sdmc_outputs();

        % Clears memory
        clear obj;
    catch ME
        % Always clear the memory if the object exists
        if exist('obj','var')
            clear obj;
        end
        rethrow(ME);
    end
    
    % Binofit statistics
    conf_alpha = 1-params.conf_level;
    [SDMCPc,SDMCPcUnc] = binofit(numHits,params.num_trials,conf_alpha);
    
    % Fill in the output table
    if isempty(trialData)
        trialData = [];
    else
        trialData = table(...
            trialData(:,1),trialData(:,2),trialData(:,3),trialData(:,4),...
            trialData(:,5),trialData(:,6),trialData(:,7),trialData(:,8),trialData(:,9),...
            trialData(:,10),trialData(:,11),trialData(:,12),trialData(:,13),trialData(:,14),...
            trialData(:,15),trialData(:,16),trialData(:,17),trialData(:,18),...
            'VariableNames',...
            {'hitIndicator','hitTime','hitTimeMiss','hitRadius',...
             'pcaTime','pcaRadius',...
             'priPosX_km',  'priPosY_km',  'priPosZ_km',...
             'priVelX_kmps','priVelY_kmps','priVelZ_kmps',...
             'secPosX_km',  'secPosY_km',  'secPosZ_km',...
             'secVelX_kmps','secVelY_kmps','secVelZ_kmps'});
    end
end


% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-13 | Initial Development
% L. Baars       | 2024-May-28 | Updated error text when system environment
%                                variable is set incorrectly, added note to
%                                the description about these variables.
% L. Baars       | 2025-Aug-06 | Minor documentation updates in preparation
%                                for public release.
% S. Shaw        | 2025-Nov-13 | Added linker flag so that the
%                                LD_LIBRARY_PATH is not necessary on linux
%                                for Matlab version 9.10+.
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================