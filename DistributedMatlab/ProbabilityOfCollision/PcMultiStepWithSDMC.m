function [Pc,out] = PcMultiStepWithSDMC(r1,v1,C1,r2,v2,C2,HBR,params)
% PcMultiStepWithSDMC - Calculate the collision probability for a
%                       conjunction using a multi-tiered algorithm,
%                       including the Simple Dynamics Monte Carlo (SDMC)
%                       calculation.
%
% Syntax: [Pc, out] = PcMultiStepWithSDMC(r1,v1,C1,r2,v2,C2,HBR);
%         [Pc, out] = PcMultiStepWithSDMC(r1,v1,C1,r2,v2,C2,HBR,params);
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
% This function builds on PcMultiStep by adding the option of executing
% Pc_SDMC after running PcMultiStep, in order to estimate the collision
% probability for a single conjunction using a 2-Body Monte Carlo
% algorithm. By default, the SDMC algorithm will only be run if there are
% usage violations encountered in all PcMultiStep algorithms or if the Pc
% returned by PcMultiStep is greater than or equal to the Pc defined in the
% PcToForceSDMCCalculation parameter. Alternatively, a user can force the
% SDMC Pc to be calculated uing the ForceSDMCCalculation parameter.
%
% Note: "BFMC" is a tool used internal to NASA CARA and will not be
%       publicly released. References to BFMC within the code should be
%       ignored by external users. Minimal documentation will be provided
%       for any BFMC functionality which may appear within the code.
%
% =========================================================================
%
% Inputs:
%
%    r1 - Primary object's ECI position vector (m)             [3x1 or 1x3]
%
%    v1 - Primary object's ECI velocity vector (m/s)           [3x1 or 1x3]
%
%    C1 - Primary object's ECI covariance matrix (m position units)   [6x6]
%
%    r2 - Secondary object's ECI position vector (m)           [3x1 or 1x3]
%
%    v2 - Secondary object's ECI velocity vector (m/s)       [3x1] or [1x3]
%
%    C2 - Secondary object's ECI covariance matrix (m position units) [6x6]  
%
%    HBR - Combined primary+secondary hard-body radii.   [1x1, 1x2, or 2x1]
%          If two elemens are passed in, it is assumed
%          each element separately represents each
%          object's individual HBR. In this case, the HBRs
%          are added together to create a single unified
%          combined HBR.
%
%    params  - (Optional) Other input parameters structure which expands on
%              the structure defined by PcMultiStep.m. The following fields
%              are specific to this function:
%
%      ForceSDMCCalculation - Boolean which forces the SDMC calculation to
%                             be performed, regardless of usage violations
%                             encountered.
%                             Defaults to false.
%
%      PreventSDMCCalculation - Boolean which prevents the SDMC calculation
%                               from running. An SDMC calculation cannot
%                               both be forced and prevented.
%                               Defaults to false.
%
%      SDMCParams - The "params" structure passed to Pc_SDMC. See Pc_SDMC.m
%                   for a definition of allowed fields within this
%                   structure.
%                   Defaults to an empty structure.
%
%      PcToForceSDMCCalculation - A threshold (between 0 and 1) which
%                                 defines a Pc returned from PcMultiStep
%                                 above which an SDMC Pc will be
%                                 calculated, regardless of usage
%                                 violations.
%                                 Defaults to 1.
%
%      ForceRecommendedTrials - Boolean which forces SDMC to use the
%                               recommended number of trials regardless of
%                               the maximum trial constraint. Note: This
%                               could lead to very long SDMC run times.
%                               Defaults to false.
%
% =========================================================================
%
% Outputs:
%
%    Pc - Recommended Pc value to use for this conjunction
%
%    out - Auxiliary output information structure which adds to the "out"
%          structure provided by PcMultiStep.m. This function adds or
%          modifies the following output fields:
%
%      SDMCPc - Pc calculated from the SDMCPc method (Pc_SDMC.m)
%
%      PcMethod - Updates PcMethod from PcMultiStep.m, description of
%                 method used to calculate the recommended Pc
%
%      PcMethodNum - Updates PcMethodNum from PcMultiStep.m, adds the
%                    following method numbers:
%        4 = SDMC Pc
%        5 = BFMC Pc
%
%      PcMethodNumAnalytical - Value of PcMethodNum that was returned by
%                              PcMultiStep.m.
%
%      PcMethodMax - Updates PcMethodMax from PcMultiStep.m, adds the
%                    following method numbers:
%        4 = SDMC Pc
%        5 = BFMC Pc
%
%      PcMethodMaxAnalytical - Value of PcMethodMax that was returned by
%                              PcMultiStep.m.
%
%      SDMCInfo - Copy of the "out" structure returned by Pc_SDMC.m. The
%                 following fields are added by this function:
%
%        PcCalculated - Boolean indicating that an SDMC Pc was calculated.
%
%        Indicators - Structure containing values calculated for SDMC
%                     usage violation checks, fields include:
%
%          Extended - Extended conjunction duration indicator (# of "hits"
%                     that are outside encounter start/stop time
%                     boundaries)
%
%        Violations - Structure of booleans showing SDMC usage violations,
%                     fields include:
%
%          Extended - True if extended indicator is greater than 0
%
%          Offset - True if any Nc-3D offset violations exist.
%
%      AnySDMCViolations - Boolean which indicates if any Exteded or Offset
%                          SDMC usage violations exist.
%
%      SDMCParams - A copy of the "params" structure that was used when
%                   calling Pc_SDMC.m.
%
% =========================================================================
%
% Disclaimer:
%
%    No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
%    WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
%    INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE
%    WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
%    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
%    INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
%    FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
%    THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
%    CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT
%    OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY
%    OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.
%    FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES
%    REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE,
%    AND DISTRIBUTES IT "AS IS."
%
%    Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
%    AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
%    SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF
%    THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES,
%    EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM
%    PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT
%    SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED
%    STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
%    PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE
%    REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL
%    TERMINATION OF THIS AGREEMENT.
%
% =========================================================================
%
% Initial version: Aug 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % Add required library paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        s = what(fullfile(p, '../Utils/PosVelTransformations')); addpath(s.path);
        pathsAdded = true;
    end

    % Default inputs
    Nargin = nargin;
    if Nargin < 8; params = []; end
    
    % Flag to force Pc calculations using the SDMC method
    params = set_default_param(params,'ForceSDMCCalculation',  false);
    params = set_default_param(params,'PreventSDMCCalculation',false);
    
    % Ensure 2D-Nc and 3D-Nc are calculated, if SDMC is forced
    if params.ForceSDMCCalculation
        params.ForceNc2DCalculation = true; params.PreventNc2DCalculation = false;
        params.ForceNc3DCalculation = true; params.PreventNc3DCalculation = false;
    end
    
    %% Call PcMultiStep to calculate the Pc2D, Nc2D and Nc3D methods
    
    % Call PcMultiStep
    [Pc,out] = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params);
    
    % Return if any data quality errors were detected
    if out.AnyDataQualityErrors
        return;
    end
    
    %% Check the SDMC run parameters
    
    % Default parameters for Pc_SDMC function
    params = out.MultiStepParams;
    params = set_default_param(params,'SDMCParams',[]);
    
    % Pc level to force using the SDMC method
    params = set_default_param(params,'PcToForceSDMCCalculation',1);
    
    % Flags to prevent SDMC method from running. 
    % Note: 2D-Nc or 3D-Nc cannot be prevented from running if SDMC is
    %       forced to run. In addition, a user cannot both prevent and
    %       force a Pc calculation.
    % SDMC must also be prevented if 3D-Nc is also prevented
    params = set_default_param(params,'PreventSDMCCalculation',params.PreventNc3DCalculation);
    
    % Forces SDMC to use the recommended number of trials instead of
    % limiting to max trials
    overrideForceRecTrials = false;
    if isfield(params,'ForceRecommendedTrials')
        overrideForceRecTrials = true;
    end
    params = set_default_param(params,'ForceRecommendedTrials',false);
    
    % Check the prevent/force parameters
    if params.PreventSDMCCalculation && params.ForceSDMCCalculation
        error('SDMC cannot be both prevented and forced to run.');
    end
    
    % Check if SDMC calculation is required
    if params.PreventSDMCCalculation
        NeedSDMCCalculation = false;
    elseif params.ForceSDMCCalculation || isnan(Pc)
        % Calculate SDMC if forced or if the output Pc value remains undefined
        NeedSDMCCalculation = true;
    elseif out.NeedNc3DCalculation && out.AnyNc3DViolations
        % Calculate SDMC if there were any 3D-Nc method usage violations
        NeedSDMCCalculation = true;
    elseif Pc >= params.PcToForceSDMCCalculation(1)
        % Calculate SDMC if semi-analytical Pc estimate is sufficiently large
        NeedSDMCCalculation = true;
    else
        NeedSDMCCalculation = false;
    end
    
    %% Perform the SDMC calculation, if required
    if NeedSDMCCalculation
        
        % Initialize parameters
        SDMCParams = params.SDMCParams;
        
        % Get the Pc value from the previous Pc calculations
        SDMCParams.InputPc = Pc;
        
        % Copy the ForceRecommendedTrials, if it was defined
        if overrideForceRecTrials
            SDMCParams.ForceRecommendedTrials = params.ForceRecommendedTrials;
        end
        
        % Make sure the generate_ca_dist_plot parameter gets transferred
        if isfield(params,'generate_ca_dist_plot')
            SDMCParams.generate_ca_dist_plot = params.generate_ca_dist_plot;
        end
        
        % Set the cross cov. correction parameters for the Nc3D calculation
        if params.apply_covXcorr_corrections
            [SDMCParams.apply_covXcorr_corrections, ...
             SDMCParams.covXcorr.sigp,SDMCParams.covXcorr.Gvecp, ...
             SDMCParams.covXcorr.sigs,SDMCParams.covXcorr.Gvecs] = ...
                get_covXcorr_parameters(params);
        else
            SDMCParams.apply_covXcorr_corrections = false;
        end
        
        % Calculate the SDMC time span and midpoint
        minPropTime = 3 * 60; % 3 minutes
        if ~isfield(SDMCParams,'span_days') || isempty(SDMCParams.span_days)
            if ~isnan(out.Nc3D) && out.Nc3DInfo.Converged
                % Time limits per Pc3D_Hall calculation
                Ta = out.Nc3DInfo.TaConj;
                Tb = out.Nc3DInfo.TbConj;
                TA = out.Nc3DInfo.Tmin_limit;
                TB = out.Nc3DInfo.Tmax_limit;

                % Expand the limits based on 3DNc violations
                if out.AnyNc3DViolations
                    expand = inf;
                else
                    expand = 3;
                end

                % Calculate the midpoint and time bounds, this generally will
                % be a short time region about the midpoint (which isn't always
                % near TCA). This calculation significantly saves on wasted
                % propagation performed by SDMC.
                tmid = (Tb+Ta)/2;
                dt = max(expand*max(Tb-Ta,1e-300),minPropTime)/2;
                ta = max(TA,tmid-dt);
                tb = min(TB,tmid+dt);
                SDMCParams.tmid = (tb+ta)/2;
                SDMCParams.span_days = (tb-ta)/2 / 86400;
            elseif ~isnan(out.Nc2D) && out.Nc2DInfo.Converged
                % Time limits per Pc2D_Hall computation
                Ta = out.Nc2DInfo.Ta;
                Tb = out.Nc2DInfo.Tb;

                % Expand the limits based on 2DNc violations
                if out.AnyNc2DViolations
                    expand = inf;
                else
                    expand = 3;
                end

                % Calculate the midpoint and time bounds, this generally will
                % be a short time region about the midpoint (which isn't always
                % near TCA). This calculation significantly saves on wasted
                % propagation performed by SDMC.
                SDMCParams.tmid = (Tb+Ta)/2;
                dt = max(expand*(Tb-Ta),minPropTime)/2;
                ta = SDMCParams.tmid-dt;
                tb = SDMCParams.tmid+dt;
                SDMCParams.span_days = (tb-ta)/2 / 86400;
            end
        else
            SDMCParams.tmid = 0;
        end
        
        % Set the retrograde orbit reorientation mode
        SDMCParams = set_default_param(SDMCParams, ...
            'RetrogradeReorientation',params.RetrogradeReorientation);

        % Save SDMC parameters
        out.SDMCParams = SDMCParams;
        
        % Don't try to call SDMC if Nc3D didn't converge and an invalid
        % covariance has been detected
        if (~out.Nc3DInfo.Converged && out.DataQualityError.invalidCov6x6)
            out.SDMCPc = nan;
        else
            % Calculate the SDMC value
            try
                % Update max attempted Pc calculation method
                % (1 = 2D-Pc, 2 = 2D-Nc, 3 = 3D-Nc, 4 = SDMC, 5 = BFMC)
                out.PcMethodMaxAnalytical = out.PcMethodMax;
                out.PcMethodMax = 4;
                [out.SDMCPc, out.SDMCInfo] = ...
                    Pc_SDMC(r1,v1,C1,r2,v2,C2,HBR,SDMCParams);
            catch ME
                if strcmp(ME.identifier,'SDMC:invalid_6x6_cov')
                    out.SDMCPc = nan;
                else
                    rethrow(ME);
                end
            end
        end
    
        % Process SDMC method usage violations
        if isnan(out.SDMCPc)

            % Convergence violation
            out.SDMCInfo.PcCalculated = false;
            out.AnySDMCViolations = true;

        else

            out.SDMCInfo.PcCalculated = true;

            % Get the encounter segment
            ta = out.SDMCInfo.ta;
            tb = out.SDMCInfo.tb;

            % Get the time points that are within 5% of the encounter
            % start/stop times
            pct5 = (tb - ta) * 0.05;
            t_min = ta + pct5;
            t_max = tb - pct5;

            % Check to see if any hits occur outside of the min/max limits
            trialData = out.SDMCInfo.trialData;
            if ~isempty(trialData)
                numHits = sum(trialData.hitIndicator == 1);
                if numHits > 0
                    trialHits = trialData(trialData.hitIndicator == 1,:);
                    numExtendedHits = sum(trialHits.hitTime < t_min) + ...
                        sum(trialHits.hitTime > t_max);
                else
                    numExtendedHits = 0;
                end
            else
                numExtendedHits = 0;
            end

            % Extended conjunction duration usage violation indicator
            out.SDMCInfo.Indicators.Extended = numExtendedHits;

            % Extended conjunction violation occurs if any extended hits are
            % encountered
            out.SDMCInfo.Violations.Extended = out.SDMCInfo.Indicators.Extended > 0;

            % Offset violation occurs if Nc3D had an offset violation
            out.SDMCInfo.Violations.Offset = ~isfield(out.Nc3DInfo, 'Violations') || ...
                out.Nc3DInfo.Violations.Offset;

            % Check if any violations occured
            out.AnySDMCViolations = out.SDMCInfo.Violations.Extended || ...
                out.SDMCInfo.Violations.Offset;
            
            % Check for a BFMC run to substitute for the SDMC run
            out.SDMCInfo.BFMCsubstitution = false;
            if isfield(params,'BFMCpath') && ~isempty(params.BFMCpath)
                % Update the max Pc estimation method attempted so far
                out.PcMethodMax = 5; % Indicates BFMC has been attempted
                % Process BFMC VCM mode data
                BFMClisfile = fullfile(params.BFMCpath,'bfmc_vcm.lis');
                BFMCoutfile = fullfile(params.BFMCpath,'bfmc_vcm.out');
                if ~exist(BFMClisfile,'file')
                    warning('bfmc_vcm.lis not found in specified path; not substituting BFMC data');
                elseif ~exist(BFMCoutfile,'file')
                    warning('bfmc_vcm.out not found in specified path; not substituting BFMC data');
                else
                    % Load BFMC conjunction data from .lis vile
                    BFMClisfile = fullfile(params.BFMCpath,'bfmc_vcm.lis');
                    [lisdata] = BFMC_lis_reader(BFMClisfile);
                    %
                    % % Test transformation from the BFMC ref. frame to
                    % the SDMC ref. frame: TEME2J2K
                    %
                    % % TCA epoch
                    % TCAUTC = datestr(lisdata.TCA,'yyyy-mm-dd HH:MM:SS.FFF');
                    %
                    % r1BFMC = [lisdata.PriX    lisdata.PriY    lisdata.PriZ   ] * 1e3;
                    % v1BFMC = [lisdata.PriXdot lisdata.PriYdot lisdata.PriZdot] * 1e3;
                    % 
                    % [r1TRAN,v1TRAN] = PosVelConvert(r1BFMC,v1BFMC, ...
                    %     TCAUTC,'TEME2J2K','4terms'); % 4terms found to be better than 106terms !!!!
                    % [norm(r1-r1TRAN) norm(v1-v1TRAN)]
                    %
                    % Load BFMC trial data from .out file
                    disp('--- LOADING BFMC SUBSTITUTION DATA...')
                    SDMCInfo = out.SDMCInfo;
                    BFMCdata = array2table(load(BFMCoutfile,"-ascii"));
                    BFMCdata.Properties.VariableNames = out.SDMCInfo.trialData.Properties.VariableNames;
                    % Trim BFMC data to a round number of trials
                    numTrials = numel(BFMCdata.hitIndicator);
                    if numTrials < 1e4
                        warning('--- Fewer than 1e4 BFMC trials; not substituting BFMC data');
                    else
                        strTrials = smart_exp_format(numTrials,4);
                        rndTrials = str2double(strTrials);
                        if rndTrials < numTrials
                            disp(['--- ROUNDING BFMC TRIALS FROM ' num2str(numTrials) ' to ' strTrials])
                            BFMCdata = BFMCdata(1:rndTrials,:);
                        end
                        % Substitute BFMC data for SDMC data
                        disp('--- PROCESSING BFMC SUBSTITUTION DATA')
                        out.SDMCInfo.BFMCsubstitution = true;
                        out.SDMCInfo.numTrials = numel(BFMCdata.hitIndicator);
                        out.SDMCInfo.numHits = sum(BFMCdata.hitIndicator);
                        ndx = BFMCdata.hitRadius <= 1.01*2*sqrt(2)*HBR;
                        out.SDMCInfo.trialData = TransformBFMCdata(BFMCdata(ndx,:),lisdata.TCA);
                        out.SDMCInfo.overallPicTrialData = TransformBFMCdata(BFMCdata(1:1e4,:),lisdata.TCA);
                        [~,Pvalue] = binomial_prop_test( ...
                            out.SDMCInfo.numHits/out.SDMCInfo.numTrials,out.SDMCInfo.numTrials, ...
                                SDMCInfo.numHits/    SDMCInfo.numTrials,    SDMCInfo.numTrials);
                        if Pvalue < 1e-3
                            warning(['Possible BFMC vs SDMC Pc difference: Pvalue = ' ...
                                smart_exp_format(Pvalue,3)]);
                        end
                        [out.SDMCInfo.Pc,out.SDMCInfo.PcUnc] = binofit( ...
                            out.SDMCInfo.numHits,out.SDMCInfo.numTrials);
                        out.SDMCPc = out.SDMCInfo.Pc;
                        out.SDMCInfo.Violations.Extended = 0;
                        out.SDMCInfo.Violations.Offset = 0;
                    end
                end
            end

        end

        % Update the reported Pc value
        if out.SDMCInfo.PcCalculated
            Pc = out.SDMCPc;
            if out.SDMCInfo.BFMCsubstitution
                out.PcMethod = 'BFMC';
                out.PcMethodNumAnalytical = out.PcMethodNum;
                out.PcMethodNum = 5;
            else
                out.PcMethod = 'SDMC';
                out.PcMethodNumAnalytical = out.PcMethodNum;
                out.PcMethodNum = 4;
            end
        end
        
    end
 
end

%% ========================================================================
function data = TransformBFMCdata(data,TCADN)

    % Transform BFMC state data into same ref. frame as SDMC state data
    % Specifically: TEME2J2K with 4terms (found to be better than 106terms)
    
    N = size(data,1);
    disp(['---- TRANSFORMING STATES FOR ' num2str(N) ' TRIALS']);
    
    % Extract states for parfor processing
    priPosX_km = data.priPosX_km; priVelX_kmps = data.priVelX_kmps;
    priPosY_km = data.priPosY_km; priVelY_kmps = data.priVelY_kmps;
    priPosZ_km = data.priPosZ_km; priVelZ_kmps = data.priVelZ_kmps;
    secPosX_km = data.secPosX_km; secVelX_kmps = data.secVelX_kmps;
    secPosY_km = data.secPosY_km; secVelY_kmps = data.secVelY_kmps;
    secPosZ_km = data.secPosZ_km; secVelZ_kmps = data.secVelZ_kmps;
    
    % Hit times relative to TCA in days
    hitTimeDays = data.hitTime/86400;
    
    % Transform the states (TEME2J2K)
    parfor n=1:N
        % Epoch for this trial
        EpochUTC = datestr(TCADN+hitTimeDays(n),'yyyy-mm-dd HH:MM:SS.FFF');
        % Primary TCA states
        [r,v] = PosVelConvert( ...
            [priPosX_km(n)   priPosY_km(n)   priPosZ_km(n)  ], ...
            [priVelX_kmps(n) priVelY_kmps(n) priVelZ_kmps(n)], ...
            EpochUTC,'TEME2J2K','4terms');
        priPosX_km(n)   = r(1);
        priPosY_km(n)   = r(2);
        priPosZ_km(n)   = r(3);
        priVelX_kmps(n) = v(1);
        priVelY_kmps(n) = v(2);
        priVelZ_kmps(n) = v(3);
        % Secondary TCA states
        [r,v] = PosVelConvert( ...
            [secPosX_km(n)   secPosY_km(n)   secPosZ_km(n)  ], ...
            [secVelX_kmps(n) secVelY_kmps(n) secVelZ_kmps(n)], ...
             EpochUTC,'TEME2J2K','4terms');
        secPosX_km(n)   = r(1);
        secPosY_km(n)   = r(2);
        secPosZ_km(n)   = r(3);
        secVelX_kmps(n) = v(1);
        secVelY_kmps(n) = v(2);
        secVelZ_kmps(n) = v(3);
    end
    
    % Store transformed states
    data.priPosX_km = priPosX_km; data.priVelX_kmps = priVelX_kmps;
    data.priPosY_km = priPosY_km; data.priVelY_kmps = priVelY_kmps;
    data.priPosZ_km = priPosZ_km; data.priVelZ_kmps = priVelZ_kmps;
    data.secPosX_km = secPosX_km; data.secVelX_kmps = secVelX_kmps;
    data.secPosY_km = secPosY_km; data.secVelY_kmps = secVelY_kmps;
    data.secPosZ_km = secPosZ_km; data.secVelZ_kmps = secVelZ_kmps;
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2023-AUG-21 | Initial Development.
% L. Baars       | 2023-AUG-24 | Fixed invalidCov6x6 field reference
% L. Baars       | 2024-JAN-11 | Standardized "out" structure
% D. Hall        | 2024-SEP-16 | Added support for BFMC data processing
% D. Hall        | 2024-DEC-11 | Added PcMethod*Analytical outputs and
%                                adjusted PcMethod* outputs.
% L. Baars       | 2025-AUG-22 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
