function [Pc, out] = Pc_SDMC(r1,v1,C1,r2,v2,C2,HBR,params)
% Pc_SDMC - Calculate the Simple Dynamics Monte Carlo (SDMC) approximation
%           for the probability of collision between two satellites, given
%           input states and covariances at the nominal time of closest
%           approach (TCA).
%
% Syntax: [Pc, out] = Pc_SDMC(r1,v1,C1,r2,v2,C2,HBR);
%         [Pc, out] = Pc_SDMC(r1,v1,C1,r2,v2,C2,HBR,params);
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
%   The SDMC method locally approximates both satellite trajectories using
%   either two-body or rectilinear equations of motion, depending on user
%   inputs. Being a Monte Carlo approach, the method will run a set number
%   of trials with small deviations to the object states. Each trial will
%   determine if the deviated states yield a close-approach event which
%   results in a miss distance less than the combined hard-body radius
%   (HBR) of the conjunction (a hit). The sum of the number of hits divided
%   by the total number of trials is the SDMC Pc estimate. This method is
%   fully explained in Hall (2018).
%
%   The majority of the processing is performed within a compiled
%   Linux-only shared object library, libsdmctask.so. Pc_SDMC.m servers as
%   a Matlab wrapper around the library, converting the inputs/outputs of
%   the library between the Matlab and library expected formats. Since the
%   library is a Linux-only library, the wrapper will exit with errors if
%   it is not run on Linux.
%
%   All input/output units are meters/seconds/radians, unless otherwise
%   noted.
%
%   IMPORTANT NOTE:
%     In order to run call_SDMC, a system environment variable must be set
%     using the absolute path to the following path:
%
%       SDK/DistributedMatlab/ProbabilityOfCollision/SDMC_Utils/lib
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
% =========================================================================
%
% Input:
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
%    params  - (Optional) Auxilliary input parameter structure, described
%              in detail in function "default_params_Pc_SDMC".
%
% =========================================================================
%
% Output:
%
%  Pc - The estimated conjunction Pc value.
%
%  out - An auxilliary output structure which contains quantities from the
%        SDMC method calculation, within the following fields (listed
%        roughly in decreasing order of importance):
%
%    Pc = SDMC estimated Pc                                           [1x1]
%
%    PcUnc = Uncertainty of the Pc to the confidence level of         [1x2]
%            params.conf_level
%
%    numHits = The number of Monte Carlo trials with close approach   [1x1]
%              events less than the HBR
%
%    numTrials = The total number of Monte Carlo trials run           [1x1]
%
%    numWorkers = The number of workers used in the parallel pool     [1x1]
%                 when SDMC was called
%
%    maxOutputTrials = The maximum number in trialData                [1x1]
%
%    trialData = Trial specific data for all trials with a CA       [table]
%                distance less than or equal to max_radius. See
%                call_SDMC.m documentation for a full description.
%
%    overallPicTrialData = Trial specific data for 10,000 trials    [table]
%                          with a CA distance less than or equal to
%                          1e9 meters. This field is only filled in
%                          when params.generate_ca_dist_plot is set
%                          to true. See call_SDMC.m documentation
%                          for a full description.
%
%    covXcorr_corrections_applied = Indicates if covariance cross-    [1x1]
%                                   correlation corrections were
%                                   applied.
%
%    r1 = Primary object's ECI position vector (m)             [3x1 or 1x3]
%
%    v1 = Primary object's ECI velocity vector (m/s)           [3x1 or 1x3]
%
%    C1 = Primary object's ECI covariance matrix (m^2, m^2/s,         [6x6]
%         m^2/s^2)
%
%    r2 = Secondary object's ECI position vector (m)           [3x1 or 1x3]
%
%    v2 = Secondary object's ECI velocity vector (m/s)         [3x1 or 1x3]
%
%    C2 = Secondary object's ECI covariance matrix (m^2, m^2/s,       [6x6]
%         m^2/s^2)
%
%    HBR = Combined primary+secondary hard-body radii (m)             [1x1]
%
%    params = A copy of the parameters structure used for the      [struct]
%             calculation
%
%    RetrogradeReorientation = Boolean indicating whether or not the  [1x1]
%                              oribt was reoriented
%
%    ta,tb = Start/stop propagation bounds (in seconds from           [1x1]
%            params.tmid) for the SDMC run
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
% Initial version: Feb 2023;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------

    %% Add necessary paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'Utils')); addpath(s.path);
        s = what(fullfile(p,'SDMC_Utils')); addpath(s.path);
        s = what(fullfile(p,'../Utils/CovarianceTransformations')); addpath(s.path);
        s = what(fullfile(p,'../Utils/LoggingAndStringReporting')); addpath(s.path);
        s = what(fullfile(p,'../Utils/OrbitTransformations')); addpath(s.path);
        pathsAdded = true;
    end

    %% Check inputs and parameters
    
    % Initializations and defaults
    Nargin = nargin;

    % Non-verbose processing by default
    if (Nargin < 8) || isempty(params)
        params.verbose = false;
    end

    % Set any parameters that remain to be defined to their defaults, but leave
    % any that are already defined unchanged.
    params = default_params_Pc_SDMC(params);
    
    % Copy parameters to the output structure
    out.params = params;
    
    % Process the input HBR
    N_HBR = numel(HBR);

    if (N_HBR == 2)
        % Process individual primary and secondary HBR values
        %  HBR(1) = Primary   HBR
        %  HBR(2) = Secondary HBR
        if (min([HBR(1) HBR(2)]) < 0)
           error('Both HBR values must be nonnegative.');
        end  
        HBR = HBR(1)+HBR(2);
    elseif (N_HBR ~= 1)
        error('Input HBR must have one or two elements.');
    end

    if (HBR <= 0)
        error('Combined HBR value must be positive.');
    end

    out.HBR = HBR;
    
    % Make the max radius slightly bigger than the HBR so that we know the
    % hit/miss data will include all hits
    params.max_radius = HBR * 1.001;
    
    %% Retrograde orbit processing

    % Reorient to handle conjunctions involving retrograde orbits, which occur
    % for the Alfano (2009) test cases, but have never been encountered with
    % actual CARA conjunctions. This eliminates the occurence of retrograde
    % orbits by reorienting the (x,y,z) axes of the reference frame.
    if params.RetrogradeReorientation > 0
        [r1,v1,C1,r2,v2,C2,RRout] = ...
            RetrogradeReorientation(r1,v1,C1,r2,v2,C2,params);
        out.RetrogradeReorientation = RRout.Reoriented;
    else
        out.RetrogradeReorientation = false;
    end
    
    % Copy states into output structure
    out.r1 = r1;
    out.v1 = v1;
    out.C1 = C1;
    out.r2 = r2;
    out.v2 = v2;
    out.C2 = C2;
    
    %% Get covariance cross correlation parameters
    if params.apply_covXcorr_corrections
        [XCprocessing,sigp,Gp,sigs,Gs] = get_covXcorr_parameters(params);
    else
        XCprocessing = false;
    end
    if XCprocessing
        % Reorient the 1x6 sensitivity vectors, if required
        if out.RetrogradeReorientation
            Gp = (RRout.M6 * Gp')';
            Gs = (RRout.M6 * Gs')';
        end
    else
        % Null DCP sigma values and sensitivity vectors
        sigp = 0;
        sigs = 0;
        Gp = zeros(1,6);
        Gs = zeros(1,6);
    end
    out.covXcorr_corrections_applied = XCprocessing;
    
    params.apply_covXcorr_corrections = XCprocessing;
    params.covXcorr.sigp = sigp;
    params.covXcorr.Gvecp = Gp;
    params.covXcorr.sigs = sigs;
    params.covXcorr.Gvecs = Gs;
    
    %% Calculate the number of trials needed
    if isempty(params.num_trials)
        inputPc = params.InputPc;
        if isempty(inputPc) || isnan(inputPc) || inputPc == 0
            if params.warning_level > 0
                warning(['Could not determine Pc; setting number of trials to ' ...
                    smart_exp_format(params.default_num_trials)]);
            end
            params.num_trials = params.default_num_trials;
        else
            if params.verbose
                disp(['InputPc = ' num2str(inputPc,'%.3e')]);
            end
            if numel(params.Target95pctPcAccuracy) == 1
                % Process single target accuracy
                params.num_trials = EstimateRequiredSamples(inputPc, ...
                    params.Target95pctPcAccuracy, 0.95);
            else
                % Procsess using target accuracy table
                log10Pc = log10(inputPc);
                if log10Pc <= params.Target95pctPcAccuracy(1,1)
                    tgtacc = params.Target95pctPcAccuracy(1,2);
                elseif log10Pc >= params.Target95pctPcAccuracy(end,1)
                    tgtacc = params.Target95pctPcAccuracy(end,2);
                else
                    tgtacc = interp1(                      ...
                        params.Target95pctPcAccuracy(:,1), ...
                        params.Target95pctPcAccuracy(:,2), ...
                        log10Pc);
                end
                if isnan(tgtacc) || tgtacc <= 0 || tgtacc >= 1
                    error('Target accuracy table processing failure');
                end
                params.num_trials = EstimateRequiredSamples(inputPc, ...
                    tgtacc, 0.95);
            end
            if params.num_trials > params.max_num_trials
                if isfield(params,'ForceRecommendedTrials') && params.ForceRecommendedTrials
                    if params.warning_level > 0
                        warning(['Forcing the recommended number of trials of ' ...
                            smart_exp_format(params.num_trials) newline ...
                            'SDMC may take a long time to run!']);
                    end
                else
                    if params.warning_level > 1
                        warning(['Calculated number of trials exceeds max. allowed; limiting to ' ...
                            smart_exp_format(params.max_num_trials)]);
                    end
                    params.num_trials = params.max_num_trials;
                end
            end
        end
    end
    
    %% Get the SDMC time span and midpoint
    if isempty(params.span_days) || isempty(params.tmid)
        % Calculate the orbital periods of the primary and secondary
        % since we couldn't get a more accurate time span from data passed
        % in
        period1 = orbit_period(r1,v1);
        period2 = orbit_period(r2,v2);
        periodMin = min(period1,period2);
        halfOrbitPeriodMin = periodMin / 2;
        tmid = 0;
        ta = -halfOrbitPeriodMin / 2;
        tb = halfOrbitPeriodMin / 2;
        params.span_days = (tb-ta)/2 / 86400;
    else
        tmid = params.tmid;
    end
    % Override the TCA to be the midpoint time. This helps to reduce the
    % amount of SDMC propagation when Nc3D time bounds are used.
    params.TCA = params.TCA + tmid/86400;
    out.ta = tmid - params.span_days * 86400;
    out.tb = tmid + params.span_days * 86400;
    
    %% Convert several components to RIC coordinates for use with SDMC
    % Convert the covariance matrices
    C1_RIC = ECI2RIC(C1,r1,v1);
    C2_RIC = ECI2RIC(C2,r2,v2);
    
    % Calculate the RIC coordinates of covariance cross-correlation vectors
    if XCprocessing
        r1vecECI = Gp(1:3);
        v1vecECI = Gp(4:6);
        r2vecECI = Gs(1:3);
        v2vecECI = Gs(4:6);
        r1vecRIC = ECI2RIC(r1vecECI,r1,v1);
        v1vecRIC = ECI2RIC(v1vecECI,r1,v1);
        r2vecRIC = ECI2RIC(r2vecECI,r2,v2);
        v2vecRIC = ECI2RIC(v2vecECI,r2,v2);
        params.covXcorr.Gvecp = [r1vecRIC v1vecRIC];
        params.covXcorr.Gvecs = [r2vecRIC v2vecRIC];
    end

    %% Run SDMC
    if ~params.generate_ca_dist_plot
        % Generic SDMC run
        tic
        params.max_output_trials = 1e4;
        [Pc,PcUnc,numHits,trialData,numWorkers] = ...
            call_SDMC_parallel(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,params);
        if params.verbose
            toc
        end
    else
        % Run SDMC twice, once to get the zoomed in view and another time
        % to get the overall picture. This mode overrides the max_raidus,
        % max_output_trials, and num_trials parameters. Negative max_radius
        % values indicate that position/velocity data for the primary and
        % secondary will be added to the trialData table. This data is
        % required in order to generate CA distribution plots.
        
        % First, run with 10,000 trials and a really large max radius in
        % order to get the overall picture
        tic
        origNumTrials = params.num_trials;
        params.num_trials = 1e4;
        params.max_output_trials = params.num_trials;
        params.max_radius = -1e9;
        [Pc,PcUnc,numHits,overallPicTrialData,numWorkers] = ...
            call_SDMC_parallel(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,params);
        % Fix trial data for time offsets imposed by tmid
        if ~isempty(overallPicTrialData)
            overallPicTrialData.hitTime = overallPicTrialData.hitTime + tmid;
            overallPicTrialData.hitTimeMiss = overallPicTrialData.hitTimeMiss + tmid;
            overallPicTrialData.pcaTime = overallPicTrialData.pcaTime + tmid;
        end
        % Now run with the original number of trials and a max radius of
        % 3*HBR, this allows the CA distribution plot to produce a zoomed
        % in view around the CA distribution. Slightly increased the
        % max_num_output_trials in order to make sure we get all the hits
        % and misses around the HBR.
        if origNumTrials > params.num_trials
            % Only run another trial if the number of trials needed exceeds
            % the number of trials from the overall pic (this will happen
            % most of the time)
            params.num_trials = origNumTrials;
            params.max_output_trials = 5e4; % This value cannot exceed 50,000
            % Max radius is just a little larger than 1/2 the diagonal of a
            % square with a side length of HBR*4
            params.max_radius = -HBR*2*sqrt(2)*1.001;
            if params.verbose
                disp('Calling SDMC...');
            end
            [Pc,PcUnc,numHits,trialData,numWorkers] = ...
                call_SDMC_parallel(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,params);
            % A special case occurs when the number of hits does not match
            % the hitIndicators within the trialData. This can occur when
            % the HBR is really large (such as the ISS) and the number of
            % trials is large. In order to fix this, we have to rerun SDMC
            % with a much smaller HBR and then replace the original hits
            % with new hits and then also replace some misses with the
            % remaining hits.
            if ~isempty(trialData) && numHits > 0 && numHits ~= sum(trialData.hitIndicator)
                % Calculate the number of hits that were missing from the
                % original run
                numMissing = numHits - sum(trialData.hitIndicator);
                % Remove the original hits from the trialData
                trialData(trialData.hitIndicator == 1,:) = [];
                % Randomly remove misses from the original run equal to the
                % numMissing parameter
                removeIdx = randperm(height(trialData),numMissing);
                trialData(removeIdx,:) = [];
                % Rerun SDMC (using the original seed number to ensure we
                % have the same outputs) but with a much smaller
                % max_radius. This should ensure that we don't get missing
                % hits data.
                params.max_radius = -HBR*1.001;
                [~,~,numHits2,trialData2,~] = ...
                    call_SDMC_parallel(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,params);
                % Remove any misses from this new trial
                trialData2(trialData2.hitIndicator == 0,:) = [];
                % Make sure that the current and previous number of hits
                % are the same
                if numHits ~= numHits2
                    % We should never get here
                    error('Random number generator trials did not match!');
                end
                % Provide a warning if we still to have a match
                if numHits ~= height(trialData2)
                    if params.warning_level > 0
                        warning('Detailed trial data numbers do not match the summary number of hits, SDMC plotting will be incorrect!');
                    end
                end
                % Append the hit data to the previous trial data (which
                % does not have any hits in it)
                trialData = [trialData; trialData2];
                mixOrder = randperm(height(trialData))';
                trialData = trialData(mixOrder,:);
            end
            if params.verbose
                toc
            end
        else
            % Use the overall pic data if it has better accuracy than the
            % original number of trials
            trialData = overallPicTrialData;
            % Unadjust the offsets because they will be fixed below
            trialData.hitTime = trialData.hitTime - tmid;
            trialData.hitTimeMiss = trialData.hitTimeMiss - tmid;
            trialData.pcaTime = trialData.pcaTime - tmid;
        end
    end
    
    %% Fix trial data for time offsets imposed by tmid
    if ~isempty(trialData)
        trialData.hitTime = trialData.hitTime + tmid;
        trialData.hitTimeMiss = trialData.hitTimeMiss + tmid;
        trialData.pcaTime = trialData.pcaTime + tmid;
    end
    
    %% Copy the output data to the output data structure
    out.Pc = Pc;
    out.PcUnc = PcUnc;
    out.numHits = numHits;
    out.numTrials = params.num_trials;
    out.numWorkers = numWorkers;
    out.maxOutputTrials = params.max_output_trials;
    out.trialData = trialData;
    if params.generate_ca_dist_plot
        out.overallPicTrialData = overallPicTrialData;
    end
    
    if params.verbose
        disp(['SDMC Hits/Trials = ' smart_exp_format(numHits) ...
            '/' smart_exp_format(params.num_trials)]);
        disp(['SDMC Pc = ' smart_exp_format(Pc,params.PcNsf,[false true])]);
        disp(['SDMC 95% Pc Conf = ' ...
            smart_exp_format(PcUnc(1),params.PcNsf,[false true]) ' to ' ...
            smart_exp_format(PcUnc(2),params.PcNsf,[false true])]);
    end
end

%% Function which calls SDMC in parallel (when allowed)
function [Pc,PcUnc,numHits,trialData,numWorkers] = call_SDMC_parallel(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,params)
    % Check to see if parallel runs should be used
    runParallel = true;
    if params.num_trials <= 1e4 || ~params.use_parallel
        runParallel = false;
        numWorkers = 1;
    else
        try
            % If this call fails, then either the parallel toolbox is not
            % licensed or something went wrong with starting the parallel
            % pool. In either case we can't run in parallel.
            p = gcp;
            numWorkers = p.NumWorkers;
        catch
            runParallel = false;
            numWorkers = 1;
            warning('Parallel processing is not available, Pc_SDMC can take a very long time to run');
        end
    end
    
    if ~runParallel
        if ~isempty(params.num_workers)
            warning('num_workers parameter is not used when parallel processing is not available');
        end
        % Just run SDMC serially
        [Pc,PcUnc,numHits,trialData] = ...
            call_SDMC(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,params);
    else
        % Run SDMC in parallel, splitting the samples evenly among workers
        if ~isempty(params.num_workers)
            poolsize = params.num_workers;
            numWorkers = params.num_workers;
        else
            poolsize = p.NumWorkers;
        end
        numPerWorker = floor(params.num_trials/poolsize);
        if params.num_trials/poolsize ~= numPerWorker
            numPerWorker = numPerWorker + 1;
        end
        nSamp = ones(poolsize,1) * numPerWorker;
        % Make sure the number of samples is exactly equal to the
        % num_trials passed in
        sum_nSamp = sum(nSamp);
        if sum_nSamp > params.num_trials
            amountOver = sum_nSamp - params.num_trials;
            nSamp(1:amountOver) = nSamp(1:amountOver) - 1;
        end
        % Set the seed for Matlab's random number generator (used here
        % because external calls are using random number generators)
        rng(abs(params.seed));
        % The seeds used by SDMC are integers in the range -1 to
        % -2147483648
        seeds = (params.seed:-1:params.seed-poolsize)';
        idx = seeds < -2147483648;
        seeds(idx) = seeds(idx) + 2147483648;
        % Setup output variables for each worker
        wNumHits = nan(poolsize,1);
        wTrialData = cell(poolsize,1);
        % Run SDMC across all workers
        parfor i=1:poolsize
            tempParams = params;
            tempParams.num_trials = nSamp(i);
            tempParams.seed = seeds(i);
            [~,~,wNumHits(i),wTrialData{i}] = ...
                call_SDMC(r1,v1,C1_RIC,r2,v2,C2_RIC,HBR,tempParams);
        end
        % Determine the Pc and Pc Uncertainty based on the total number of
        % hits and total number of trials
        numHits = sum(wNumHits);
        conf_alpha = 1-params.conf_level;
        [Pc,PcUnc] = binofit(numHits,params.num_trials,conf_alpha);
        % Combine the trial data
        trialData = [];
        for i=1:poolsize
            if ~isempty(wTrialData{i})
                trialData = cat(1,trialData,wTrialData{i});
            end
        end
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
% L. Baars       | 2023-02-10  | Initial Development
% D. Hall        | 2023-07-27  | Added retrograde orbit processing
% L. Baars       | 2023-08-31  | Changed the default number of trials to
%                                process when an analytical Pc cannot be
%                                calculated.
% L. Baars       | 2024-05-28  | Added note to the description of this
%                                function about system environment
%                                variables which need to be set.
% L. Baars       | 2024-09-23  | Updated the add path section to include a
%                                missing path.
% D. Hall        | 2025-04-10  | Added check to determine if the structure
%                                "overallPicTrialData" is empty before
%                                adjusting contents, to avoid error
% L. Baars       | 2025-08-25  | Updated for public release of code.
% L. Baars       | 2025-09-03  | Added support for num_workers parameter.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
