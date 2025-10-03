function [figInfo] = Pc_Time_Plot(params, pcInfo)
% Pc_Time_Plot - Produces temporal plots of several Pc theories against
%                each other.
%
% Syntax: [figInfo] = Pc_Time_Plot(params, pcInfo);
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
%   This function creates temporal plots of several Pc theories against
%   each other. The outputs that are generated include two plots and some
%   informational text.
%
%   The first plot that is generated graphs the cummulative Pc against
%   time from TCA. The Pc theories available for this plot include:
%     PcCircle  - labeled as "2D-Pc Method"
%     Pc2D_Hall - labeled as "2D-Nc Method"
%     Pc3D_Hall - labeled as "3D-Nc Method"
%     Pc_SDMC   - labeled as "SDMC Method" and "MC 95% confidence"
%
%   Note that a confidence interval is included for Pc_SDMC since those
%   values can be determined from the results of Monte Carlo processing.
%
%   The second plot that is generated graphs the Pc rate against time from
%   TCA. The Pc theories available for this plot include Pc3D_Hall and
%   Pc_SDMC. In addition to these theories, vertical lines are also
%   included depicting the conjunction bounds as determined by Pc2D_Hall
%   (labeled as "2D-Nc") and Pc3D_Hall (labeled as "3D-Nc").
%
%   The informational text includes general conjunction information,
%   calculated Pc values, Pc usage violation checks, and conjunction bounds
%   information. If a usage violation check fails for a particular Pc
%   calculation theory, both the violation check and the reported Pc value
%   will be displayed in a red font.
%
%   It is expected that all Pc calculations have occurred before calling
%   this function and these calculations have been used to populate the
%   pcInfo structure with the needed data.
%
%   Note: "BFMC" is a tool used internal to NASA CARA and will not be
%         publicly released. References to BFMC within the code should be
%         ignored by external users. Minimal documentation will be provided
%         for any BFMC functionality which may appear within the code.
%
% =========================================================================
%
% Input:
%
%   params - Input parameter structure with fields defined from the
%            default_params_Pc_Time_Plot.m function. See that function for
%            a full description of the parameters.
%
%   pcInfo - Input structure defining information specific to the
%            conjunction to be graphed. A portion of the fields included in
%            the structure are optional based on the graphing parameters
%            supplied in the params structure. However, in general, the
%            "out" structure provided by PcMultiStepWithSDMC.m should be
%            used as pcInfo, with the following fields added to the
%            structure:
%
%       r1  - Primary object's ECI position vector (m)         [3x1 or 1x3]
%
%       v1  - Primary object's ECI velocity vector (m/s)       [3x1 or 1x3]
%
%       C1  - Primary object's ECI covariance matrix                  [6x6]
%
%       r2  - Secondary object's ECI position vector (m)       [3x1 or 1x3]
%
%       v2  - Secondary object's ECI velocity vector (m/s)     [3x1 or 1x3]
%
%       C2  - Secondary object's ECI covariance matrix                [6x6]
%
% =========================================================================
%
% Output:
%
%   figInfo - Information pertaining to the figure generated. If
%             params.plot_save_loc is true, then this value will contain
%             the name of the file that was saved (including the full file
%             path), otherwise, it will contain a reference to the figure
%             handle for the plot that was created.
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    %% Setup paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        addpath(fullfile(p,'../LoggingAndStringReporting'));
        addpath(fullfile(p,'Utils'));
        pathsAdded = true;
    end
    
    %% Check Input Parameters
    figInfo = [];
    % Load the default parameters
    params = default_params_Pc_Time_Plot(params);
    % Check if a time plot needs to be generated
    if ~params.generate_time_plot
        return;
    end
    
    allFieldsExist = true;
    
    % Check for conjunction parameters
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'r1') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'v1') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'C1') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'r2') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'v2') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'C2') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'HBR') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'AnyDataQualityErrors') && allFieldsExist;
    allFieldsExist = dispMissingField(pcInfo.DataQualityError, 'pcInfo.DataQualityError', 'invalidCov6x6') && allFieldsExist;
    if ~allFieldsExist
        error('Cannot plot Pc_Time_Plot since basic conjunction parameters are missing');
    end
    
    % Check for a data quality error
    if pcInfo.AnyDataQualityErrors && ~pcInfo.DataQualityError.invalidCov6x6
        warning('Conjunction had data quality errors, cannot plot results!');
        return;
    end
    
    % Check required parameters for plotting PcCircle data
    if ~isfield(params,'plot_PcCircle'); params.plot_PcCircle = []; end
    if isempty(params.plot_PcCircle); params.plot_PcCircle = false; end
    if params.plot_PcCircle
        allFieldsExist = true;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Pc2D') && allFieldsExist;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Pc2DViolations') && allFieldsExist;
        if isfield(pcInfo,'Pc2DViolations')
            allFieldsExist = dispMissingField(pcInfo.Pc2DViolations, 'pcInfo.Pc2DViolations', 'NPD') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Pc2DViolations, 'pcInfo.Pc2DViolations', 'Extended') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Pc2DViolations, 'pcInfo.Pc2DViolations', 'Offset') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Pc2DViolations, 'pcInfo.Pc2DViolations', 'Inaccurate') && allFieldsExist;
        end
        if ~allFieldsExist
            % error('Cannot plot PcCircle data due to missing input data!');
            warning('Cannot plot PcCircle data due to missing input data; not making time plot');
            params.plot_PcCircle = false;
        end
    end
    
    % Check required parameters for plotting Pc2D_Hall data
    if ~isfield(params,'plot_Pc2D_Hall'); params.plot_Pc2D_Hall = []; end
    if isempty(params.plot_Pc2D_Hall); params.plot_Pc2D_Hall = false; end
    if params.plot_Pc2D_Hall
        allFieldsExist = true;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Nc2D') && allFieldsExist;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Nc2DInfo') && allFieldsExist;
        if isfield(pcInfo,'Nc2DInfo')
            allFieldsExist = dispMissingField(pcInfo.Nc2DInfo, 'pcInfo.Nc2DInfo', 'Violations') && allFieldsExist;
            if isfield(pcInfo.Nc2DInfo,'Violations')
                allFieldsExist = dispMissingField(pcInfo.Nc2DInfo.Violations, 'pcInfo.Nc2DInfo.Violations', 'Extended') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.Nc2DInfo.Violations, 'pcInfo.Nc2DInfo.Violations', 'Offset') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.Nc2DInfo.Violations, 'pcInfo.Nc2DInfo.Violations', 'Inaccurate') && allFieldsExist;
            end
            allFieldsExist = dispMissingField(pcInfo.Nc2DInfo, 'pcInfo.Nc2DInfo', 'TMeanRate') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc2DInfo, 'pcInfo.Nc2DInfo', 'TSigmaRate') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc2DInfo, 'pcInfo.Nc2DInfo', 'TQ0min') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc2DInfo, 'pcInfo.Nc2DInfo', 'SigmaQ0min') && allFieldsExist;
        end
        if ~allFieldsExist
            if pcInfo.DataQualityError.invalidCov6x6
                warning('Invalid 6x6 covariance caused missing Pc2D_Hall parameters, turning off Pc2D_Hall plot');
                params.plot_Pc2D_Hall = false;
            else
                % error('Cannot plot Pc2D_Hall data due to missing input data!');
                warning('Cannot plot Pc2D_Hall data due to missing input data; not making time plot');
                params.plot_Pc2D_Hall = false;
            end
        end
    end
    
    % Check required parameters for plotting Pc3D_Hall data
    if ~isfield(params,'plot_Pc3D_Hall'); params.plot_Pc3D_Hall = []; end
    if isempty(params.plot_Pc3D_Hall); params.plot_Pc3D_Hall = false; end
    if params.plot_Pc3D_Hall
        allFieldsExist = true;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Nc3D') && allFieldsExist;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Nc3DInfo') && allFieldsExist;
        if isfield(pcInfo,'Nc3DInfo')
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'Violations') && allFieldsExist;
            if isfield(pcInfo.Nc3DInfo,'Violations')
                allFieldsExist = dispMissingField(pcInfo.Nc3DInfo.Violations, 'pcInfo.Nc3DInfo.Violations', 'Extended') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.Nc3DInfo.Violations, 'pcInfo.Nc3DInfo.Violations', 'Offset') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.Nc3DInfo.Violations, 'pcInfo.Nc3DInfo.Violations', 'Inaccurate') && allFieldsExist;                
            end
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'Teph') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'Nc') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'Ncdot') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'Nccum') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'TaConj') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Nc3DInfo, 'pcInfo.Nc3DInfo', 'TbConj') && allFieldsExist;
        end
        if ~allFieldsExist
            if pcInfo.DataQualityError.invalidCov6x6
                warning('Invalid 6x6 covariance caused missing Pc3D_Hall parameters, turning off Pc3D_Hall plot');
                params.plot_Pc3D_Hall = false;
            else
                % error('Cannot plot Pc3D_Hall data due to missing input data!');
                warning('Cannot plot Pc3D_Hall data due to missing input data; not making time plot');
                params.plot_Pc3D_Hall = false;
            end
        end
    end
    
    % Check required parameters for plotting Pc_SDMC data
    if ~isfield(params,'plot_Pc_SDMC'); params.plot_Pc_SDMC = []; end
    if isempty(params.plot_Pc_SDMC); params.plot_Pc_SDMC = false; end
    if params.plot_Pc_SDMC
        allFieldsExist = true;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'SDMCPc') && allFieldsExist;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'SDMCInfo') && allFieldsExist;
        if isfield(pcInfo,'SDMCInfo')
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'Violations') && allFieldsExist;
            if isfield(pcInfo.SDMCInfo,'Violations')
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.Violations, 'pcInfo.SDMCInfo.Violations', 'Extended') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.Violations, 'pcInfo.SDMCInfo.Violations', 'Offset') && allFieldsExist;
            end
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'Pc') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'PcUnc') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'numHits') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'numTrials') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'maxOutputTrials') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'trialData') && allFieldsExist;
            if isfield(pcInfo.SDMCInfo,'trialData') && ~isempty(pcInfo.SDMCInfo.trialData)
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'hitIndicator') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'hitTime') && allFieldsExist;
            end
        end
        if ~allFieldsExist
            if pcInfo.DataQualityError.invalidCov6x6
                warning('Invalid 6x6 covariance caused missing Pc_SDMC parameters, turning off Pc_SDMC plot');
                params.plot_Pc_SDMC = false;
            else
                % error('Cannot plot Pc_SDMC data due to missing input data!');
                % warning('Cannot plot Pc_SDMC data due to missing input data, not making time plot');
                params.plot_Pc_SDMC = false;
            end
        end
    end
    
    % Check to make sure at least one individual plot parameter has been
    % turned on
    if ~params.plot_PcCircle && ~params.plot_Pc2D_Hall && ...
            ~params.plot_Pc3D_Hall && ~params.plot_Pc_SDMC
        warning('No plots generated since none of the individual plot parameters were turned on!');
        return;
    end
    
    %% Prepare basic plot data
    
    % Calculate the alpha confidence value
    conf_alpha = 1 - params.conf_level;
    
    % General conjunction parameters
    [conjID, conjInfoStr, conjInfoStr2] = GetConjInfo(params, pcInfo);
    if ~isempty(params.plot_conjID_string)
        conjID = params.plot_conjID_string;
    end
    
    Pc2D = [];
    Pc2DViolations = {};
    numPc2DViolations = 0;
    if params.plot_PcCircle
        % Get Pc and Violation text
        Pc2D = pcInfo.Pc2D;
        maxNumViolations = 4;
        Pc2DViolations = cell(maxNumViolations,1);
        Pc2DViolations{1} = ['NPD Check:        ' getPassedFailed(pcInfo.Pc2DViolations.NPD)];
        Pc2DViolations{2} = ['Extended Check:   ' getPassedFailed(pcInfo.Pc2DViolations.Extended)];
        Pc2DViolations{3} = ['Offset Check:     ' getPassedFailed(pcInfo.Pc2DViolations.Offset)];
        Pc2DViolations{4} = ['Inaccuracy Check: ' getPassedFailed(pcInfo.Pc2DViolations.Inaccurate)];
        for i = 1:maxNumViolations
            if endsWith(Pc2DViolations{i},'Failed')
                numPc2DViolations = numPc2DViolations + 1;
            end
        end
        
        if params.verbose
            disp(params.verboseTextSep);
            disp(['PcCircle Pc = ' smart_exp_format(Pc2D,params.PcNsf,[false true])]);
            disp('PcCircle Violations:');
            for i = 1:maxNumViolations
                disp(['  ' Pc2DViolations{i}]);
            end
        end
    end
    
    Nc2D = [];
    Nc2DViolations = {};
    Ta_Nc2D = [];
    Tb_Nc2D = [];
    numNc2DViolations = 0;
    if params.plot_Pc2D_Hall
        % Get Pc and Violation text
        Nc2D = pcInfo.Nc2D;
        maxNumViolations = 3;
        Nc2DViolations = cell(maxNumViolations,1);
        Nc2DViolations{1} = ['Extended Check:   ' getPassedFailed(pcInfo.Nc2DInfo.Violations.Extended)];
        Nc2DViolations{2} = ['Offset Check:     ' getPassedFailed(pcInfo.Nc2DInfo.Violations.Offset)];
        Nc2DViolations{3} = ['Inaccuracy Check: ' getPassedFailed(pcInfo.Nc2DInfo.Violations.Inaccurate)];
        for i = 1:maxNumViolations
            if endsWith(Nc2DViolations{i},'Failed')
                numNc2DViolations = numNc2DViolations + 1;
            end
        end
        
        % Extract peak Ncdot time and sigma
        TMeanRate  = pcInfo.Nc2DInfo.TMeanRate;
        TSigmaRate = pcInfo.Nc2DInfo.TSigmaRate;
        if isnan(TMeanRate) || isnan(TSigmaRate)
            TMeanRate  = pcInfo.Nc2DInfo.TQ0min;
            TSigmaRate = pcInfo.Nc2DInfo.SigmaQ0min;
        end

        % Calculate the conjunction time bounds using the Nc2D information
        gamma = 1e-6;    
        dtau = TSigmaRate * sqrt(2)*erfcinv(gamma);
        Ta_Nc2D = TMeanRate - dtau;
        Tb_Nc2D = TMeanRate + dtau;
        
        if params.verbose
            disp(params.verboseTextSep);
            disp(['Pc2D_Hall Pc = ' smart_exp_format(Nc2D,params.PcNsf,[false true])]);
            if ~isnan(Ta_Nc2D) && ~isnan(Tb_Nc2D)
                disp(['  Pc2D_Hall Conjunction Bounds = ' smart_exp_format(Ta_Nc2D,params.TcNsf,[false true]) ' to ' smart_exp_format(Tb_Nc2D,params.TcNsf,[false true])]);
            end
            disp('Pc2D_Hall Violations:');
            for i = 1:maxNumViolations
                disp(['  ' Nc2DViolations{i}]);
            end
        end
    end
    
    Nc3D = [];
    Tc3DHall = [];
    Rc3DHall = [];
    PcCum3DHall = [];
    Nc3DViolations = {};
    Ta_Nc3D = [];
    Tb_Nc3D = [];
    numNc3DViolations = 0;
    if params.plot_Pc3D_Hall
        % Get Pc and Violation text
        Nc3D = pcInfo.Nc3D;
        maxNumViolations = 3;
        Nc3DViolations = cell(maxNumViolations,1);
        Nc3DViolations{1} = ['Extended Check:   ' getPassedFailed(pcInfo.Nc3DInfo.Violations.Extended)];
        Nc3DViolations{2} = ['Offset Check:     ' getPassedFailed(pcInfo.Nc3DInfo.Violations.Offset)];
        Nc3DViolations{3} = ['Inaccurate Check: ' getPassedFailed(pcInfo.Nc3DInfo.Violations.Inaccurate)];
        for i = 1:maxNumViolations
            if endsWith(Nc3DViolations{i},'Failed')
                numNc3DViolations = numNc3DViolations + 1;
            end
        end
        
        % Get conjunction info
        Tc3DHall = pcInfo.Nc3DInfo.Teph;
        Pc3DHall = pcInfo.Nc3DInfo.Nc;
        Rc3DHall = pcInfo.Nc3DInfo.Ncdot;
        PcCum3DHall = pcInfo.Nc3DInfo.Nccum;
        Ta_Nc3D = pcInfo.Nc3DInfo.TaConj;
        Tb_Nc3D = pcInfo.Nc3DInfo.TbConj;
        
        % Reduce the data set just to the actual Rc data regions in case
        % Pc3D_Hall overruns
        firstIdx = find(Rc3DHall,1,'first');
        lastIdx = find(Rc3DHall,1,'last');
        if ~isempty(firstIdx)
            % Include one zero both before and after the actual data as
            % long as it doesn't go past the array bounds
            if firstIdx > 1; firstIdx = firstIdx - 1; end
            if lastIdx < numel(Rc3DHall); lastIdx = lastIdx + 1; end
            Tc3DHall = Tc3DHall(firstIdx:lastIdx);
            Rc3DHall = Rc3DHall(firstIdx:lastIdx);
            PcCum3DHall = PcCum3DHall(firstIdx:lastIdx);
        end
        
        % Get peak information
        [Rc3DHallPeak, iPeak] = max(Rc3DHall);
        if Rc3DHallPeak > 0
            Tc3DHallPeak = Tc3DHall(iPeak);
        end
        
        if params.verbose
            disp(params.verboseTextSep);
            disp(['Pc3D_Hall Pc = ' smart_exp_format(Nc3D,params.PcNsf,[false true])]);
            disp(['  NcRate peak = ' smart_exp_format(Rc3DHallPeak,params.RcNsf,[false true])]);
            if Rc3DHallPeak > 0
                tscale = Pc3DHall / Rc3DHallPeak;
                disp(['  NcRate peak time relative to TCA (s) = ' smart_exp_format(Tc3DHallPeak,params.TcNsf,[false true])]);
                disp(['  Effective Duration (s) = Nc/NcdotPeak = ' smart_exp_format(tscale,params.TcNsf,[false true])]);
                if ~isnan(Ta_Nc3D) && ~isnan(Tb_Nc3D)
                    disp(['  Pc3D_Hall Conjunction Bounds = ' smart_exp_format(Ta_Nc3D,params.TcNsf,[false true]) ' to ' smart_exp_format(Tb_Nc3D,params.TcNsf,[false true])]);
                end
            end
            disp('Pc3D_Hall Violations:');
            for i = 1:maxNumViolations
                disp(['  ' Nc3DViolations{i}]);
            end
        end
    end
    
    SDMCPc = [];
    SDMCPcUnc = [];
    SDMCNumHits = [];
    SDMCNumTrials = [];
    SDMCTrials = [];
    SDMCViolations = {};
    numSDMCViolations = 0;
    if params.plot_Pc_SDMC
        % Get Pc and Violation text
        maxNumViolations = 2;
        SDMCViolations = cell(maxNumViolations,1);
        SDMCViolations{1} = ['Extended Check:   ' getPassedFailed(pcInfo.SDMCInfo.Violations.Extended)];
        SDMCViolations{2} = ['Offset Check:     ' getPassedFailed(pcInfo.SDMCInfo.Violations.Offset)];
        for i = 1:maxNumViolations
            if endsWith(SDMCViolations{i},'Failed')
                numSDMCViolations = numSDMCViolations + 1;
            end
        end
        
        % Get conjunction info
        SDMCNumHits = pcInfo.SDMCInfo.numHits;
        SDMCNumTrials = pcInfo.SDMCInfo.numTrials;
        SDMCTrials = pcInfo.SDMCInfo.trialData;
        
        % Recalculate the Pc and uncertainty
        [SDMCPc, SDMCPcUnc] = binofit(SDMCNumHits, SDMCNumTrials, conf_alpha);
        if SDMCPc ~= pcInfo.SDMCPc
            warning('SDMC Pc mismatch detected');
        end
        if ~isequal(SDMCPcUnc,pcInfo.SDMCInfo.PcUnc)
            warning('SDMC Pc uncertainty mismatch detected, confidence level inputs are likely different');
        end
        
        % Reduce reported trial data just to actual hits and then get
        % min/max hit times
        if ~isempty(SDMCTrials)
            SDMCTrials = SDMCTrials(SDMCTrials.hitIndicator == 1,:);
            Nhit = height(SDMCTrials);
        else
            Nhit = 0;
        end
        if Nhit ~= SDMCNumHits
            warning('Mismatch detected between SDMC summary and trial data, SDMC graph outputs will be inaccurate!');
        end
        if Nhit > 0
            min_thit = min(SDMCTrials.hitTime);
            max_thit = max(SDMCTrials.hitTime);
            [ycdf, tcdf] = ecdf(SDMCTrials.hitTime);
            ncdf = ycdf * numel(SDMCTrials.hitTime);
            Pcdf = ncdf / SDMCNumTrials;
            [~, UUcdf] = binofit(round(ncdf),SDMCNumTrials);
        else
            min_thit = Inf;
            max_thit = -Inf;
            tcdf = [];
            Pcdf = 0;
            [~, UUcdf] = binofit(0,SDMCNumTrials);
        end
        
        if params.verbose
            disp(params.verboseTextSep);
            disp(['Pc_SDMC Pc = ' smart_exp_format(SDMCPc,params.PcNsf,[false true])]);
            disp(['  ' num2str(params.conf_level*100) '% Pc Conf = ' smart_exp_format(SDMCPcUnc(1),params.PcNsf,[false true]) ' to ' smart_exp_format(SDMCPcUnc(2),params.PcNsf,[false true])]);
            disp(['  Num Hits = ' smart_exp_format(SDMCNumHits)]);
            disp(['  Num Trials = ' smart_exp_format(SDMCNumTrials)]);
            if Nhit > 0
                disp(['  Min hit time relative to TCA (s) = ' smart_exp_format(min_thit,params.TcNsf,[false true])]);
                disp(['  Max hit time relative to TCA (s) = ' smart_exp_format(max_thit,params.TcNsf,[false true])]);
            else
                disp('  Cannot display min/max hit times since no hit data is available');
            end
            disp('Pc_SDMC Violations:');
            for i = 1:maxNumViolations
                disp(['  ' SDMCViolations{i}]);
            end
        end
    end
    
    %% Determine x-range values
    if params.plot_Pc2D_Hall || params.plot_Pc3D_Hall
        xMin = min([Ta_Nc2D Tc3DHall Ta_Nc3D]);
        xMax = max([Tb_Nc2D Tc3DHall Tb_Nc3D]);
        xRangeDataFound = true;
        if isnan(xMin) || isnan(xMax)
            xRangeDataFound = false;
        end
    else
        xMin = Inf;
        xMax = -Inf;
        xRangeDataFound = false;
    end
    
    %% Bin the SDMC data for the rate plots
    tfc_bins.N = 0;
    tfc_bins.x1 = nan;
    tfc_bins.x2 = nan;
    if params.plot_Pc_SDMC
        tfc_bins.N = 50;
        if xRangeDataFound
            tfc_bins.x1 = xMin;
            tfc_bins.x2 = xMax;
            if Nhit > 0
                tfc_bins.x1 = min([xMin min_thit]);
                tfc_bins.x2 = max([xMax max_thit]);
            end
        else
            % Add 1 second wrapper around cases where the number of hits is 0
            % or 1, otherwise use the min/max hit time values
            if Nhit == 0
                tfc_bins.x1 = -0.5;
                tfc_bins.x2 = 0.5;
            elseif Nhit == 1
                tfc_bins.x1 = min_thit-0.5;
                tfc_bins.x2 = max_thit+0.5;
            else
                tfc_bins.x1 = min_thit;
                tfc_bins.x2 = max_thit;
            end
        end
        
        if isnan(tfc_bins.x1) || isnan(tfc_bins.x2) || ...
                isinf(tfc_bins.x1) || isinf(tfc_bins.x2)
            error('Bad Monte Carlo time bin limits generated');
        end
        
        % Setup N evenly spaced bins between the x1 and x2 values
        tfc_bins.wid = tfc_bins.x2 - tfc_bins.x1;
        tfc_bins.del = tfc_bins.wid / tfc_bins.N;
        tfc_bins.xhi = tfc_bins.x1 + (1:tfc_bins.N)*tfc_bins.del;
        tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
        tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;
        
        % Get overall x min and max values
        xmn0 = min(min(tfc_bins.xlo),xMin);
        xmx0 = max(max(tfc_bins.xhi),xMax);
        
        % Bin the hits
        if isempty(SDMCTrials) || height(SDMCTrials) == 0
            nbin = zeros(size(tfc_bins.xmd));
        else
            binranges = [tfc_bins.xlo tfc_bins.xhi(end)];
            nbin = histcounts(SDMCTrials.hitTime,binranges);
            % Get the nbin into the correct size format
            if ~isequal(size(nbin),size(tfc_bins.xmd))
                nbin = nbin';
                if ~isequal(size(nbin),size(tfc_bins.xmd))
                    error('Nbin size mismatch');
                end
            end
        end
        
        % Create rate confidence intervals for each bin
        rbin = zeros(size(nbin));
        rtop = rbin;
        rbot = rbin;
        nsum = 0;
        Nsamp = pcInfo.SDMCInfo.numTrials;
        for n = 1:tfc_bins.N
            % Rate curve and confidence interval
            [a,b] = binofit(nbin(n),Nsamp,conf_alpha);
            rbin(n) = a / tfc_bins.del;
            rtop(n) = b(2) / tfc_bins.del;
            rbot(n) = b(1) / tfc_bins.del;
            nsum = nsum + nbin(n);
        end
    else
        % Set defaults for the rate plot if SDMC isn't being plotted
        if xRangeDataFound
            xmn0 = xMin;
            xmx0 = xMax;
        else
            xmn0 = -0.5;
            xmx0 = 0.5;
        end
        tfc_bins.xmd = [];
        rbin = [];
        rtop = [];
        rbot = [];
    end
    
    %% Get final x and y range values
    xrng = plot_range([xmn0 xmx0],params.fig.xpad);
    
    ycmb = 0;
    if params.plot_PcCircle; ycmb = [ycmb Pc2D]; end
    if params.plot_Pc2D_Hall; ycmb = [ycmb Nc2D]; end
    if params.plot_Pc3D_Hall; ycmb = [ycmb Nc3D]; end
    if params.plot_Pc_SDMC; ycmb = [ycmb UUcdf(end,2)]; end
    yrng = plot_range(ycmb,params.fig.ypad);
    
    ycmb2 = 0;
    if params.plot_Pc3D_Hall; ycmb2 = [ycmb2 Rc3DHallPeak]; end
    if params.plot_Pc_SDMC; ycmb2 = [ycmb2 max(rtop)]; end
    yrng2 = plot_range(ycmb2,params.fig.ypad);
    
    %% Setup figure and figure parameters
    fh = figure_set_up(12/8,'on');
    
    %% Fix the SDMC data to span the entire x-range
    if params.plot_Pc_SDMC
        if isempty(tcdf)
            tcdf = xmn0;
        end
        [xVals,yVals] = FitToRange([xmn0 xmx0],tcdf',Pcdf');
        [~,yLoVals] = FitToRange([xmn0 xmx0],tcdf',UUcdf(:,1)');
        [~,yHiVals] = FitToRange([xmn0 xmx0],tcdf',UUcdf(:,2)');
        [xSDMC, ySDMC] = stairs(xVals', yVals');
        [~, yLoSDMC] = stairs(xVals', yLoVals');
        [~, yHiSDMC] = stairs(xVals', yHiVals');
    else
        xSDMC = [];
        ySDMC = [];
        yLoSDMC = [];
        yHiSDMC = [];
    end
    
    %% Plot the cumulative Pc
    subplot(2,2,1);
    hold on;
    PlotErrorBand(xSDMC', yLoSDMC', yHiSDMC', params);
    PlotFormattedLine('Pc_SDMC', xSDMC, ySDMC, params);
    PlotFormattedLine('PcCircle', [xmn0 xmx0], Pc2D, params);
    PlotFormattedLine('Pc2D_Hall', [xmn0 xmx0], Nc2D, params);
    PlotFormattedLine('Pc3D_Hall', Tc3DHall, PcCum3DHall, params);
    xlim(xrng);
    ylim(yrng);
    SetLabelsAndWeights('Time from TCA (s)','P_{c} Cumulative',params);
    
    %% Plot the Pc rate (but only if there is data to plot)
    pcRatePlotted = false;
    if numel(ycmb2) > 1
        subplot(2,2,3);
        hold on;
        PlotErrorBand(tfc_bins.xmd, rbot, rtop, params);
        PlotFormattedLine('Pc_SDMC', tfc_bins.xmd, rbin, params);
        PlotFormattedLine('Pc3D_Hall', Tc3DHall, Rc3DHall, params);
        PlotFormattedLine('Pc2D_Hall_ConjBound',[Ta_Nc2D Ta_Nc2D],[min(ycmb2) max(ycmb2)], params);
        PlotFormattedLine('Pc2D_Hall_ConjBound',[Tb_Nc2D Tb_Nc2D],[min(ycmb2) max(ycmb2)], params);
        PlotFormattedLine('Pc3D_Hall_ConjBound',[Ta_Nc3D Ta_Nc3D],[min(ycmb2) max(ycmb2)], params);
        PlotFormattedLine('Pc3D_Hall_ConjBound',[Tb_Nc3D Tb_Nc3D],[min(ycmb2) max(ycmb2)], params);
        xlim(xrng);
        ylim(yrng2);
        SetLabelsAndWeights('Time from TCA (s)','P_{c} Rate (s^{-1})',params);
        pcRatePlotted = true;
    end
    
    %% Add title, legend, and notes on the right side of the plot
    subplot(2,2,2);
    % Title
    axis off;
    hold on;
    nt=0;
    nt=nt+1; titl{nt} = strrep(conjID,'_','\_');
    if pcInfo.covXcorr_corrections_applied
        titl{nt} = [titl{nt} ' (w/ CCC)'];
    end
    nt=nt+1; titl{nt} = 'Conjunction Collision Probability Temporal Plot';
    nt=nt+1; titl{nt} = strrep(conjInfoStr,'deg','\circ');
    nt=nt+1; titl{nt} = conjInfoStr2;
    text(0.45,1.1,titl,'FontSize',params.fig.tfsz,'FontWeight',params.fig.tfwt, ...
        'HorizontalAlignment','center','VerticalAlignment','top');
    
    % Legend
    Nlgnd = 10;
    hlgnd = zeros(1,Nlgnd);
    slgnd = cell(1,Nlgnd);
    n = 0;
    [hlgnd, slgnd, n] = AddLineToLegend('PcCircle',hlgnd,slgnd,n,params,numPc2DViolations, ...
        [' = ' smart_exp_format(Pc2D,params.PcNsf,[false true])]);
    [hlgnd, slgnd, n] = AddLineToLegend('Pc2D_Hall',hlgnd,slgnd,n,params,numNc2DViolations, ...
        [' = ' smart_exp_format(Nc2D,params.PcNsf,[false true])]);
    [hlgnd, slgnd, n] = AddLineToLegend('Pc3D_Hall',hlgnd,slgnd,n,params,numNc3DViolations, ...
        [' = ' smart_exp_format(Nc3D,params.PcNsf,[false true])]);
    [hlgnd, slgnd, n] = AddLineToLegend('Pc_SDMC',hlgnd,slgnd,n,params,numSDMCViolations, ...
        [' = ' smart_exp_format(SDMCPc,params.PcNsf,[false true]) ...
        ' (' smart_exp_format(SDMCNumHits,10) ' hits / ' ...
        smart_exp_format(SDMCNumTrials,10) ' trials)']);
    if isfield(pcInfo,'SDMCInfo') && pcInfo.SDMCInfo.BFMCsubstitution
        slgnd{n} = strrep(slgnd{n},'SDMC','BFMC');
    end
    [hlgnd, slgnd, n] = AddErrorBandToLegend(hlgnd,slgnd,n,params,SDMCPc,SDMCPcUnc,numSDMCViolations);
    if pcRatePlotted
        [hlgnd, slgnd, n] = AddLineToLegend('Pc2D_Hall_ConjBound',hlgnd,slgnd,n,params);
        [hlgnd, slgnd, n] = AddLineToLegend('Pc3D_Hall_ConjBound',hlgnd,slgnd,n,params);
    end
    axpos = get(gca,'OuterPosition');
    lft = axpos(1);
    bot = axpos(2);
    wid = axpos(3);
    hgt = axpos(4);
    lgpos(1) = lft+0.05*wid;
    lgpos(2) = bot+0.20*hgt;
    lgpos(3) = 0.90*wid;
    lgpos(4) = 0.50*hgt;
    legend(hlgnd(1:n),slgnd(1:n),'Position',lgpos, ...
        'FontSize',params.fig.lfsz,'FontWeight',params.fig.lfwt);
    legend('boxoff');
    
    % Violation checks
    subplot(2,2,4);
    axis off;
    hold on;
    txt = {};
    nt = 0;
    nt = nt + 1; txt{nt} = 'Usage Violation Checks';
    [txt, nt] = DisplayViolations(params.plot_PcCircle, '2D-Pc Method', Pc2DViolations, txt, nt);
    [txt, nt] = DisplayViolations(params.plot_Pc2D_Hall, '2D-Nc Method', Nc2DViolations, txt, nt);
    [txt, nt] = DisplayViolations(params.plot_Pc3D_Hall, '3D-Nc Method', Nc3DViolations, txt, nt);
    if isfield(pcInfo,'SDMCInfo') && pcInfo.SDMCInfo.BFMCsubstitution
        [txt, ~]  = DisplayViolations(params.plot_Pc_SDMC, 'BFMC Method', SDMCViolations, txt, nt);
    else
        [txt, ~]  = DisplayViolations(params.plot_Pc_SDMC, 'SDMC Method', SDMCViolations, txt, nt);
    end
    text(-0.25,1.4,txt,'FontSize',params.fig.tfsz,'FontWeight',params.fig.tfwt,'VerticalAlignment','top');
    
    % Conjunction bounds
    txt = {};
    nt = 0;
    if params.plot_Pc2D_Hall || params.plot_Pc3D_Hall || params.plot_Pc_SDMC
        % nt = nt + 1; txt{nt} = 'Conjunction Bounds';
        nt = nt + 1; txt{nt} = 'Conjunction Duration Bounds';
        nt = nt + 1; txt{nt} = ' ';
        if params.plot_Pc2D_Hall
            nt = nt + 1; txt{nt} = ['    2D-Nc = ' smart_exp_format(Ta_Nc2D,params.TcNsf,[false true]) ' to ' smart_exp_format(Tb_Nc2D,params.TcNsf,[false true]) ' sec'];
        end
        if params.plot_Pc3D_Hall
            if ~isnan(Ta_Nc3D) && ~isnan(Tb_Nc3D)
                nt = nt + 1; txt{nt} = ['    3D-Nc = ' smart_exp_format(Ta_Nc3D,params.TcNsf,[false true]) ' to ' smart_exp_format(Tb_Nc3D,params.TcNsf,[false true]) ' sec'];
            end
            if Rc3DHallPeak > 0
                nt = nt + 1; txt{nt} = ['    Peak Pc rate = ' smart_exp_format(Tc3DHallPeak,params.TcNsf,[false true]) ' sec'];
            end
        end   
        if params.plot_Pc_SDMC
            if Nhit > 0
                nt = nt + 1; txt{nt} = ['    SDMC hits = ' smart_exp_format(min_thit,params.TcNsf,[false true]) ' to ' smart_exp_format(max_thit,params.TcNsf,[false true]) ' sec'];
            end
        end
        nt = nt + 1; txt{nt} = ' ';
        nt = nt + 1; txt{nt} = '    (All times relative to TCA)';
    end
    text(0.45,1.4,txt,'FontSize',params.fig.tfsz,'FontWeight',params.fig.tfwt,'VerticalAlignment','top');
    
    %% Save the figure and close it if the figure has been saved
    [figSaved, fileName] = SaveConjFigure(fh, conjID, 'PcTime', params);
    if figSaved
        close(fh);
        figInfo = fileName;
    else
        figInfo = fh;
    end
end

%% Displays violations with colors
function [txt, nt] = DisplayViolations(plotEnabled, dispText, violations, txt, nt)
    if plotEnabled
        nt = nt + 1; txt{nt} = ' ';
        nt = nt + 1; txt{nt} = ['\color{black}    ' dispText ' Violations'];
        for i = 1:length(violations)
            nt = nt + 1;
            violations{i} = regexprep(violations{i},' +',' ');
            if endsWith(violations{i},'Failed')
                txt{nt} = ['\color{red}        ' violations{i}];
            else
                txt{nt} = ['\color{black}        ' violations{i}];
            end
        end
    end
end

%% Fits x and y values to the range determined by xrng
function [x, y] = FitToRange(xrng, x, y)
    if xrng(1) < x(1)
        x = [xrng(1) x];
        y = [y(1)    y];
    end
    if xrng(end) > x(end)
        x = [x xrng(end)];
        y = [y y(end)];
    end
end

%% Adds the error band to the legend
function [hlgnd, slgnd, n] = AddErrorBandToLegend(hlgnd,slgnd,n,params,PcMC,PcMC_unc,numViolations)
    if params.plot_Pc_SDMC
        n = n + 1;
        [~,xstr1,xstr2] = smart_error_range(PcMC,PcMC_unc(1),PcMC_unc(2));
        slgnd{n} = ['MC ' smart_exp_format(params.conf_level*100,10) '% Confidence ' ...
            xstr1 ' \leq Pc \leq ' xstr2];
        if numViolations > 0
            slgnd{n} = ['\color{red}' slgnd{n} '\color{black}'];
        end
        xplt = nan(1,2); yplt = xplt;
        xband = [xplt flip(xplt,2)];
        yband = [yplt flip(yplt,2)];
        hlgnd(n) = fill(xband,yband,params.fig.col_MC_err,'EdgeColor',params.fig.col_MC_err);
    end
end

%% Plots the error band
function PlotErrorBand(xvals, ybot, ytop, params)
    if params.plot_Pc_SDMC
        xband = [xvals flip(xvals,2)];
        yband = [ybot, flip(ytop,2)];
        fill(xband,yband,params.fig.col_MC_err,'EdgeColor',params.fig.col_MC_err);
    end
end

%% Plots a formatted line
function PlotFormattedLine(plotType,xVals,yVals,params)
    [pltEnabled, mrkr, mcol, msiz, lcol, lnst, lnwd] = GetPlotFormat(plotType,params);
    if pltEnabled
        if isequal(size(xVals),[1 2]) && isequal(size(yVals),[1 1])
            yVals = yVals + zeros(1,2);
        end
        plot(xVals,yVals,...
        'LineStyle',lnst,...
        'LineWidth',lnwd,...
        'Color',lcol,...
        'Marker',mrkr,...
        'MarkerFaceColor',mcol,...
        'MarkerEdgeColor',mcol,...
        'MarkerSize',msiz);
    end
end

%% Adds legend data
function [hlgnd, slgnd, n] = AddLineToLegend(plotType,hlgnd,slgnd,n,params,numViolations,extraText)
    if nargin == 5
        numViolations = 0;
        extraText = '';
    elseif nargin == 6
        extraText = '';
    end
    [pltEnabled, mrkr, mcol, msiz, lcol, lnst, lnwd, name] = GetPlotFormat(plotType,params);
    if pltEnabled
        n = n + 1;
        slgnd{n} = [name extraText];
        if numViolations > 0
            slgnd{n} = ['\color{red}' slgnd{n} '\color{black}'];
        end
        xplt = nan(2,1); yplt = xplt;
        hlgnd(n) = plot(xplt,yplt, ...
            'LineStyle',lnst, ...
            'LineWidth',lnwd, ...
            'Color',lcol, ...
            'Marker',mrkr, ...
            'MarkerFaceColor',mcol, ...
            'MarkerEdgeColor',mcol, ...
            'MarkerSize',msiz);
    end
end

%% Returns passed or failed text based on boolean value
function [passFailText] = getPassedFailed(value)
    if ~value
        passFailText = 'Passed';
    else
        passFailText = 'Failed';
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
% L. Baars       | 2023-Mar-24 | Initial Development
% L. Baars       | 2023-Apr-20 | Added SDMC offset violation.
% L. Baars       | 2023-May-19 | Added data quality and invalid 6x6 checks.
% L. Baars       | 2023-Jul-25 | Added figInfo output
% D. Hall        | 2023-Aug-30 | Changed errors to warnings when plots
%                                cannot be made. Added plot_conjID_string
%                                input parameter. Updated text displayed in
%                                plot outputs.
% L. Baars       | 2024-Jan-11 | Changed references to PcConjPlaneCircle
%                                into PcCircle. Also, tool no longer
%                                returns early if a plot cannot be made.
% D. Hall        | 2024-Sep-16 | Added support for plotting BFMC data.
% E. Toumey      | 2025-Mar-02 | Moved file for new directory structure.
% L. Baars       | 2025-Aug-29 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
