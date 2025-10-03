function [figInfo] = CA_Dist_Plot(params, pcInfo)
% CA_Dist_Plot - Produces close approach distribution plots of several Pc
%                theories against each other.
%
% Syntax: [figInfo] = CA_Dist_Plot(params, pcInfo);
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
%   This function creates close approach distribution plots of several Pc
%   theories against each other. The outputs that are generated include
%   three plots and some informational text.
%
%   The first plot (top left corner) graphs the CA distribution using a
%   square axis centered on the primary object. The min/max x and y axis
%   values are based on the largest absolute value of any single x or y
%   value that has been plotted. This is the "unzoomed" plot of the
%   conjunction.
%
%   The second plot (across the bottom) graphs the CA distribution without
%   making the x and y axis equal. The graph displays all of the data
%   points using a smallest possible bounding box. This is the "zoomed best
%   fit" plot of the conjunction.
%
%   The third plot (top right corner) graphs the CA distribution zoomed in
%   to 2x the combined HBR. This is the "zoomed HBR" plot of the
%   conjunction.
%
%   The informational text in the upper middle includes general conjunction
%   information, calculated Pc values, and a summary of Pc usage violation
%   checks. If a usage violation check fails for a particular Pc
%   calculation theory, both the violation check and the reported Pc value
%   will be displayed in a red font.
%
%   Data is only plotted if either or both the PcCircle or Pc_SDMC data is
%   selected for display. Other Pc options provide data to display in the
%   informational text, but do not plot any data. Available Pc theories
%   include:
%     PcCircle  - labeled as "2D-Pc Method"
%     Pc2D_Hall - labeled as "2D-Nc Method"
%     Pc3D_Hall - labeled as "3D-Nc Method"
%     Pc_SDMC   - labeled as "SDMC Method" and "SDMC 95% confidence"
%
%   Note that a confidence interval is included for Pc_SDMC since those
%   values can be determined from the results of Monte Carlo processing.
%
%   It is expected that all Pc calculations have occurred before calling
%   this function and these calculations have been used to populate the
%   pcInfo structure with the needed data.
%
%   All plots can have the following elements displayed (depending on the
%   zoom level of the plot and actual SDMC results):
%
%     - Magenta ellipse - Indicates the 3-sigma combined covariance
%                         projected onto the nominal conjunction plane
%                         (only provided when "2D Pc" plots are enabled)
%
%     - Blue dots - Indicates the location of individual SDMC trials with a
%                   PCA distance greater than the combined HBR
%
%     - Red dots - Indicates the location of individual SDMC trials with a
%                  PCA distance less than or equal to the combined HBR
%
%     - Green circle - Indicates the size of the combined HBR centered on
%                      the primary object
%
%     - Black diamond - Indicates the nominal (unperterbed) miss distance
%                       of the secondary object at PCA. By default, this
%                       will be on the positive x-axis.
%
%     - Dashed lines - Indicates the "0" values of the x and y axis, their
%                      crossing indicates the position of the primary
%                      object in the default view.
%
%   The nominal conjunction plane definition used in all three graphs is
%   defined as follows:
%
%     - z-axis is along the negative of the relative velocity (from the
%       primary to the secondary object). In other words, the relative
%       velocity is going down into the paper while the z-axis is pointing
%       straight up from the paper.
%
%     - y-axis is the cross product of the relative postion (from the
%       primary to the secondary) and the relative velocity vector.
%
%     - x-axis is the cross product of the y-axis and the z-axis, making a
%       right handed coordinate system triad. Thus placing the nominal
%       secondary position on the positive x-axis at the nominal miss
%       distance.
%
%   An alternative conjunction plane is defined when the params.alt_ca_dist
%   is set to 1 or 2. If 1, the combined covariance is placed at
%   the origin. The HBR and nominal miss distance are placed on the
%   positive x-axis. The z-axis is the opposite of the nominal conjunction
%   plane, and the y-axis makes up the right handed coordinate system
%   triad. This conjunction plane is a 180 degree rotation about the y-axis
%   of the nominal conjunction plane and is translated in the positive
%   x-direction by the nominal miss distance. If 2, then the conj. plane is
%   rotated to align the axes with with principal axes of the 2D-Pc
%   ellipse.
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
%            default_params_CA_Dist_Plot.m function. See that function for
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
%       C1  - Primary object's ECI covariance matrix (meter units)    [6x6]
%
%       r2  - Secondary object's ECI position vector (m)       [3x1 or 1x3]
%
%       v2  - Secondary object's ECI velocity vector (m/s)     [3x1 or 1x3]
%
%       C2  - Secondary object's ECI covariance matrix (meter units)  [6x6]
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
    persistent pathsAdded researchCodeAvailable
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        addpath(fullfile(p,'Utils'));
        addpath(fullfile(p,'../LoggingAndStringReporting'));
        try
            s = what(fullfile(p, '../../../ResearchCode/OCMDB_PcMultiStep/src')); addpath(s.path);
            s = what(fullfile(p, '../../../ResearchCode/Utils/Plotting')); addpath(s.path);
            researchCodeAvailable = true;
        catch
            researchCodeAvailable = false;
        end
        pathsAdded = true;
    end
    
    %% Check Input Parameters
    figInfo = [];
    % Load the default parameters
    params = default_params_CA_Dist_Plot(params);
    % Check if a time plot needs to be generated
    if ~params.generate_ca_dist_plot
        return;
    end
    if ~researchCodeAvailable && params.AuxCAdistContour ~= 0
        error('The selected params.AuxCAdistContour value is not currently supported!');
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
        error('Cannot plot CA_Dist_Plot since basic conjunction parameters are missing');
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
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'Pc2DInfo') && allFieldsExist;
        if isfield(pcInfo,'Pc2DInfo')
            allFieldsExist = dispMissingField(pcInfo.Pc2DInfo, 'pcInfo.Pc2DInfo', 'EigV1') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Pc2DInfo, 'pcInfo.Pc2DInfo', 'EigL1') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.Pc2DInfo, 'pcInfo.Pc2DInfo', 'EigL2') && allFieldsExist;
        end
        if ~allFieldsExist
            % error('Cannot plot PcCircle data due to missing input data!');
            warning('Cannot plot PcCircle data due to missing input data; not making CA dist. plot');
            params.plot_PcCircle = false;
        end
    end
    
    % Check required parameters for displaying Pc2D_Hall data
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
        end
        if ~allFieldsExist
            if pcInfo.DataQualityError.invalidCov6x6
                warning('Invalid 6x6 covariance caused missing Pc2D_Hall parameters, turning off Pc2D_Hall plot');
                params.plot_Pc2D_Hall = false;
            else
                % error('Cannot display Pc2D_Hall data due to missing input data!');
                warning('Cannot display Pc2D_Hall data due to missing input data');
                params.plot_Pc2D_Hall = false;
            end
        end
    end
    
    % Check required parameters for displaying Pc3D_Hall data
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
        end
        if ~allFieldsExist
            if pcInfo.DataQualityError.invalidCov6x6
                warning('Invalid 6x6 covariance caused missing Pc3D_Hall parameters, turning off Pc3D_Hall plot');
                params.plot_Pc3D_Hall = false;
            else
                % error('Cannot display Pc3D_Hall data due to missing input data!');
                warning('Cannot display Pc3D_Hall data due to missing input data; not making CA dist. plot');
                params.plot_Pc3D_Hall = false;
            end
        end
    end
    
    % Check required parameters for plotting Pc_SDMC data
    if ~isfield(params,'plot_Pc_SDMC'); params.plot_Pc_SDMC = []; end
    if isempty(params.plot_Pc_SDMC); params.plot_Pc_SDMC = false; end
    if params.plot_Pc_SDMC
        allFieldsExist = true;
        allFieldsExist = dispMissingField(pcInfo, 'pcInfo', 'SDMCInfo') && allFieldsExist;
        if isfield(pcInfo,'SDMCInfo')
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'Violations') && allFieldsExist;
            if isfield(pcInfo.SDMCInfo,'Violations')
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.Violations, 'pcInfo.SDMCInfo.Violations', 'Extended') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.Violations, 'pcInfo.SDMCInfo.Violations', 'Offset') && allFieldsExist;
            end
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'numHits') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'numTrials') && allFieldsExist;
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'trialData') && allFieldsExist;
            if isfield(pcInfo.SDMCInfo,'trialData') && ~isempty(pcInfo.SDMCInfo.trialData)
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'hitIndicator') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'priPosX_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'priPosY_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'priPosZ_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'priVelX_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'priVelY_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'priVelZ_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'secPosX_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'secPosY_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'secPosZ_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'secVelX_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'secVelY_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.trialData, 'pcInfo.SDMCInfo.trialData', 'secVelZ_kmps') && allFieldsExist;
            end
            allFieldsExist = dispMissingField(pcInfo.SDMCInfo, 'pcInfo.SDMCInfo', 'overallPicTrialData') && allFieldsExist;
            if isfield(pcInfo.SDMCInfo,'overallPicTrialData') && ~isempty(pcInfo.SDMCInfo.overallPicTrialData)
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'hitIndicator') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'priPosX_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'priPosY_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'priPosZ_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'priVelX_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'priVelY_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'priVelZ_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'secPosX_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'secPosY_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'secPosZ_km') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'secVelX_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'secVelY_kmps') && allFieldsExist;
                allFieldsExist = dispMissingField(pcInfo.SDMCInfo.overallPicTrialData, 'pcInfo.SDMCInfo.overallPicTrialData', 'secVelZ_kmps') && allFieldsExist;
            end
        end
        if ~allFieldsExist
            if pcInfo.DataQualityError.invalidCov6x6
                warning('Invalid 6x6 covariance caused missing Pc_SDMC parameters, turning off Pc_SDMC plot');
                params.plot_Pc_SDMC = false;
            else
                % error('Cannot plot Pc_SDMC data due to missing input data!');
                % warning('Cannot plot Pc_SDMC data due to missing input data; not making CA dist. plot');
                params.plot_Pc_SDMC = false;
            end
        end
    end
    
    % Check to make sure at least one individual plot parameter has been
    % turned on
    if ~params.plot_PcCircle && ~params.plot_Pc_SDMC
        warning('No plots generated since neither the plot_PcCircle nor the plot_Pc_SDMC parameters were turned on!');
        return;
    end
    
    % General conjunction parameters
    [conjID, conjInfoStr, conjInfoStr2] = GetConjInfo(params, pcInfo);
    if ~isempty(params.plot_conjID_string)
        conjID = params.plot_conjID_string;
    end
    
    conjInfo = strsplit(conjInfoStr,' RelV');
    conjInfo{2} = ['RelV' conjInfo{2}];
    if ~isempty(conjInfoStr2)
        conjInfo2 = strsplit(conjInfoStr2,' Sec');
        conjInfo2{2} = ['Sec' conjInfo2{2}];
    else
        conjInfo2{1} = '';
        conjInfo2{2} = '';
    end
    
    Pc2D = [];
    if params.plot_PcCircle
        Pc2D = pcInfo.Pc2D;
    end
    Nc2D = [];
    if params.plot_Pc2D_Hall
        Nc2D = pcInfo.Nc2D;
    end
    Nc3D = [];
    if params.plot_Pc3D_Hall
        Nc3D = pcInfo.Nc3D;
    end
    SDMCPc = [];
    SDMCPcUnc = [];
    trialData = [];
    overallPicTrialData = [];
    if params.plot_Pc_SDMC
        conf_alpha = 1 - params.conf_level;
        SDMCNumHits = pcInfo.SDMCInfo.numHits;
        SDMCNumTrials = pcInfo.SDMCInfo.numTrials;
        [SDMCPc, SDMCPcUnc] = binofit(SDMCNumHits, SDMCNumTrials, conf_alpha);
        trialData = pcInfo.SDMCInfo.trialData;
        overallPicTrialData = pcInfo.SDMCInfo.overallPicTrialData;
    end
    
    %% Plot the data
    fh = figure_set_up(1.5,'on');
    PlotTrialData(params, pcInfo, overallPicTrialData, 1);
    PlotTrialData(params, pcInfo, trialData, 2);
    PlotTrialData(params, pcInfo, overallPicTrialData, 3);
    
    %% Add title and notes
    subplot(2,3,2);
    axis off;
    hold on;
    
    % Title
    titl = {};
    nt=0;
    nt=nt+1; titl{nt} = strrep(conjID,'_','\_');
    if pcInfo.covXcorr_corrections_applied
        titl{nt} = [titl{nt} ' (w/ CCC)'];
    end
    nt=nt+1; titl{nt} = 'Conjunction Close Approach Distribution Plot';
    nt=nt+1; titl{nt} = conjInfo{1};
    nt=nt+1; titl{nt} = strrep(conjInfo{2},'deg','\circ');
    nt=nt+1; titl{nt} = conjInfo2{1};
    nt=nt+1; titl{nt} = conjInfo2{2};
    text(0.45,1.15,titl,'FontSize',params.fig.tfsz,'FontWeight',params.fig.tfwt,...
        'HorizontalAlignment','center','VerticalAlignment','top');
    
    txt = {};
    txt2 = {};
    nt = 0;
    nt2 = 0;
    [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('PcCircle',params);
    if pltEnabled
        [nv, violText] = GetNumberOfViolations(1, pcInfo, name);
        nt=nt+1; txt{nt} = [name ' = ' smart_exp_format(Pc2D,params.PcNsf,[false true])];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        nt2=nt2+1; txt2{nt2} = violText;
    end
    [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('Pc2D_Hall',params);
    if pltEnabled
        [nv, violText] = GetNumberOfViolations(2, pcInfo, name);
        nt=nt+1; txt{nt} = [name ' = ' smart_exp_format(Nc2D,params.PcNsf,[false true])];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        nt2=nt2+1; txt2{nt2} = violText;
    end
    [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('Pc3D_Hall',params);
    if pltEnabled
        [nv, violText] = GetNumberOfViolations(3, pcInfo, name);
        nt=nt+1; txt{nt} = [name ' = ' smart_exp_format(Nc3D,params.PcNsf,[false true])];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        nt2=nt2+1; txt2{nt2} = violText;
    end
    [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('Pc_SDMC',params);
    if pltEnabled
        [nv, violText] = GetNumberOfViolations(4, pcInfo, name);
        if pcInfo.SDMCInfo.BFMCsubstitution
            name = strrep(name,'SDMC','BFMC');
        end
        [nameparts,~] = string_parts(name);
        nt=nt+1; txt{nt} = [name ' = ' ...
            smart_exp_format(SDMCPc,params.PcNsf,[false true])];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        nt=nt+1; txt{nt} = ['  (' smart_exp_format(SDMCNumHits,10) ' hits / ' ...
            smart_exp_format(SDMCNumTrials,10) ' trials)'];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        [~,xstr1,xstr2] = smart_error_range(SDMCPc,SDMCPcUnc(1),SDMCPcUnc(2));
        nt=nt+1; txt{nt} = [nameparts{1} ' ' smart_exp_format(params.conf_level*100,10) '% Confidence:'];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        nt=nt+1; txt{nt} = ['  ' xstr1 ' \leq Pc \leq ' xstr2];
        if nv > 0; txt{nt} = ['\color{red}' txt{nt} '\color{black}']; end
        nt2=nt2+1; txt2{nt2} = violText;
    end
    text(0,0.6,txt,'FontSize',params.fig.tfsz,'FontWeight',params.fig.tfwt,'VerticalAlignment','top');
    text(-0.1,0.0,txt2,'FontSize',params.fig.tfsz,'FontWeight',params.fig.tfwt,'VerticalAlignment','top');
    
    %% Save the figure and close it if the figure has been saved
    altTxt = '';
    if params.alt_ca_dist ~= 0
        altTxt = ['Alt' num2str(params.alt_ca_dist) ];
    end
    auxTxt = '';
    if params.AuxCAdistContour ~= 0
        auxTxt = ['Aux' num2str(params.AuxCAdistContour) ];
    end
    [figSaved, fileName] = SaveConjFigure(fh, conjID, ['CADist' altTxt auxTxt], params);
    if figSaved
        close(fh);
        figInfo = fileName;
    else
        figInfo = fh;
    end

end

% Plots actual trial data
function PlotTrialData(params, pcInfo, trialData, plotType)
    %% Calculate the coordinate frame
    X1 = [pcInfo.r1 pcInfo.v1]'; X2 = [pcInfo.r2 pcInfo.v2]';
    [~,X1CA,X2CA] = FindNearbyCA(X1,X2);
    r1CA = X1CA(1:3)'; v1CA = X1CA(4:6)';
    r2CA = X2CA(1:3)'; v2CA = X2CA(4:6)';
    [Xca0, Yca0, phat, pnrm] = CalcCloseApproachCoordinates(r1CA, v1CA, r2CA, v2CA, params);
    if params.alt_ca_dist == 1
        Xca0 = Xca0 + pnrm;
    end
    
    if ~isempty(trialData)
        numTrials = height(trialData);
        hitInd = trialData.hitIndicator == 1;
    else
        numTrials = 0;
        hitInd = [];
    end
    Xca = nan(numTrials,1);
    Yca = nan(numTrials,1);
    for i = 1:numTrials
        % Get sampling trial data
        r1 = [trialData.priPosX_km(i)   trialData.priPosY_km(i)   trialData.priPosZ_km(i)];
        v1 = [trialData.priVelX_kmps(i) trialData.priVelY_kmps(i) trialData.priVelZ_kmps(i)];
        r2 = [trialData.secPosX_km(i)   trialData.secPosY_km(i)   trialData.secPosZ_km(i)];
        v2 = [trialData.secVelX_kmps(i) trialData.secVelY_kmps(i) trialData.secVelZ_kmps(i)];
        % Relative state for the sampling trial
        [Xca(i), Yca(i)] = CalcCloseApproachCoordinates(r1, v1, r2, v2, params, phat);
    end
    % Convert values from km to m
    Xca = Xca * 1e3;
    Yca = Yca * 1e3;
    if params.alt_ca_dist == 1
        Xca = Xca + pnrm;
    end
    
    %% Calculate the ellipse depicting the 2D-Pc projected covariance
    xellipse = [];
    yellipse = [];
    angleCAdist = 0;
    if params.plot_PcCircle
        
        % Process the primary and secondary ellipse information
        if plotType == 3 && isfield(pcInfo.Pc2DInfo,'EigV1Pri')
            % Calculate the primary and secondary covariance ellipses for
            % the unequally zoomed-in plot (plotType = 3)
            PlotPriSecEllipses = true;
            % Get the semi-major eigenvector and eigenvalues from the projected
            % covariance
            eigV1 = pcInfo.Pc2DInfo.EigV1Pri;
            eigL1 = pcInfo.Pc2DInfo.EigL1Pri;
            eigL2 = pcInfo.Pc2DInfo.EigL2Pri;
            % Get the semi-major and semi-minor axis 3-sigma values
            a = sqrt(eigL1)*3;
            b = sqrt(eigL2)*3;
            % Calculate the clock angle in the conjunction plane coordinates
            % (includes a rotation from the conjunction plane used by 2D Pc
            % into the nominal conjunction plane used in the graphs)
            eigV1 = eigV1/norm(eigV1);
            clockAng = atan2(-eigV1(2),eigV1(1)) * 180/pi;
            % Adjust the position based on where the secondary should be
            % positioned
            xOffset = pnrm;
            % Adjust the clock angle and offset if we're using the alternative
            % CA dist plane
            if params.alt_ca_dist == 1
                clockAng = 180 - clockAng;
                xOffset = 0;
            end
            % Generate an ellipsed sized by the 3-sigma values, rotated by the
            % clock angle, and centered on the nominal miss distance (which
            % should always be on the negative x-axis)
            [xellipsePri, yellipsePri] = GenerateEllipsePoints(a, b, clockAng, [xOffset 0]);
            % Get the semi-major eigenvector and eigenvalues from the projected
            % covariance
            eigV1 = pcInfo.Pc2DInfo.EigV1Sec;
            eigL1 = pcInfo.Pc2DInfo.EigL1Sec;
            eigL2 = pcInfo.Pc2DInfo.EigL2Sec;
            % Get the semi-major and semi-minor axis 3-sigma values
            a = sqrt(eigL1)*3;
            b = sqrt(eigL2)*3;
            % Calculate the clock angle in the conjunction plane coordinates
            % (includes a rotation from the conjunction plane used by 2D Pc
            % into the nominal conjunction plane used in the graphs)
            eigV1 = eigV1/norm(eigV1);
            clockAng = atan2(-eigV1(2),eigV1(1)) * 180/pi;
            % Adjust the position based on where the secondary should be
            % positioned
            xOffset = pnrm;
            % Adjust the clock angle and offset if we're using the alternative
            % CA dist plane
            if params.alt_ca_dist == 1
                clockAng = 180 - clockAng;
                xOffset = 0;
            end
            % Generate an ellipsed sized by the 3-sigma values, rotated by the
            % clock angle, and centered on the nominal miss distance (which
            % should always be on the negative x-axis)
            [xellipseSec, yellipseSec] = GenerateEllipsePoints(a, b, clockAng, [xOffset 0]);
        else
            PlotPriSecEllipses = false;
        end
        
        % Get the semi-major eigenvector and eigenvalues from the projected
        % covariance
        eigV1 = pcInfo.Pc2DInfo.EigV1;
        eigL1 = pcInfo.Pc2DInfo.EigL1;
        eigL2 = pcInfo.Pc2DInfo.EigL2;
        
        % Get the semi-major and semi-minor axis 3-sigma values
        a = sqrt(eigL1)*3;
        b = sqrt(eigL2)*3;
        
        % Calculate the clock angle in the conjunction plane coordinates
        % (includes a rotation from the conjunction plane used by 2D Pc
        % into the nominal conjunction plane used in the graphs)
        eigV1 = eigV1/norm(eigV1);
        clockAng = atan2(-eigV1(2),eigV1(1)) * 180/pi;
        
        % Adjust the position based on where the secondary should be
        % positioned
        xOffset = pnrm;
        
        % Adjust the clock angle and offset if we're using the alternative
        % CA dist plane
        if params.alt_ca_dist == 1
            clockAng = 180 - clockAng;
            xOffset = 0;
        end
        
        % Generate an ellipsed sized by the 3-sigma values, rotated by the
        % clock angle, and centered on the nominal miss distance (which
        % should always be on the negative x-axis)
        [xellipse, yellipse] = GenerateEllipsePoints(a, b, clockAng, [xOffset 0]);
        
        % Save the clock angle for the rotated CA dist., and rotate the
        % ellipse points
        if params.alt_ca_dist == 2
            angleCAdist = atan2(eigV1(2),eigV1(1));
            [xellipse, yellipse] = rotateCApoints( ...
                angleCAdist, xellipse, yellipse);
            if PlotPriSecEllipses
                [xellipsePri, yellipsePri] = rotateCApoints(...
                    angleCAdist, xellipsePri, yellipsePri);
                [xellipseSec, yellipseSec] = rotateCApoints(...
                    angleCAdist, xellipseSec, yellipseSec);
            end
        end
        
    end
    
    % Rotate conj. plane points calculated earlier
    if angleCAdist ~= 0
        [Xca0, Yca0] = rotateCApoints(angleCAdist, Xca0, Yca0);
        [Xca, Yca] = rotateCApoints(angleCAdist, Xca, Yca);
    end
    
    %% Calculate the contour depicting the 2D-Nc Nsigma contour
    xcontour2DNc = []; ycontour2DNc = [];
    % Determine if 3D-Nc contour can be calculated and plotted
    if params.plot_Pc2D_Hall
        if params.AuxCAdistContour ~= 3
            AuxCAdistContour = 0;
        elseif params.alt_ca_dist == 1
            AuxCAdistContour = 0;
            warning('Combined AuxCAdistContour and alt_ca_dist==1 options are currently unavailable');
        else
            [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('Pc2D_Hall',params);
            if pltEnabled
                [Nc2DViolations, ~] = GetNumberOfViolations(2, pcInfo, name);
            else
                Nc2DViolations = 0;
            end
            AuxCAdistContour = Nc2DViolations == 0;
            % if pcInfo.covXcorr_corrections_applied && AuxCAdistContour
            %     AuxCAdistContour = 0;
            %     warning('2D-Nc method contour not available with covariance cross corrections');
            % end
        end
        % Generate 2D-Nc contour if required
        if AuxCAdistContour
            % Nominal (2D-Pc) conj plane (CP) x axis along rel miss vector
            xhat0 = pcInfo.Nc2DInfo.Xmean20(1:3)-pcInfo.Nc2DInfo.Xmean10(1:3);
            xhat0 = xhat0/sqrt(xhat0'*xhat0);
            % Use 2D-Nc effective state to define 2D-Nc CP ref frame axes
            zhat = pcInfo.Nc2DInfo.vCAeff; zhat = zhat/sqrt(zhat'*zhat);
            xhat = xhat0-zhat*(zhat'*xhat0); xhat = xhat/sqrt(xhat'*xhat);
            yhat = cross(zhat,xhat);
            Mrot = [xhat yhat zhat]';
            % Miss location vector and cov matrix in CP frame
            rmiss = Mrot*pcInfo.Nc2DInfo.rCAeff;
            Amiss = Mrot*pcInfo.Nc2DInfo.covCAeff*Mrot';
            [VV,LL] = eig(Amiss(1:2,1:2)); LL = diag(LL);
            % Get the semi-major and semi-minor axis 3-sigma values
            aa = sqrt(LL(1))*3;
            bb = sqrt(LL(2))*3;
            % Calculate the clock angle in the conjunction plane coordinates
            % (includes a rotation from the conjunction plane used by 2D Pc
            % into the nominal conjunction plane used in the graphs)
            V1 = VV(:,1); V1 = V1/sqrt(V1'*V1);
            cAng = atan2(-V1(2),V1(1)) * 180/pi;
            % Generate an ellipsed sized by the 3-sigma values, rotated by the
            % clock angle, and centered on the nominal miss distance (which
            % should always be on the negative x-axis)
            [xcontour2DNc, ycontour2DNc] = ...
                GenerateEllipsePoints(aa, bb, cAng, [rmiss(1) -rmiss(2)]);
        end
        % Rotate conj. plane points
        if angleCAdist ~= 0 && AuxCAdistContour
            [xcontour2DNc,ycontour2DNc] = rotateCApoints( ...
                angleCAdist,xcontour2DNc,ycontour2DNc);
        end
        
    end    
    
    %% Calculate the contour depicting the 3D-Nc Nsigma contour
    xcontour3DNc = []; ycontour3DNc = [];
    % Determine if 3D-Nc contour can be calculated and plotted
    if params.plot_Pc3D_Hall && isempty(xcontour3DNc)
        if params.AuxCAdistContour == 0
            AuxCAdistContour = 0;
        elseif params.alt_ca_dist == 1
            AuxCAdistContour = 0;
            warning('Combined AuxCAdistContour and alt_ca_dist==1 options are currently unavailable');
        else
            [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('PcCircle',params);
            if pltEnabled
                [Pc2DViolations, ~] = GetNumberOfViolations(1, pcInfo, name);
            else
                Pc2DViolations = 0;
            end
            [pltEnabled, ~, ~, ~, ~, ~, ~, name] = GetPlotFormat('Pc3D_Hall',params);
            if pltEnabled
                [Nc3DViolations, ~] = GetNumberOfViolations(3, pcInfo, name);
            else
                Nc3DViolations = 0;
            end
            if params.AuxCAdistContour == 2
                AuxCAdistContour = Nc3DViolations == 0;
            else % if params.AuxCAdistContour == 1
                AuxCAdistContour = Pc2DViolations > 0 & Nc3DViolations == 0;
            end
            if pcInfo.covXcorr_corrections_applied && AuxCAdistContour
                AuxCAdistContour = 0;
                warning('3D-Nc method contour not available with covariance cross corrections');
            end
        end
        % Generate 3D-Nc contour if required
        if AuxCAdistContour
            [xcontour3DNc,ycontour3DNc] = CADistNsigmaContour( ...
                pcInfo.Nc3DInfo.Xmean10(1:3)*1e3,      ...
                pcInfo.Nc3DInfo.Xmean10(4:6)*1e3,      ...
                pcInfo.Nc3DInfo.Pmean10*1e6,           ...
                pcInfo.Nc3DInfo.Xmean20(1:3)*1e3,      ...
                pcInfo.Nc3DInfo.Xmean20(4:6)*1e3,      ...
                pcInfo.Nc3DInfo.Pmean20*1e6);
        end
        % Rotate conj. plane points
        if angleCAdist ~= 0 && AuxCAdistContour
            [xcontour3DNc,ycontour3DNc] = rotateCApoints( ...
                angleCAdist,xcontour3DNc,ycontour3DNc);
        end
        
    end
    
    %% Get the minimum and maximum extents for the graphs
    HBR = pcInfo.HBR;
    % Set axis limits
    if params.alt_ca_dist == 1
        % Limits for the alternative CA dist plane
        xmn0 = min([-HBR+pnrm min(Xca) min(xellipse)]);
        xmx0 = max([ HBR+pnrm max(Xca) max(xellipse)]);
        ymn0 = min([HBR min(Yca) min(yellipse)]);
        ymx0 = max([HBR max(Yca) max(yellipse)]);    
    else
        if ~isempty(xcontour2DNc) || ~isempty(xcontour3DNc) 
            % Set limits excluding the MC positions, for repeatability
            xmn0 = min([-HBR min(xellipse) min(xcontour2DNc) min(xcontour3DNc)]);
            xmx0 = max([ HBR max(xellipse) max(xcontour2DNc) max(xcontour3DNc)]);
            ymn0 = min([-HBR min(yellipse) min(ycontour2DNc) min(ycontour3DNc)]);
            ymx0 = max([ HBR max(yellipse) max(ycontour2DNc) max(ycontour3DNc)]);
        else
            % Set limits including the MC positions
            xmn0 = min([-HBR min(Xca) min(xellipse)]);
            xmx0 = max([ HBR max(Xca) max(xellipse)]);
            ymn0 = min([-HBR min(Yca) min(yellipse)]);
            ymx0 = max([ HBR max(Yca) max(yellipse)]);
        end
    end

    %% Set parameters for each plot type
    if plotType == 1
        % Square axis - no zoom
        subplot(2,3,1);
        xymx = max([xmx0 ymx0 abs(xmn0) abs(ymn0)]);
        xrng = plot_range([-xymx xymx],params.fig.xpad);
        yrng = xrng;
        yAxisLoc = 'left';
        showLegend = false;
    elseif plotType == 2
        % Square axis - Zoom in on HBR
        subplot(2,3,3);
        xrng = [-HBR HBR] * 2;
        yrng = xrng;
        if params.alt_ca_dist == 1
            xrng = xrng + pnrm;
        end
        yAxisLoc = 'right';
        showLegend = false;
    elseif plotType == 3
        % Rectangular axis - no zoom
        subplot(2,3,4:6);
        xrng = plot_range([xmn0 xmx0],params.fig.xpad);
        yrng = plot_range([ymn0 ymx0],params.fig.ypad);
        yAxisLoc = 'left';
        showLegend = true;
    else
        error('Undefined plot type entered');
    end

    %% Scales
    if max(abs(xrng)) >= 1000
        Lxscl = 1e3;
        Lxstr = 'km';
    else
        Lxscl = 1;
        Lxstr = 'm';
    end
    if max(abs(yrng)) >= 1000
        Lyscl = 1e3;
        Lystr = 'km';
    else
        Lyscl = 1;
        Lystr = 'm';
    end

    %% Start graphing data
    hold on;
    SetLabelsAndWeights(['X_{CA} (' Lxstr ')'], ['Y_{CA} (' Lystr ')'], params);
    set(gca,'YAxisLocation',yAxisLoc);
    box on;
    nlgnd = 0;
    lgndObj = [];
    lgndTxt = {};
    
    if params.AuxCAdistContour == 3 && ~isempty(xcontour2DNc)
        AuxCAdistContour = true;
        acontour = '2D-Nc';
        xcontour = xcontour2DNc;
        ycontour = ycontour2DNc; 
    elseif ~isempty(xcontour3DNc)
        AuxCAdistContour = true;
        acontour = '3D-Nc';
        xcontour = xcontour3DNc;
        ycontour = ycontour3DNc; 
    else
        AuxCAdistContour = false;
        acontour = '';
    end
    
    if params.plot_PcCircle
        % Shade in the projected 2D Pc covariance ellipse
        magenta = [1 0 1];
        inten = 0.2;
        col_shade_ellipse = inten*magenta + (1-inten)*[1 1 1];
        nlgnd=nlgnd+1;
        lgndObj(nlgnd) = patch(xellipse/Lxscl,yellipse/Lyscl,col_shade_ellipse, ...
            'EdgeColor',magenta,'LineWidth',params.fig.lnwd1);
        lgndTxt{nlgnd} = '2D-Pc method 3\sigma ellipse';
        if AuxCAdistContour
            lgndTxt{nlgnd} = cat(2,lgndTxt{nlgnd},'   ');
        end
    end
    
    if AuxCAdistContour
        % Shade in the projected auxiliary covariance contour
        gold = [1 0.667 0];
        inten = 0.2;
        col_shade_contour = inten*gold + (1-inten)*[1 1 1];
        nlgnd=nlgnd+1;
        lgndObj(nlgnd) = patch(xcontour/Lxscl,ycontour/Lyscl,col_shade_contour, ...
            'EdgeColor',gold,'LineWidth',params.fig.lnwd1);
        lgndTxt{nlgnd} = [acontour ' method 3\sigma contour'];
    end

    % Plot lines through the origin
    zrng = zeros(size(xrng));
    plot(xrng/Lxscl,zrng/Lyscl,'-','Color',params.fig.gray);
    plot(zrng/Lxscl,yrng/Lyscl,'-','Color',params.fig.gray);
    
    % if params.alt_ca_dist == 2
    %     plot(xrng/Lxscl,[Yca0 Yca0]/Lyscl,'--','Color',params.fig.gray);
    %     plot([Xca0 Xca0]/Lxscl,yrng/Lyscl,'--','Color',params.fig.gray);
    % end

    if params.plot_Pc_SDMC
        % Plot all misses in blue
        mrkr = '.';
        msiz = params.fig.msiz+(1/params.fig.scaleFactor);
        msiz2 = msiz*2;
        mcol = 'b';
        nlgnd=nlgnd+1;
        plot(Xca(~hitInd)/Lxscl,Yca(~hitInd)/Lyscl, ...
            'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
            'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);
        % 2nd "plot" makes sure the marker is visibile in the legend
        lgndObj(nlgnd) = plot(nan,nan, ...
            'Marker',mrkr,'MarkerSize',msiz2,'LineStyle','none', ...
            'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);
        lgndTxt{nlgnd} = 'MC misses';

        % Plot all hits in red
        mcol = 'r';
        nlgnd=nlgnd+1;
        plot(Xca(hitInd)/Lxscl,Yca(hitInd)/Lyscl, ...
            'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
            'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);
        % 2nd "plot" makes sure the marker is visibile in the legend
        lgndObj(nlgnd) = plot(nan,nan, ...
            'Marker',mrkr,'MarkerSize',msiz2,'LineStyle','none', ...
            'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);
        lgndTxt{nlgnd} = 'MC hits';
    end

    % Plot the projected covariance ellipse in magenta
    if params.plot_PcCircle
        
        if PlotPriSecEllipses
            % lcol=[1 0 1 0.4];
            lcol = [0 0 0];
            plot(xellipsePri/Lxscl,yellipsePri/Lyscl,'--',...
                'Color',lcol, ...
                'LineWidth',params.fig.lnwd1+1);
            lcol = [0 0 0];
            plot(xellipseSec/Lxscl,yellipseSec/Lyscl,':',...
                'Color',lcol, ...
                'LineWidth',params.fig.lnwd1+1);
        end
        
        lcol=[1 0 1 0.4];
        plot(xellipse/Lxscl,yellipse/Lyscl,'-',...
            'Color',lcol, ...
            'LineWidth',params.fig.lnwd1);

    end
    
    if AuxCAdistContour
        lcol=[1 0.667 0 0.4];
        plot(xcontour/Lxscl,ycontour/Lyscl,'-',...
            'Color',lcol, ...
            'LineWidth',params.fig.lnwd1);
    end
    
    % Plot the HBR boundary in green
    if params.alt_ca_dist ~= 1
        [xcirc,ycirc] = GenerateEllipsePoints(HBR);
    else
        [xcirc,ycirc] = GenerateEllipsePoints(HBR, HBR, 0, [pnrm 0]);
    end
    nlgnd=nlgnd+1;
    lgndObj(nlgnd) = plot(xcirc/Lxscl,ycirc/Lyscl,'-g', ...
        'LineWidth',params.fig.lnwd1);
    lgndTxt{nlgnd} = 'HBR';
    
    % Plot nominal CA as black diamond
    mrkr = 'd';
    msiz = params.fig.msiz;
    mcol = 'k';
    nlgnd=nlgnd+1;
    lgndObj(nlgnd) = plot(Xca0/Lxscl,Yca0/Lyscl, ...
        'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
        'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);
    lgndTxt{nlgnd} = 'Nominal miss';

    if showLegend
        legend(lgndObj,lgndTxt, ...
            'FontSize',params.fig.lfsz,'FontWeight',params.fig.lfwt, ...
            'Location','NorthOutside','Orientation','horizontal');
        legend boxoff;
    end

    % Set the limits and square the axes, if needed
    xlim(xrng/Lxscl);
    ylim(yrng/Lyscl);
    if isequal(xrng,yrng)
        axis square;
    end
    
    % Reverse the y-axis, if needed
    if params.alt_ca_dist == 2 && Yca0 < 0
        set(gca,'Ydir','reverse');
    end
    
end

%% Calculate the close approach coordinates in the CA coordinate system
function [Xca, Yca, phat, pnrm] = CalcCloseApproachCoordinates(r1, v1, r2, v2, params, phat)
    if nargin == 5
        % pvec is normally the relative position from the primary to the
        % secondary
        pvec = (r2 - r1)';
        if params.alt_ca_dist == 1
            % If we're using the alternative CA dist plane, then
            % pvec is the opposite of the normal placement
            pvec = -pvec;
        end
        pnrm = sqrt(pvec(1)^2 + pvec(2)^2 + pvec(3)^2);
        if pnrm < 1e-10
            error('Nominal TCA miss distance cannot be zero (if possible, make << HBR, but finite)');
        end
        phat = pvec'/pnrm;
    end
    
    % Z-axis is along the negative of the relative velocity between the
    % primary and the secondary
    zvec = -(v2 - v1)';
    if params.alt_ca_dist == 1
        % If we're using the alternative CA dist plane, then the
        % z-axis is in the opposite direction
        zvec = -zvec;
    end
    znrm = sqrt(zvec(1)^2 + zvec(2)^2 + zvec(3)^2);
    if znrm < 1e-10
        error('Relative velocity cannot be zero');
    end
    zhat = zvec / znrm;
    
    % Y-axis is normal to the plane made by the relative position vector
    % and the relative velocity vector
    yhat = cross(phat,-zhat);
    ynrm = sqrt(yhat(1)^2 + yhat(2)^2 + yhat(3)^2);
    if ynrm < 1e-10
        error('Nominal TCA miss distance and relative velocity cannot be parallel');
    end
    yhat = yhat / ynrm;
    
    % X-axis makes up the triad, which should normally place the secondary
    % on the positive x-axis
    xhat = cross(yhat, zhat);
    
    % Conjunction plane (X,Y) coordinates for this sampling trial
    dr = (r2 - r1);
    Xca = dr * xhat';
    Yca = dr * yhat';
end

%% Count up the number of violations for each Pc type
function [nv, numViolationsTxt] = GetNumberOfViolations(pcType, pcInfo, name)
    nv = 0;
    if pcType == 1
        if pcInfo.Pc2DViolations.NPD;             nv=nv+1; end
        if pcInfo.Pc2DViolations.Extended;        nv=nv+1; end
        if pcInfo.Pc2DViolations.Offset;          nv=nv+1; end
        if pcInfo.Pc2DViolations.Inaccurate;      nv=nv+1; end
    elseif pcType == 2
        if pcInfo.Nc2DInfo.Violations.Extended;   nv=nv+1; end
        if pcInfo.Nc2DInfo.Violations.Offset;     nv=nv+1; end
        if pcInfo.Nc2DInfo.Violations.Inaccurate; nv=nv+1; end
    elseif pcType == 3
        if pcInfo.Nc3DInfo.Violations.Extended;   nv=nv+1; end
        if pcInfo.Nc3DInfo.Violations.Offset;     nv=nv+1; end
        if pcInfo.Nc3DInfo.Violations.Inaccurate; nv=nv+1; end
    elseif pcType == 4
        if pcInfo.SDMCInfo.Violations.Extended;   nv=nv+1; end
        if pcInfo.SDMCInfo.Violations.Offset;     nv=nv+1; end
        if pcInfo.SDMCInfo.BFMCsubstitution
            name = strrep(name,'SDMC','BFMC');
        end
    else
        error('Invalid pcType passed in');
    end
    if nv == 0
        numViolationsTxt = [name ' Violations: None'];
    else
        numViolationsTxt = ['\color{red}' name ' Violations: ' num2str(nv) '\color{black}'];
    end
end

%% Rotate conjunction plane points 

function [xrot, yrot] = rotateCApoints(angle, x, y)

% Rotate 2-D CA distribution points

% Return immediately for zero rotation
if angle == 0
    xrot = x;
    yrot = y;
    return;
end

% Check for equal dimensions
S = size(x);
if ~isequal(S,size(y))
    error('Unequal x and y array dimensions');
elseif numel(angle) ~= 1
    error('Rotation angle must be a scalar');
elseif ~isreal(angle)
    error('Rotation angle must be real');
end

% Return for zero points
N = numel(x); 
if N == 0
    xrot = x;
    yrot = y;
    return;
end

% Cosine and sine values
c = cos(angle); s = sin(angle);

% Calculate rotated x and y arrays
xrot = c*x - s*y;
yrot = s*x + c*y;

return
end

%%
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-24 | Initial Development
% L. Baars       | 2023-Apr-11 | Added call to FindNearbyCA before running
%                                CalcCloseApproachCoordinates.
% L. Baars       | 2023-Apr-20 | Added alt_ca_dist option to graph an
%                                alternate distribution used by Ops.
% L. Baars       | 2023-May-19 | Added data quality and invalid 6x6 checks.
% L. Baars       | 2023-Jul-25 | Added figInfo output
% D. Hall        | 2023-Aug-30 | Changed errors to warnings when plots
%                                cannot be made. Added plot_conjID_string
%                                input parameter. Updated text displayed in
%                                plot outputs.
% L. Baars       | 2024-Jan-11 | Changed references to PcConjPlaneCircle
%                                into PcCircle. Also, tool no longer
%                                returns early if a plot cannot be made.
% D. Hall        | 2024-Feb-07 | Added code to plot the primary and
%                                secondary covariance ellipses for
%                                the unequally-zoomed CA distribution
%                                (i.e., plotType = 3); only performed if
%                                required information is available in 
%                                pcInfo.Pc2DInfo structure (i.e., EigV1Pri,
%                                EigL1Pri, etc.) as produced using the
%                                PcCircle function with the flag
%                                "PriSecCovProcessing" flag set to true.
% D. Hall        | 2024-Sep-18 | Added support for plotting BFMC data.
% D. Hall        | 2024-Nov-18 | Added code to allow a new alternative
%                                CA distribution plane orientation. When 
%                                params.alt_ca_dist is set to 1 or 2.
%                                If 1, the combined covariance is placed at
%                                the origin. If 2, then the plane is
%                                rotated to align the axes with with 
%                                principal axes of the 2D-Pc ellipse.
%                                Added params.AuxCAdistContour parameter to
%                                plot a 3D-Nc or 2D-Nc method 3-sigma
%                                contour, in addition to the 2D-Pc ellipse.
% E. Toumey      | 2025-Mar-02 | Moved file for new directory structure.
% L. Baars       | 2025-Aug-28 | Updated code for public release.
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
