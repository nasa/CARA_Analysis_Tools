function [Pc,out] = PcMultiStepWithPlots(r1,v1,C1,r2,v2,C2,HBR,params)
% PcMultiStepWithPlots - Executes the function PcMultiStepWithSDMC and
%                        makes temporal and close approach distribution
%                        plots.
%
% Syntax: [Pc, out] = PcMultiStepWithPlots(r1,v1,C1,r2,v2,C2,HBR);
%         [Pc, out] = PcMultiStepWithPlots(r1,v1,C1,r2,v2,C2,HBR,params);
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
% This tool will call PcMultiStepWithSDMC and then will plot the temporal
% and close approach distribution plots for the conjunction. See Hall et.
% al. (2018) for a discussion of the temporal plots and Hall et. al. (2023)
% for a discussion of the CA distribution plots.
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
%              the structure defined by PcMultiStepWithSDMC.m. The
%              following fields are specific to this function:
%
%      generate_conj_plots - Boolean array which         [1x1, 1x2, or 2x1]
%                            indicates the plots to be generated. The first
%                            element turns on/off temporal conjunction
%                            plots and the second element turns on/off
%                            close approach (CA) distribution plots. If a
%                            single boolean values is passed in, then that
%                            value is used for both plots.
%                            Defaults to [true true].
%
%      figure_scale_factor - Scales the fonts and plot elements by
%                            1/scaleFactor for OS screen scales that aren't
%                            set to 100%. For example, if the screen scale
%                            is 125%, set the scaleFactor to 1.25 to best
%                            replicate nominal CARA plots.
%                            Defaults to 1.0.
%
%      ForceSDMCCalculation - This parameter has the same meaning as in
%                             PcMultiStepWithSDMC, except with a new
%                             default.
%                             Defaults to true.
%
%      PreventSDMCCalculation - If both the "Force" and "Prevent" SDMC
%                               parameters are set, the "Prevent" parameter
%                               will take precedence. This will allow a
%                               user easily turn off SDMC graphing if they
%                               want to speed up the generation process.
%                               However, it is recommended to set this
%                               value to false if any usage violations are
%                               found. Note: The "Prevent" parameters for
%                               the other Pc methods are ignored within
%                               this function.
%                               Defaults to false.
%
%      plot_conjID_string - The conjunction ID string to display at the top
%                           of the plots. It will also be the first part of
%                           the saved file name. If this field is empty,
%                           then no conjunction ID is displayed within the
%                           plots and the saved file name will begin with
%                           'UndefinedConj'.
%                           Defaults to ''.
%
%      plot_save_tag - An optional tag that is appended to the
%                      plot_conjID_string when saving the file name. This
%                      tag is not displayed in the plot title.
%                      Defaults to ''.
%
%      plot_save_loc - Location where plots will be saved.
%                      Defaults to the current working directory.
%
%      plot_save_ext - Extension used for saving the figure outputs. Valid
%                      possible values are documented in Matlab's
%                      "Extension" section of the "filename" option for the
%                      "saveas" function.
%                      Defaults to '.png'.
%
%      priLastObsAge - The amount of time (in days) between the last
%                      observation of the primary object and TCA.
%                      Defaults to NaN.
%
%      secLastObsAge - The amount of time (in days) between the last
%                      observation of the secondary object and TCA.
%                      Defaults to NaN.
%
%      alt_ca_dist - Allows the generation of the CA distribution plot with
%                    alternate configurations.
%                      0 = The HBR is placed at the origin and the nominal
%                          miss distance is placed on the positive x-axis.
%                          The z-axis is the opposite of the relative
%                          velocity vector and the y-axis makes up the
%                          right handed coordinate system triad. This is
%                          considered the nominal conjunction plane.
%                      1 = The nominal miss is placed at the origin and the
%                          HBR is placed on the positive x-axis. The z-axis
%                          is the opposite of the nominal conjunction
%                          plane, and the y-axis makes up the right handed
%                          coordinate system triad. This conjunction plane
%                          is a 180 degree rotation about the y-axis of the
%                          nominal conjunction plane and is translated in
%                          the positive x-direction by the nominal miss
%                          distance. This configuration is the format for
%                          historical CARA CA distribution plots.
%                      2 = In this configuration, the conjunction plane is
%                          rotated to align the axes with the principal
%                          axes of the 2D-Pc ellipse.
%                    Defaults to 0.
%
%      AuxCAdistContour - Provides an option to turn on a beta-version
%                         analytical visualization of the CA distribution
%                         contour.
%                           0 = Never plot auxiliary contours for the 3D-Nc
%                               or 2D-Nc methods
%                           1 = Plot 3D-Nc contour, but only for 2D-Pc
%                               usage violations
%                           2 = Always plot 3D-Nc contour (the "bananoid"),
%                               when possible (for testing and R&D)
%                           3 = Always plot 2D-Nc contour, when possible
%                               (for testing and R&D)
%                         Note: As of the current release, options 1-3 are
%                               disabled in the public-release version of
%                               the software.
%
% =========================================================================
%
% Outputs:
%
%    Pc - Recommended Pc value to use for this conjunction
%
%    out - Auxiliary output information structure which adds to the "out"
%          structure provided by PcMultiStepWithSDMC.m. This function adds
%          the following output fields:
%
%      pcTimePlotFile - The name of the temporal plot file created by this
%                       process.
%
%      caDistPlotFile - The name of the CA distribution plot file created
%                       by this process.
%
% =========================================================================
%
% References:
%
%    D.Hall, S.Casali, L.Johnson, B.Skrehart, and L.Baars (2018) "High
%    Fidelity Collision Probabilities Estimated using Brute Force Monte
%    Carlo Simulations" AAS 18-244.
%
%    D.Hall, L.Baars, and S.Casali (2023) "A Multistep Probability of
%    Collision Algorithm" AAS 23-398.
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
    persistent pathsAdded researchCodeAvailable
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        s = what(fullfile(p, '../Utils/Plotting')); addpath(s.path);
        s = what(fullfile(p, '../Utils/Plotting/Utils')); addpath(s.path);
        try
            s = what(fullfile(p, '../../ResearchCode/OCMDB_PcMultiStep/src')); addpath(s.path);
            s = what(fullfile(p, '../../ResearchCode/Utils/Plotting')); addpath(s.path);
            researchCodeAvailable = true;
        catch
            researchCodeAvailable = false;
        end
        pathsAdded = true;
    end

    % Default inputs
    Nargin = nargin;
    if Nargin < 8; params = []; end
    
    % Flags to make temporal and close approach (CA) distribution plots
    params = set_default_param(params,'generate_conj_plots',[true true]);
    if numel(params.generate_conj_plots) == 1
        params.generate_conj_plots = repmat(params.generate_conj_plots,[1 2]);
    end
    
    % If making the conjunction plots, set the required parameters
    any_conj_plots = any(params.generate_conj_plots);
    if any_conj_plots
        params.generate_time_plot    = params.generate_conj_plots(1);
        params.generate_ca_dist_plot = params.generate_conj_plots(2);        
        % Always force 2D-Nc and 3D-Nc when plotting
        params.ForceNc2DCalculation = true; params.PreventNc2DCalculation = false;
        params.ForceNc3DCalculation = true; params.PreventNc3DCalculation = false;
        % Instead of always forcing SDMC, use the specified parameters but
        % give the PreventSDMCCalculation parameter precedence
        params = set_default_param(params,'PreventSDMCCalculation',false);
        if params.PreventSDMCCalculation
            params.ForceSDMCCalculation = false;
        else
            params = set_default_param(params,'ForceSDMCCalculation',true);
        end
        % Other plot params
        params = set_default_param(params,'figure_scale_factor',1.0);
        params = set_default_param(params,'plot_save_loc',pwd);
        params = set_default_param(params,'plot_save_tag','');
        params = set_default_param(params,'plot_save_ext','.png');
        params = set_default_param(params,'plot_conjID_string','');
        params = set_default_param(params,'priLastObsAge',NaN);
        params = set_default_param(params,'secLastObsAge',NaN);
        params = set_default_param(params,'alt_ca_dist',0);
        % Default for plotting the 3D-Nc method CA dist contour
        % 0 = Never; 1 = Only for 2D-Pc usage violations; 2 = Always when possible
        AuxCAdistContourDefault = 0; 
        params = set_default_param(params,'AuxCAdistContour',AuxCAdistContourDefault);
        if ~(isequal(params.AuxCAdistContour,0) || ...
             isequal(params.AuxCAdistContour,1) || ...
             isequal(params.AuxCAdistContour,2) || ...
             isequal(params.AuxCAdistContour,3))
            warning(['Invalid AuxCAdistContour value; setting to default = ' num2str(AuxCAdistContourDefault)]);
            params.AuxCAdistContour = AuxCAdistContourDefault;
        end
        if ~researchCodeAvailable && params.AuxCAdistContour ~= AuxCAdistContourDefault
            error('The selected params.AuxCAdistContour value is not currently supported!');
        end
    else
        params.generate_ca_dist_plot = false;
        params.generate_time_plot    = false;
    end
    
    % Path to find BFMC data to substitute for SDMC data
    params = set_default_param(params,'BFMCpath','');
    
    %% Use PcMultiStepWithSDMC to calculate Pc2D, Nc2D, Nc3D and PcSDMC
    
    % Call PcMultiStep
    [Pc,out] = PcMultiStepWithSDMC(r1,v1,C1,r2,v2,C2,HBR,params);
    
    % Return now if not making plots
    if ~any_conj_plots
        return;
    end
    
    % Return if any data quality errors were detected
    if out.AnyDataQualityErrors
        warning('PcMultiStepWithPlots plots not created due to data quality');
        return;
    end
    
    %% Make the plots
    
    if any_conj_plots

        % Copy the output structure for plotting
        PcMSWout = out;

        % Set up plotting parameters
        if isfield(PcMSWout,'SDMCParams')
            PlotPars = PcMSWout.SDMCParams;
        else
            PlotPars = struct();
        end
        PlotPars.generate_time_plot = params.generate_time_plot;
        PlotPars.generate_ca_dist_plot = params.generate_ca_dist_plot;
        PlotPars.fig.scaleFactor = params.figure_scale_factor;
        PlotPars.plot_PcCircle = true;
        PlotPars.plot_Pc2D_Hall = true;
        PlotPars.plot_Pc3D_Hall = true;
        PlotPars.plot_Pc_SDMC = ~params.PreventSDMCCalculation;
        PlotPars.AuxCAdistContour = params.AuxCAdistContour;
        PlotPars.plot_save_loc = params.plot_save_loc;
        PlotPars.plot_save_tag = params.plot_save_tag;
        PlotPars.plot_save_ext = params.plot_save_ext;
        PlotPars.plot_conjID_string = params.plot_conjID_string;
        PlotPars.priLastObsAge = params.priLastObsAge;
        PlotPars.secLastObsAge = params.secLastObsAge;
        PlotPars.alt_ca_dist = params.alt_ca_dist;

        % If a plot_save_tag is used and it doesn't begin with '_'
        if ~isempty(PlotPars.plot_save_tag) && ~startsWith(PlotPars.plot_save_tag,'_')
            % Add the '_' to the beginning so that the tag appears in the
            % appropriate spot in the file name when figures are saved
            PlotPars.plot_save_tag = ['_' PlotPars.plot_save_tag];
        end
        if PcMSWout.covXcorr_corrections_applied
            PlotPars.plot_save_tag = ...
                cat(2,PlotPars.plot_save_tag,'_wCCC');
        end
        if isfield(PcMSWout,'SDMCInfo')
            if PcMSWout.SDMCInfo.BFMCsubstitution
                PlotPars.plot_save_tag = ...
                    cat(2,PlotPars.plot_save_tag,'_wBFMC');
            else
                PlotPars.plot_save_tag = ...
                    cat(2,PlotPars.plot_save_tag,'_wSDMC');
            end
        end

        % Conjunction states and covariances
        PcMSWout.r1 = r1;
        PcMSWout.v1 = v1;
        PcMSWout.C1 = C1;
        PcMSWout.r2 = r2;
        PcMSWout.v2 = v2;
        PcMSWout.C2 = C2;

        % Make the temporal plot
        if params.generate_time_plot
            out.pcTimePlotFile = Pc_Time_Plot(PlotPars, PcMSWout);
        else
            out.pcTimePlotFile = '';
        end

        % Make the CA distribution plot
        if params.generate_ca_dist_plot
            out.caDistPlotFile = CA_Dist_Plot(PlotPars, PcMSWout);
        else
            out.caDistPlotFile = '';
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
% D. Hall        | 2023-AUG-21 | Initial Development.
% L. Baars       | 2024-JAN-11 | Standardized "out" structure
% D. Hall        | 2024-SEP-16 | Added support for BFMC data processing
% D. Hall        | 2024-OCT-19 | Enabled user control of SDCM plotting by
%                                adding the prevent of force options.
% D. Hall        | 2024-NOV-19 | Added options for alternative CA dist and
%                                aux CA dist contour plot
% L. Baars       | 2025-AUG-25 | Updated code for public release. Added
%                                figure_scale_factor parameter.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
