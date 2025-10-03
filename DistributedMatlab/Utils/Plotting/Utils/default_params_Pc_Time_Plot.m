function params = default_params_Pc_Time_Plot(params)
% default_params_Pc_Time_Plot - Add and/or set the defaults for the
%                               parameters used by the function
%                               Pc_Time_Plot.
%
% Syntax: params = default_params_Pc_Time_Plot;
%         params = default_params_Pc_Time_Plot(params);
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
%   params = (Optional) Empty or partially populated structure containing
%            the parameters used by the Pc_Time_Plot function, including
%            the following fields:
%
%     generate_time_plot = Boolean which controls whether or not the Pc
%                          time plot is generated (default = false)
%
%     plot_PcCircle = Plot data from the conjunction plane circle Pc
%                     calculation (default = false)
%
%     plot_Pc2D_Hall = Plot data from the Pc2D_Hall Pc calculation
%                      (default = false)
%
%     plot_Pc3D_Hall = Plot data from the Pc3D_Hall Pc calculation
%                      (default = false)
%
%     plot_Pc_SDMC = Plot data from the Pc_SDMC Pc calculation
%                    (default = false)
%
%     conf_level = Confidence level for the error bands of the Pc_SDMC
%                  plots (default = 0.95)
%
%     PcNsf = Number of significant figures for Pc related labels
%             (default = 4)
%
%     RcNsf = Number of significant figures for Pc rate related labels
%             (default = 3)
%
%     TcNsf = Number of significant figures for time related labels
%             (default = 3)
%
%     verbose = Display verbose outputs when generating the graphs
%               (default = false)
%
%     verboseTextSep = Line separator used between verbose output sections
%
%     trajectory_mode = Indicates the trajectory mode used by SDMC
%                         0 = 2-body motion
%                         1 = full rectilinear motion
%                         2 = rectilinear motion (position deviations only)
%                       (default = 0)
%
%     plot_save_loc = Save location for plot outputs (if blank, the plots
%                     are not saved)
%                     (default = '')
%
%     plot_save_tag = Set a save tag for plot outputs, the tag is placed
%                     within the file name to uniquely identify the saved
%                     file (if blank, the conjunction ID is used as the
%                     tag)
%                     (default = '')
%
%     fig = Sub-structure used to define figure specific
%           parameters, including the following fields:
%
%       dispConjID = Display the conjunction ID in the title of plots
%                    (default = true)
%
%       xpad = Amount of padding (as a percentage of graphed x-values)
%              added to the x-axis (default = 0.05)
%
%       ypad = Amount of padding (as a percentage of graphed y-values)
%              added to the y-axis (default = params.fig.xpad)
%
%       scaleFactor = Scales the fonts and plot elements by 1/scaleFactor
%                     for OS screen scales that aren't set to 100%. For
%                     example, if the screen scale is 125%, set the
%                     scaleFactor to 1.25 to best replicate nominal CARA
%                     plots. (default = 1.0)
%
%       afsz = Axis font size (default = 12)
%
%       afwt = Axis font weight (default = 'bold')
%
%       alwd = Axis line width (default = 2)
%
%       xfsz = x-axis label font size (default = params.fig.afsz)
%
%       xfwt = x-axis label font weight (default = params.fig.afwt)
%
%       yfsz = y-axis label font size (default = params.fig.afsz)
%
%       yfwt = y-axis label font weight (default = params.fig.afwt)
%
%       tfsz = graph title font size (default = params.fig.afsz+1)
%
%       tfwt = graph title font weight (default = 'bold')
%
%       lfsz = legend font size (default = params.fig.afsz)
%
%       lfwt = legend font weight (default = params.fig.afwt)
%
%       lnwd0 = "thin" line width for line plots (default = 2)
%
%       lnwd1 = "medium" line width for line plots
%               (default = params.fig.lnwd0+0.5)
%
%       lnwd2 = "thick" line width for line plots
%               (default = params.fig.lnwd1+1)
%
%       mrkr = marker type for line plots (default = 'none')
%
%       msiz = marker size for line plots (default = 5)
%
%       darkgreen = dark green color (default = [62 150 81]/255)
%
%       col_MC = line color for Monte Carlo data, only used by Pc_SDMC
%                plots (default = "magenta")
%
%       inten_err = error band color intensity as a percentage of col_MC,
%                   only used by Pc_SDMC plots, valid range is 0.0 to 1.0
%                   (default = 0.2)
%
%       col_MC_err = error band color for Monte Carlo data, only used by
%                    Pc_SDMC plots
%                    (default = inten_err percentage of col_MC)
%
% =========================================================================
%
% OUTPUT:
%
%   params = Fully-populated structure containing parameters used by the 
%            function Pc_Time_Plot.
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % Initializions and defaults
    Nargin = nargin;
    if Nargin < 1; params = []; end
    
    % Generate time plot
    params = set_default_param(params,'generate_time_plot',false);
    
    % Individual plots to include in the graphs
    params = set_default_param(params,'plot_PcCircle',false);
    params = set_default_param(params,'plot_Pc2D_Hall',false);
    params = set_default_param(params,'plot_Pc3D_Hall',false);
    params = set_default_param(params,'plot_Pc_SDMC',false);
    
    % Confidence level for the error bands of the Pc_SDMC plots
    params = set_default_param(params,'conf_level',0.95);
    
    % Number of significant digits to display for each data type
    params = set_default_param(params,'PcNsf',4);                          % Pc values
    params = set_default_param(params,'RcNsf',3);                          % Rc values
    params = set_default_param(params,'TcNsf',3);                          % Time values
    
    % Controls for verbose outputs
    params = set_default_param(params,'verbose',false);                    % Verbose outputs are displayed to STDOUT
    params = set_default_param(params,'verboseTextSep',repelem('-',1,56)); % Line separator between verbose output sections
    
    % Default trajectory mode
    %   0 = 2-body motion
    %   1 = full rectilinear motion
    %   2 = rectilinear motion (position deviations only)
    params = set_default_param(params,'trajectory_mode',0);
    if params.trajectory_mode ~= 0 && params.trajectory_mode ~= 1 && ...
            params.trajectory_mode ~= 2
        error('Supplied trajectory_mode must be 0, 1, or 2');
    end
    
    % Set a save location for plot ouputs (if blank, the plots are not
    % saved)
    params = set_default_param(params,'plot_save_loc','');
    
    % Set a save tag for plot outputs, the tag is placed within the file
    % name to uniquely identify the saved file (if blank, the conjunction
    % ID is used as the tag)
    params = set_default_param(params,'plot_save_tag','');
    
    % Figure specific parameters
    params = set_default_param(params,'fig',[]);                           % "fig" structure within params
    params.fig = set_default_param(params.fig,'dispConjID',true);              % display the conjunction ID in plots
    params.fig = set_default_param(params.fig,'xpad',0.05);                    % percentage of graph data to pad in x-axis
    params.fig = set_default_param(params.fig,'ypad',params.fig.xpad);         % percentage of graph data to pad in y-axis
    params.fig = set_default_param(params.fig,'scaleFactor',1.0);              % scale factor for fonts and plot elements
    params.fig = set_default_param(params.fig,'afsz',12);                      % axis font size
    params.fig = set_default_param(params.fig,'afwt','bold');                  % axis font weight
    params.fig = set_default_param(params.fig,'alwd',2);                       % axis line width
    params.fig = set_default_param(params.fig,'xfsz',params.fig.afsz);         % x-axis labels font size
    params.fig = set_default_param(params.fig,'xfwt',params.fig.afwt);         % x-axis labels font weight
    params.fig = set_default_param(params.fig,'yfsz',params.fig.afsz);         % y-axis labels font size
    params.fig = set_default_param(params.fig,'yfwt',params.fig.afwt);         % y-axis labels font weight
    params.fig = set_default_param(params.fig,'tfsz',params.fig.afsz+1);       % graph title font size
    params.fig = set_default_param(params.fig,'tfwt','bold');                  % graph title font weight
    params.fig = set_default_param(params.fig,'lfsz',params.fig.afsz);         % legend font size
    params.fig = set_default_param(params.fig,'lfwt',params.fig.afwt);         % legend font weight
    params.fig = set_default_param(params.fig,'lnwd0',2);                      % "thin" line width
    params.fig = set_default_param(params.fig,'lnwd1',params.fig.lnwd0+0.5);   % "medium" line width
    params.fig = set_default_param(params.fig,'lnwd2',params.fig.lnwd1+1);     % "thick" line width
    params.fig = set_default_param(params.fig,'mrkr','none');                  % marker type
    params.fig = set_default_param(params.fig,'msiz',5);                       % marker size
    params.fig = set_default_param(params.fig,'darkgreen',[62 150 81]/255);    % dark green color
    
    % Scale figure elements
    params.fig.afsz = params.fig.afsz/params.fig.scaleFactor;
    params.fig.xfsz = params.fig.xfsz/params.fig.scaleFactor;
    params.fig.yfsz = params.fig.yfsz/params.fig.scaleFactor;
    params.fig.tfsz = params.fig.tfsz/params.fig.scaleFactor;
    params.fig.lfsz = params.fig.lfsz/params.fig.scaleFactor;
    params.fig.alwd = params.fig.alwd/params.fig.scaleFactor;
    params.fig.lnwd0 = params.fig.lnwd0/params.fig.scaleFactor;
    params.fig.lnwd1 = params.fig.lnwd1/params.fig.scaleFactor;
    params.fig.lnwd2 = params.fig.lnwd2/params.fig.scaleFactor;
    params.fig.msiz = params.fig.msiz/params.fig.scaleFactor;
    
    % Figure specific parameters for the Pc_SDMC line and error bands
    params.fig = set_default_param(params.fig,'col_MC',[1 0 1]);               % Monte Carlo line color (magenta)
    params.fig = set_default_param(params.fig,'inten_err',0.2);                % Monte Carlo error band color intensity (0.0 to 1.0)
    if params.fig.inten_err < 0.0 || params.fig.inten_err > 1.0
        error('Field params.fig.inten_err must be between 0.0 and 1.0');
    end
    % Determine error band color based on intensity
    col_MC_err = params.fig.inten_err*params.fig.col_MC + (1-params.fig.inten_err)*[1 1 1];
    params.fig = set_default_param(params.fig,'col_MC_err',col_MC_err);        % Monte Carlo error band color
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
% L. Baars       | 2023-Apr-20 | Adjusted default significant digits
% L. Baars       | 2024-Jan-11 | Rename of PcConjPlanCircle to PcCircle
% L. Baars       | 2025-Aug-25 | Updated code for public release. Added
%                                scaleFactor parameter.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================