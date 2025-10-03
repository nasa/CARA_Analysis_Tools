function [params] = default_params_CA_Dist_Plot(params)
% default_params_CA_Dist_Plot - Add and/or set the defaults for the
%                               parameters used by the function
%                               CA_Dist_Plot.
%
% Syntax: params = default_params_CA_Dist_Plot;
%         params = default_params_CA_Dist_Plot(params);
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
%          parameters used by the CA_Dist_Plot function, containing the
%          following fields:
%
%   generate_ca_dist_plot = Boolean which controls whether or not the close
%                           approach distribution plot is generated.
%                           Defaults to false.
%
%   alt_ca_dist = Allows the selection of alternative CA distribution
%                 plots:
%                   0 = Nominal CA distribution
%                   1 = Combined covariance placed at the origin
%                   2 = Plane is rotated to align the axes with the
%                       principal axes of the 2D-Pc ellipse.
%                 Defaults to 0.
%
%   AuxCAdistContour - Provides an option to turn on a beta-version
%                      analytical visualization of the CA distribution
%                      contour.
%                        0 = Never plot auxiliary contours for the 3D-Nc or
%                            2D-Nc methods
%                        1 = Plot 3D-Nc contour, but only for 2D-Pc usage
%                            violations
%                        2 = Always plot 3D-Nc contour (the "bananoid"),
%                            when possible (for testing and R&D)
%                        3 = Always plot 2D-Nc contour, when possible (for
%                            testing and R&D)
%                      Note: As of the current release, options 1-3 are
%                            disabled in the public-release version of
%                            the software.
%
%   fig = Sub-structure used to define figure specific parameters,
%         including the following fields:
%
%     afsz = Axis font size (default = 14)
%
%     tfsz = graph title font size (default = params.fig.afsz)
%
%     lfsz = legend font size (default = 12)
%
%     gray = the color gray (default = [0.1 0.1 0.1])
%
%   Note: This function also calls default_params_Pc_Time_Plot to populate
%         other default parameters. See the documentation from that
%         function for a list of other parameters which will appear in the
%         params structure.
%
% =========================================================================
%
% OUTPUT:
%
%   params = Fully-populated structure containing parameters used by the 
%            function CA_Dist_Plot.
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
    params = set_default_param(params,'generate_ca_dist_plot',false);
    
    % Use the alternative CA distribution coordinate system
    params = set_default_param(params,'alt_ca_dist',0);
    if min(abs(params.alt_ca_dist-[0 1 2])) > 0
        error('Invalid alt_ca_dist parameter');
    end
    
    % Plot the auxiliary (bananaoid) contour 
    params = set_default_param(params,'AuxCAdistContour',0);
    
    % Figure parameters which override or add to the parameters from
    % default_params_Pc_Time_Plot
    params = set_default_param(params,'fig',[]);
    params.fig = set_default_param(params.fig,'afsz',14);                  % axis font size
    params.fig = set_default_param(params.fig,'tfsz',params.fig.afsz);     % title font size
    params.fig = set_default_param(params.fig,'lfsz',12);                  % legend font size
    params.fig = set_default_param(params.fig,'gray',[0.1 0.1 0.1]);
    
    % Set the defaults from Pc_Time_Plot
    params = default_params_Pc_Time_Plot(params);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-20 | Initial development.
% D. Hall        | 2023-Nov-11 | Modified params.alt_ca_dist defaults.
% D. Hall        | 2023-Nov-18 | Added params.AuxCAdistContour defaults.
% L. Baars       | 2025-Aug-25 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================