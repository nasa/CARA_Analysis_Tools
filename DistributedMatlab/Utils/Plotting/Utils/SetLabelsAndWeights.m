function SetLabelsAndWeights(xlabl, ylabl, params)
% SetLabelsAndWeights - Standardized function to set the labels and weights
%                       for the current figure axes.
%
% Syntax: SetLabelsAndWeights(xlabl, ylabl, params);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%  xlabl  - Label for the x-axis
%  ylabl  - Label for the y-axis
%  params - Auxilliary input parameter structure with the following fields:
%             fig.xfsz, fig.xfwt, fig.yfsz, fig.yfwt, fig.afsz, fig.afwt,
%             and fig.alwd
%           See the documentation within default_params_Pc_Time_Plot.m
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    xlabel(xlabl,'FontSize',params.fig.xfsz,'FontWeight',params.fig.xfwt);
    ylabel(ylabl,'FontSize',params.fig.yfsz,'FontWeight',params.fig.yfwt);
    set(gca,'FontSize',params.fig.afsz,'FontWeight',params.fig.afwt);
    set(gca,'LineWidth',params.fig.alwd);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-23 | Initial Development
% L. Baars       | 2025-Aug-25 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
