function [pltEnabled, mrkr, mcol, msiz, lcol, lnst, lnwd, name] = GetPlotFormat(plotType,params)
% GetPlotFormat - Gets plot format information for the plot type passed in.
%
% Syntax: [pltEnabled, mrkr, mcol, msiz, lcol, lnst, lnwd, name] = GetPlotFormat(plotType,params)
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
%   This function retrieves plot format information based on the plotType
%   and parameters structure passed in. This function is meant to provide a
%   common format for any probability of collision type plots made by CARA.
%
% =========================================================================
%
% Input:
%
%   plotType - Character string defining the type of plot. Valid values
%              are: 'Pc_SDMC', 'PcCircle', 'Pc2D_Hall',
%              'Pc2D_Hall_ConjBound', 'Pc3D_Hall', or
%              'Pc3D_Hall_ConjBound'.
%
%   params - Auxilliary input parameter structure with the following
%            fields:
%              plot_PcCircle, plot_Pc2D_Hall, plot_Pc3D_Hall,
%              plot_Pc_SDMC, trajectory_mode, fig.msiz, fig.lnwd1,
%              fig.lnwd2, fig.dark_green, and fig.col_MC
%            See the documentation within default_params_Pc_Time_Plot.m
%
% =========================================================================
%
% Output:
%
%   pltEnabled - Boolean indicating if the plot is enabled
%
%   mrkr - Marker indicator or 'none'
%
%   mcol - Marker color
%
%   msiz - Marker size
%
%   lcol - Line color
%
%   lnst - Line style
%
%   lnwd - Line width
%
%   name - Display name to use in legends
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % Default values
    mrkr = 'none';
    msiz = params.fig.msiz;
    lnst = '-';
    lnwd = params.fig.lnwd1;
    % Conditionals for each plot type
    if strcmp(plotType,'Pc_SDMC')
        if params.plot_Pc_SDMC
            pltEnabled = true;
        else
            pltEnabled = false;
        end
        mcol = params.fig.col_MC;
        lnwd = params.fig.lnwd2;
        % name = 'SDMC Pc';
        name = 'SDMC Method Pc';
        if params.trajectory_mode == 1
            name = [name ' (full rect.)'];
        elseif params.trajectory_mode == 2
            name = [name ' (rect.)'];
        end
    elseif strcmp(plotType,'PcCircle')
        if params.plot_PcCircle
            pltEnabled = true;
        else
            pltEnabled = false;
        end
        mcol = 'b';
        % name = '2D Pc Conj. Plane';
        name = '2D-Pc Method Pc';
    elseif strcmp(plotType,'Pc2D_Hall')
        if params.plot_Pc2D_Hall
            pltEnabled = true;
        else
            pltEnabled = false;
        end
        mcol = 'c';
        lnst = '--';
        % name = '2D Pc Unit Sphere';
        name = '2D-Nc Method Pc';
    elseif strcmp(plotType,'Pc2D_Hall_ConjBound')
        if params.plot_Pc2D_Hall
            pltEnabled = true;
        else
            pltEnabled = false;
        end
        mcol = 'g';
        lnwd = params.fig.lnwd2;
        % name = '2D Pc Unit Sphere Conj Bound (gamma = 1e-6)';
        name = '2D-Nc Duration Bounds (\gamma = 1e-6)';
    elseif strcmp(plotType,'Pc3D_Hall')
        if params.plot_Pc3D_Hall
            pltEnabled = true;
        else
            pltEnabled = false;
        end
        mcol = 'k';
        lnst = ':';
        % name = '3D Pc Unit Sphere';
        name = '3D-Nc Method Pc';
    elseif strcmp(plotType,'Pc3D_Hall_ConjBound')
        if params.plot_Pc3D_Hall
            pltEnabled = true;
        else
            pltEnabled = false;
        end
        mcol = params.fig.darkgreen;
        % name = '3D Pc Unit Sphere Conj Bound (gamma = 1e-6)';
        name = '3D-Nc Duration Bounds (\gamma = 1e-6)';
    else
        error(['Unsupported plotType passed in: ' plotType]);
    end
    lcol = mcol;
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
% L. Baars       | 2023-Apr-19 | Display name change for PcConjPlaneCircle
% D. Hall        | 2023-Aug-07 | Display name change for all methods
% L. Baars       | 2024-Jan-11 | Renamed PcConjPlaneCircle to PcCircle
% L. Baars       | 2025-Aug-25 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
