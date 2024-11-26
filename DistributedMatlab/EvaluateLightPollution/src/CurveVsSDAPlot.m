function hlegnd = CurveVsSDAPlot(xplt,yplt,nplt,Nplt,Ncurve,Ncomp,yrng,params)
% CurveVsSDAPlot - Curve versus Solar Depression Angle (SDA) plot
% Syntax: hlegnd = CurveVsSDAPlot(xplt,yplt,nplt,Nplt,Ncurve,Ncomp,yrng,params)
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: This function generates a plot of a curve versus SDA for 
% visualizing metrics such as number of satellites sunlit (Ni) or brighter 
% than threshold (Nb). 
%
% If this is the Nb plot, then plot the shaded evaluation regions of 
% varying risk levels.
%
% =========================================================================
%
% Input:
%
%    xplt       -   X-axis SDA data points
%
%    yplt       -   Y-axis data points (e.g., Ni or Nb)
%                   
%    nplt       -   Current plot identifier (1. Ni vs SDA or 2. Nb vs SDA)
%    
%    Nplt       -   Total number of plots (2)
%    
%    Ncurve     -   Curve identifier
%
%    Ncomp      -   Number of comparison datasets
%
%    yrng       -   Y-axis range
%
%    params     -   Auxilliary input parameter structure
%
% =========================================================================
%
% Output:
%
%   Plot        -   Ni vs SDA or Nb vs SDA Plot
% 
%   hlegnd      -   Legend handle
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

render_shaded_areas = (nplt == Nplt) & (Ncurve == 1);
Ncurve = abs(Ncurve);

if render_shaded_areas

    LoClr = [1.0 1.0 0.8];
    MdClr = [1.0 1.0 0.0];
    HiClr = [1.0 0.0 0.0];
    VHClr = [1.0 0.0 1.0];

    plot(NaN,NaN);
    hold on;

    if params.Plotting.NbrightShadingLevel > 0
        NSDAPoints = numel(params.Evaluation.SDAPoints);
        for n=1:NSDAPoints
            Mb = params.Evaluation.MaxBright(n);
            Mb(Mb < 0) = NaN;
            if ~isnan(Mb)
                SDA1 = params.Evaluation.SDAPoints(n);
                if n < NSDAPoints
                    SDA2 = params.Evaluation.SDAPoints(n+1);
                else
                    SDA2 = 90;
                end
                xfil = [SDA1 SDA2 SDA2 SDA1];
                if params.Plotting.NbrightShadingLevel > 2
                    % Fill the low risk color
                    Nb1 = Mb*params.Evaluation.LowToHigh(1);
                    Nb2 = Mb*params.Evaluation.LowToHigh(2);
                    yfil = [Nb1 Nb1 Nb2 Nb2];
                    fill(xfil,yfil,LoClr,'EdgeColor',LoClr);
                end
                if params.Plotting.NbrightShadingLevel > 1
                    % Fill the medium risk color
                    Nb1 = Mb*params.Evaluation.LowToHigh(2);
                    Nb2 = Mb*params.Evaluation.LowToHigh(3);
                    yfil = [Nb1 Nb1 Nb2 Nb2];
                    fill(xfil,yfil,MdClr,'EdgeColor',MdClr);
                end
                % Fill the high risk color
                Nb1 = Mb*params.Evaluation.LowToHigh(3);
                Nb2 = Mb*params.Evaluation.LowToHigh(4);
                yfil = [Nb1 Nb1 Nb2 Nb2];
                fill(xfil,yfil,HiClr,'EdgeColor',HiClr);
                % Fill the very high risk color
                Nb1 = Mb*params.Evaluation.LowToHigh(4);
                Nb2 = yrng(2);
                yfil = [Nb1 Nb1 Nb2 Nb2];
                fill(xfil,yfil,VHClr,'EdgeColor',VHClr);
            end
        end
    end

end

% Add an interpolated point to complete the curve to bottom axis,
% if required
ndx = ~(yplt == 0);
if min(yplt(ndx)) > yrng(1)
    Nyplt = numel(yplt);
    ndx = find(ndx); ndx = ndx(end);
    xm = xplt(ndx); ym = yplt(ndx);
    if ndx < Nyplt
        xp = xplt(ndx+1);
    else
        xp = xplt(ndx)+(xplt(ndx)-xplt(ndx-1));
    end
    yp = 0;
    xint = interp1([ym yp],[xm xp],yrng(1));
    xplt = cat(1,xplt,xint);
    yplt = cat(1,yplt,yrng(1));
    [xplt,nsrt] = sort(xplt);
    yplt = yplt(nsrt);
end

% Determine the line color and line style
colors = {[1 1 1]/2 [0 1 1] [0 1 0] [0 0 1]/2 [0 1 0]/2 [0 1 1]/2 };
Ncolors = numel(colors);
styles = {'-.' '--' ':'};
Nstyles = numel(styles);

if Ncurve == 1
    lsty = '-'; lcol = [0 0 0]; lwid = 4;
else
    ncol = mod(Ncurve-2,Ncolors)+1; nsty = mod(Ncurve-2,Nstyles)+1;    
    lsty = styles{nsty}; lcol = colors{ncol}; lwid = 2;
end

if Ncomp > 0
    lwid = 2;
end

% Plot the curve, capturing the handle for output
hlegnd = plot(xplt,yplt, ...
    'Color',lcol, ...
    'LineStyle',lsty, ...
    'LineWidth',lwid);

return;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================