function yrng = ConjDist_time_panel(n,xrng,xfmt,ymult,bins, ...
    ybin,ybinavg,yavg,plt,params)
% ConjDist_time_panel - Plot a panel of the Conj_Dist time plots
%
% Syntax: yrng = ConjDist_time_panel(n,xrng,xfmt,ymult,bins, ...
%                                    ybin,ybinavg,yavg,plt,params)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Plot a panel of the Conj_Dist time plots
%
% =========================================================================
%
% Input:
%
%   n       - index to plot (1-3)
%
%   xrng    - x-axis plot range                                       [1x2]
%
%   xfmt    - Date format for x-axis ticks expressed as a string 
%             recognizable by datetime. May be set as empty to 
%             plot in days elapsed
%
%   ymult   - y-value multiplier
%
%   bins    - x-axis time-bin boundaries                            3x[1xN]
%
%   ybins   - struct containing median and +/-1 sigma y-bins        3x[3xN]
%
%   ybinavg - struct containing average median and +/-1 sigma       3x[3xN]
%             y-bins
%
%   yavg    - struct containing average median and +/-1 sigma       3x[3xN]
%             y-data
%
%   plt     - plot parameter struct (See MATLAB documentation for 
%             plot function)
%
%   params  - EventRate parameters structure. See 
%             EventRate_default_params.m and/or
%             EventRate_ConjDist_default_params.m for documentation
% =========================================================================
%
% Output:
%
%   yrng - y-axis plot range                                          [1x2]
%
% =========================================================================
%
% Dependencies:
%
%   None
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

if ~isfield(plt,'wid'); plt.wid = []; end
if isempty(plt.wid); plt.wid = 1; end

if isempty(xfmt)
    xoff = min(bins.Tbin);
else
    xoff = 0;
end
xplt = bins.Tbin-xoff;
xrng = xrng-xoff;

ybin_md = ymult*ybin.md;
ybin_lo = ymult*ybin.lo;
ybin_hi = ymult*ybin.hi;

ybinavg_md = ymult*ybinavg.md;
ybinavg_lo = ymult*ybinavg.lo;
ybinavg_hi = ymult*ybinavg.hi;

yavg_md = ymult*yavg.md;
yavg_lo = ymult*yavg.lo;
yavg_hi = ymult*yavg.hi;

plot(xplt,ybin_md(n,:), ...
    'LineStyle',plt.sty, ...
    'LineWidth',plt.wid, ...
    'Color',plt.col, ...
    'Marker',plt.sym, ...
    'MarkerSize',plt.siz, ...
    'MarkerFaceColor',plt.col, ...
    'MarkerEdgeColor',plt.col);

ymn0 = min(ybin_md(n,:));
ymx0 = max(ybin_md(n,:));

if any(params.plot_err_bars)
    if params.plot_err_bars(1)
        xplo = bins.Tbinlo-xoff;
        xphi = bins.Tbinhi-xoff;
    end
    hold on;
    M = numel(ybin_md(n,:));
    for m=1:M
        if params.plot_err_bars(1) % Horizontal error bars
            plot([xplo(m) xphi(m)],[ybin_md(n,m) ybin_md(n,m)], ...
                'LineStyle','-', ...
                'LineWidth',1, ...
                'Color',plt.col, ...
                'Marker','none');
            ymn0 = min(ymn0,ybin_md(n,m));
            ymx0 = max(ymx0,ybin_md(n,m));
        end
        if params.plot_err_bars(2) % Vertical error bars
            plot([xplt(m) xplt(m)],[ybin_lo(n,m) ybin_hi(n,m)], ...
                'LineStyle','-', ...
                'LineWidth',1, ...
                'Color',plt.col, ...
                'Marker','none');
            ymn0 = min(ymn0,ybin_lo(n,m));
            ymx0 = max(ymx0,ybin_hi(n,m));
        end
    end
    hold off;
end

if any(params.plot_conf_ranges)
    hold on;
    plot(xrng,[yavg_md(n) yavg_md(n)], ...
        'LineStyle','-', ...
        'LineWidth',1, ...
        'Color',plt.col, ...
        'Marker','none');
    ymn0 = min(ymn0,yavg_md(n));
    ymx0 = max(ymx0,yavg_md(n));
    if params.plot_conf_ranges(1) % Unified bin-set confidence range
        plot(xrng,[yavg_lo(n) yavg_lo(n)], ...
            'LineStyle','-.', ...
            'LineWidth',1, ...
            'Color',plt.col, ...
            'Marker','none');
        plot(xrng,[yavg_hi(n) yavg_hi(n)], ...
            'LineStyle','-.', ...
            'LineWidth',1, ...
            'Color',plt.col, ...
            'Marker','none');
        ymn0 = min(ymn0,yavg_lo(n));
        ymx0 = max(ymx0,yavg_hi(n));
    end
    if params.plot_conf_ranges(2) % Bin-to-bin variation confidence range
        plot(xrng,[ybinavg_lo(n) ybinavg_lo(n)], ...
            'LineStyle','--', ...
            'LineWidth',1, ...
            'Color',plt.col, ...
            'Marker','none');
        plot(xrng,[ybinavg_hi(n) ybinavg_hi(n)], ...
            'LineStyle','--', ...
            'LineWidth',1, ...
            'Color',plt.col, ...
            'Marker','none');
        ymn0 = min(ymn0,ybinavg_lo(n));
        ymx0 = max(ymx0,ybinavg_hi(n));
    end
    hold off;
end

if ~isempty(xfmt)
    set(gca,'XTick',[min(xrng) mean(xrng) max(xrng)]);
    datetick('x',xfmt,'keeplimits','keepticks');
    xlim(xrng);
else
    xlim(xrng);
end

if (ymn0 == 0) && (ymx0 == 0)
    yrng = [-0.05 1.05];
else
    yrng = plot_range([ymn0 ymx0],0.05);
end
ylim(yrng);

[md,lo,hi] = smart_error_range(yavg_md(n),yavg_lo(n),yavg_hi(n));
titl{1} = ['Mean ' md ' (95% conf ' lo ' to ' hi ')'];
[~,lo,hi] = smart_error_range(ybinavg_md(n),ybinavg_lo(n),ybinavg_hi(n));
titl{2} = ['Bin variations ' lo ' to ' hi];
title(titl,'FontWeight','bold','FontSize',9);
set(gca,'FontSize',9,'FontWeight','normal','LineWidth',2);

drawnow;

return;
end
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================