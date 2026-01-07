function h = ConjDist_time_plot(nset,titl,bins,rates,fracs,params,h)
% ConjDist_time_panel - Make a 5-panel temporal plot of the OCMDB rates and 
%                       fractions
%
% Syntax: h = ConjDist_time_plot(nset,titl,bins,rates,fracs,params,h)
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
%   nset   - set index to process (equal to 1, 2, or 3)
%
%   titl   - cell array of plot titles                                {Nx1}
%
%   bins   - struct containing time bins
%
%   rates  - struct containing yellow/red event rates for 
%            SVI/non-SVI/All cases 
%
%   fracs  - struct containing yellow/red event probabilities for
%            SVI/non-SVI/All cases
%
%   params - EventRate parameters structure. See EventRate_default_params.m 
%            and/or EventRate_ConjDist_default_params.m for documentation
%
%   h      - Figure handle (optional, generated automatically if not
%            provided)
% =========================================================================
%
% Output:
%
%   h      - figure handle
%
% =========================================================================
%
% Dependencies:
%
%   ConjDist_time_panel.m
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

Nargin = nargin;
if (Nargin < 7); h = []; end
if ishandle(h)
    set(0, 'CurrentFigure', h);
else
    if params.visible_plots
        vis = 'On';
    else
        vis = 'Off';
    end
    h = figure('Visible',vis);
end

clf;

lsty = ':';
msiz = 4;

plt.all.col = params.acol; plt.all.sym = 's'; plt.all.siz = msiz-1; plt.all.sty = lsty;
plt.yel.col = params.ycol; plt.yel.sym = 'o'; plt.yel.siz = msiz;   plt.yel.sty = lsty;
plt.red.col = params.rcol; plt.red.sym = 'd'; plt.red.siz = msiz;   plt.red.sty = lsty;

xrng = plot_range([min(bins.Tbinlo) max(bins.Tbinhi)],0.05);

xlabl = 'Time';

Tspn = diff(xrng);

if (Tspn <= 100)
    xfmt = [];
    xlabl = [xlabl ' (days)'];
else
    xfmt = 'yyyymmmdd';
end

% if (Tspn <= 100)
%     xfmt = [];
%     xlabl = [xlabl ' (days)'];
% elseif (Tspn <= 200)
%     xfmt = 'ddmmmyy';
% elseif (Tspn <= 500)
%     xfmt = 'mmmyy';
% else
%     xfmt = 'yyyy';
% end

ylabl = 'Rate (day^{-1})';
ymult = 1;

subplot(3,2,1);
ConjDist_time_panel( ...
    nset,xrng,xfmt,ymult,bins, ...
    rates.binall,rates.bavgall,rates.all,plt.all,params);
ylabel(ylabl);

subplot(3,2,3);
ConjDist_time_panel( ...
    nset,xrng,xfmt,ymult,bins, ...
    rates.binyel,rates.bavgyel,rates.yel,plt.yel,params);
ylabel(ylabl);

subplot(3,2,5);
ConjDist_time_panel( ...
    nset,xrng,xfmt,ymult,bins, ...
    rates.binred,rates.bavgred,rates.red,plt.red,params);
ylabel(ylabl);

if ~isempty(xlabl)
    xlabel(xlabl);
end

ylabl = 'Fraction (%)';
ymult = 100;

subplot(3,2,4);
ConjDist_time_panel( ...
    nset,xrng,xfmt,ymult,bins, ...
    fracs.binyel,fracs.bavgyel,fracs.yel,plt.yel,params);
ylabel(ylabl);

subplot(3,2,6);
ConjDist_time_panel( ...
    nset,xrng,xfmt,ymult,bins, ...
    fracs.binred,fracs.bavgred,fracs.red,plt.red,params);
ylabel(ylabl);

if ~isempty(xlabl)
    xlabel(xlabl);
end

subplot(3,2,2);
plot(NaN,NaN);
axis off;
Ntitl = numel(titl);
dy = 1/max(6,Ntitl);
y = 1-dy*((1:Ntitl)-0.5);
for n=1:Ntitl
    text(0.5,y(n),titl{n}, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold');
end

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