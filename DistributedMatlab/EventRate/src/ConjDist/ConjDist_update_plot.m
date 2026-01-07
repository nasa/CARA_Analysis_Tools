function h = ConjDist_update_plot(titl,TCAlud,TCRseq,Pcseq,params,h)
% ConjDist_update_plot - Make a temporal update-sequence plot
%
% Syntax: h = ConjDist_update_plot(titl,TCAlud,TCRseq,Pcseq,params,h)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Make a temporal update-sequence plot
%
% =========================================================================
%
% Input:
%
%   titl   - Plot title                                         
%
%   TCAlud - TCA at last update
%
%   TCRseq - Times of OCMDB entry creation                           [1xN]
%
%   Pcseq  - Calculated Pcs                                          [1xN]
%
%   params - EventRate parameters structure. See EventRate_default_params.m 
%            and/or EventRate_ConjDist_default_params.m for documentation
%
%   h      - Optional - Figure handle 
%            (generated automatically if not provided)
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

% Initializations

Nargin = nargin;

if (Nargin < 4); h = []; end

% Set up figure

if ishandle(h)
    set(0, 'CurrentFigure', h)
else
    if params.visible_plots
        vis = 'On';
    else
        vis = 'Off';
    end
    h = figure('Visible',vis);
end
clf;

xfsz = 12;
yfsz = xfsz;
xfwt = 'bold';
yfwt = xfwt;
tfsz = xfsz;
tfwt = 'bold';
afsz = xfsz;
afwt = 'bold';
alwd = 1;

% Sort data
[TCRseq,srt] = sort(TCRseq);
Pcseq = Pcseq(srt);

% Eliminate points that are too late to plot
ndx = TCRseq-TCAlud > params.postTCA_limit_days;
if any(ndx)
   ndx = ~ndx;
   TCRseq = TCRseq(ndx);
   Pcseq = Pcseq(ndx);
end

% Last update before commit time Pc option
if strcmpi(params.redyel_event_option,'lud')
    % Last update Pc value between consider and commit times, 
    % used for modeling RMM execution rates
    PcLUDvsPcMaxOption = 1;
elseif strcmpi(params.redyel_event_option,'mud')
    % Maximum Pc value between consider and commit times, 
    % used for modeling RMM planning rates
    PcLUDvsPcMaxOption = 2;
else
    % Non standard operational mode
    PcLUDvsPcMaxOption = 0;
end

% Plot data

xplt = TCRseq-TCAlud;
xlabl = 'Time from TCA (days)';
xmn0 = min(xplt);
xmx0 = max(0,max(xplt));
if PcLUDvsPcMaxOption > 0
    xmn0 = min(xmn0,min(-params.CommitConsiderDays));
end    
xrng = plot_range([xmn0 xmx0],0.05);

yplt = Pcseq;
yplt(Pcseq < params.Pcgre) = params.Pcgre;
yplt = log10(yplt);
ylabl = 'Log_{10}(Pc)';
yrng = plot_range( ...
    [min(yplt) max(yplt) params.log10_Pcgre params.log10_Pcred],0.05);

% Plot green, yellow and red horizontal lines

lsty = '-';
lwid = 2;

plot(xrng,[params.log10_Pcgre params.log10_Pcgre], ...
    'LineStyle',lsty, ...
    'LineWidth',lwid, ...
    'Color',params.gcol);
 
hold on;

plot(xrng,[params.log10_Pcyel params.log10_Pcyel], ...
    'LineStyle',lsty, ...
    'LineWidth',lwid, ...
    'Color',params.ycol);

plot(xrng,[params.log10_Pcred params.log10_Pcred], ...
    'LineStyle',lsty, ...
    'LineWidth',lwid, ...
    'Color',params.rcol);

% Plot a vertical line at the remediation maneuver commit time

if PcLUDvsPcMaxOption > 0
    lsty = ':';
    lcol = 'k';
    for nn=1:numel(params.CommitConsiderDays)
        plot(-[params.CommitConsiderDays(nn) params.CommitConsiderDays(nn)],yrng, ...
            'LineStyle',lsty, ...
            'LineWidth',lwid, ...
            'Color',lcol);
        % if numel(titl) >= 3
        %     titl{3} = ['\color{magenta}' titl{3}];
        % end
    end
end    

% Plot a line connecting data points
 
lsty = '-';
lwid = 1;
lcol = 'k';
mrkr = 'none';
mcol = lcol;
msiz = 6;

plot(xplt,yplt, ...
    'LineStyle',lsty, ...
    'LineWidth',lwid, ...
    'Color',lcol, ...
    'Marker',mrkr, ...
    'MarkerSize',msiz, ...
    'MarkerFaceColor',mcol, ...
    'MarkerEdgeColor',mcol);

% Plot the symbols

lsty = 'none';
lwid = 1;
msiz0 = msiz;

Nplt = numel(xplt);

for n=1:Nplt
    
    if Pcseq(n) <= 0
        mrkr = 'o';
        msiz = msiz0;
    elseif Pcseq(n) <= params.Pcgre
        mrkr = 'v';
        msiz = msiz0-1;
    else
        mrkr = 'd';
        msiz = msiz0;
    end
    
    plot(xplt(n),yplt(n), ...
        'LineStyle',lsty, ...
        'LineWidth',lwid, ...
        'Color',lcol, ...
        'Marker',mrkr, ...
        'MarkerSize',msiz, ...
        'MarkerFaceColor',mcol, ...
        'MarkerEdgeColor',mcol);       

end

hold off;

xlim(xrng);
ylim(yrng);

xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
title(titl,'FontSize',tfsz,'FontWeight',tfwt);

set(gca,'FontSize',afsz,'FontWeight',afwt);
set(gca,'LineWidth',alwd);

return
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