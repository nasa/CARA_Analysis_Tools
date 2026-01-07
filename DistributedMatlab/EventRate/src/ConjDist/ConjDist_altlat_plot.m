function [xrng,yrng,zrng,h] = ConjDist_altlat_plot( ...
    plotmode,PcLUDvsPcMaxOption,eset,titl,xrng,yrng,zrng,CDout,params,h)
% ConjDist_altlat_plot - Plot the altitude-latitude distribution for a set
%                        of events, or update-sequences of conjunctions
% Syntax: [xrng,yrng,zrng,h] = ...
%         ConjDist_altlat_plot( plotmode,PcLUDvsPcMaxOption,eset,titl, ...
%                               xrng,yrng,zrng,CDout,params,h)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Plot the altitude-latitude distribution for a set of events, 
% or update-sequences of conjunctions
%
% =========================================================================
%
% Input:
%
%   plotmode           - Integer value ranging from 0-3
%
%   PcLUDvsPcMaxOption - Pc-max vs Pc-latest-update flag
%
%   eset               - Event set indices to plot
%
%   titl               - Plot title (string)
%
%   xrng               - Plotting range for x-values                  [1x2]
%                        Will be automatically generated if empty
%
%   xrng               - Plotting range for y-values                  [1x2]
%                        Will be automatically generated if empty
%
%   xrng               - Plotting range for z-values                  [1x2]
%                        Will be automatically generated if empty
%
%   CDout              - Conjunction data struct with the following fields:
%
%       TCAmed: median TCA for each event                        [neventx1]
%
%       DB: OCM DB table                                        [nconjx277]
%
%       event: array of event numbers corresponding to each       [nconjx1]
%       conjunction            
%                       
%       Pclud: Pc at last update for each event                  [neventx1]
%       (only used if lud option is selected)
%
%       Pcmud: Pc at max update for each event                   [neventx1]
%       (only used if lud option is selected)
%
%       Pcmax: Maximum Pc for each event                         [neventx1]
%       (only used if neither option is selected)
%
%       USElud: array of bools indicating if lud was calculated  [neventx1]
%       for each event. (only used if neither option is selected)
%
%       USEmud: array of bools indicating if mud was calculated  [neventx1]
%       for each event. (only used if neither option is selected)
%
%   params             - EventRate parameters structure. See
%                        EventRate_default_params.m and/or 
%                        EventRate_ConjDist_default_params.m for
%                        documentation
%
%   h                  - Figure handle to plot to. Will be set
%                        automatically if unset
%
% =========================================================================
%
% Output:
%
%   xrng               - Plotting range for x-values                  [1x2]
%
%   yrng               - Plotting range for y-values                  [1x2]
%
%   zrng               - Plotting range for z-values                  [1x2]
%
%   h                  - Figure handle
%
% =========================================================================
%
% Dependencies:
%
%   None
%
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
if (Nargin < 9); h = []; end

Re = 6378.137;

% Check for valid modes

if min(abs(plotmode-(0:3))) ~= 0
    error('Invalid altlat plot mode');
end

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

xfsz = 11;
xfwt = 'bold';
yfsz = xfsz;
yfwt = xfwt;
zfsz = xfsz;
zfwt = xfwt;
tfsz = xfsz-1;
tfwt = 'bold';
tfan = 'italic';
afsz = xfsz;
afwt = 'bold';
alwd = 1;
cfsz = tfsz-2;
cfwt = 'bold';
lfsz = cfsz;
lfwt = 'bold';

% Plot the altitude-latitiude distribution

subplot(3,20,[(28:40) (48:60)]);

% Process plotting modes:
%  plotmode = 0  =>  gray histogram, black diamonds, no   colorbar
%  plotmode = 1  =>  gyr  histogram, gyr   diamonds, no   colorbar
%  plotmode = 2  =>  gray histogram, cbar  diamonds, time colorbar
%  plotmode = 3  =>  gray histogram, cbar  bubbles  , time colorbar

% mrkr = 'x'; % Default marker
% msiz = 5;   % Default marker size

mrkr = 'd'; % Default marker
msiz = 3;   % Default marker size

if plotmode == 3
    % Pc-sized bubble markers
    msiz_Pcyel = 3;
    msiz_Pcred = 11;
    msiz_delta = (msiz_Pcred-msiz_Pcyel) ...
        / (params.log10_Pcred-params.log10_Pcyel);
end

altlat_time_cbar = plotmode > 1; % Colorbar for modes 2 and 3

% Define the z range for the colorbar

zlabl = [];
if altlat_time_cbar
    % Z axis for colorbar is time in days
    if isempty(zrng)
        zmin = min(CDout.TCAmed(eset));
        zmax = max(CDout.TCAmed(eset));
        zrng = [zmin zmax];
    end
    zdel = diff(zrng);
    if (zdel <= 3)
        zlabl = ['Time from ' datestr(zrng(1),'yymmm')];
        datefmt = 'dd HH:MM';
    else
        datefmt = 'yyyymmmdd';
    end
    % Colormap
    Nmap = 256;
    % Parula colormap
    cmap = colormap(parula(256));
    Nmapm1 = Nmap-1;
end

% Calculate latitude range
using_PcTable = ~isempty(params.PcTable);
if using_PcTable

    % Get altitudes & latitudes stored during PcTable construction
    B = CDout.DB(:,params.PcTable.ExpandedDBColumns);

else

    % ECI coordinates for the conjunctions in this event
    X = CDout.DB(:,172);
    Y = CDout.DB(:,173);
    Z = CDout.DB(:,174);

    % Radius (R), altitude (A), and latitude (B).
    R = sqrt(X.^2+Y.^2+Z.^2);
    % A = R-Re;
    B = asin(Z./R)*180/pi;
    B(B < -90) = -90; B(B > +90) = +90;

end

% Median values
Bmax = max(abs(B));

% Bins for the latitude distribution

NbinTarget = 24;
Bwd = 2*Bmax/NbinTarget;
Nb = max(2,round(180/Bwd));
    
bwd = 180/Nb;
bwh = bwd/2;
bhi = -90+bwd*(1:Nb);
blo = bhi-bwd;
bmd = bhi-bwh;
bct = zeros(size(bmd));
bcy = bct; bcr = bct;

% Find all of the conjunctions for the current event set

evset  = find(eset);
Nevset = numel(evset);

% Sort by date

[~,evsrt] = sort(CDout.TCAmed(evset));

% Pc-max vs Pc-latest-update

if PcLUDvsPcMaxOption == 1
    Pcplt = CDout.Pclud;
    Pcuse = CDout.USElud;
elseif PcLUDvsPcMaxOption == 2
    Pcplt = CDout.Pcmud;
    Pcuse = CDout.USEmud;
else
    Pcplt = CDout.Pcmax;
    Pcuse = true(size(Pcplt));
end

% Plot one update-sequence event per loop in sorted order

initializing_plot = true;
        
for nnne=1:Nevset
    
    % Only plot usable events
    
    if Pcuse(nnne)

        % Indices for this update-sequence event
        nne = evsrt(nnne);
        ne = evset(nne);
        ndx = find(ne == CDout.event);
        
        % MarkerFaceColor
        if altlat_time_cbar % plot modes 2 or 3 have time colorbar
            % MarkerFaceColor from colorbar
            znow = CDout.TCAmed(ne);
            nmap = 1+round(Nmapm1*(znow-zrng(1))/diff(zrng));
            nmap = max(1,min(Nmap,nmap));
            mcol = cmap(nmap,:);
        else
            if plotmode == 0
                % MarkerFaceColor black for mode 0
                mcol = [0 0 0];
            else
                % MarkerFaceColor red, yellow or green for mode 1
                if Pcplt(ne) >= params.Pcred
                    mcol = params.rcol;
                elseif Pcplt(ne) >= params.Pcyel
                    mcol = params.ycol;
                else
                    mcol = params.gcol;
                end
            end
        end

        if using_PcTable

            % Get altitudes & latitudes stored during PcTable construction
            A = CDout.DB(ndx,params.PcTable.ExpandedDBColumns-1)/1e3;
            B = CDout.DB(ndx,params.PcTable.ExpandedDBColumns);

        else

            % ECI coordinates for the conjunctions in this event
            X = CDout.DB(ndx,172);
            Y = CDout.DB(ndx,173);
            Z = CDout.DB(ndx,174);

            % Radius (R), altitude (A), and latitude (B).
            R = sqrt(X.^2+Y.^2+Z.^2);
            A = R-Re;
            B = asin(Z./R)*180/pi;
            B(B < -90) = -90; B(B > +90) = +90;
            
        end

        % Median values
        Amed = median(A);
        Bmed = median(B);

        % Latitude bin
        ibin = find((blo <= Bmed) & (Bmed <= bhi));
        ibin = ibin(1);
        bct(ibin) = bct(ibin)+1;

        if (Pcplt(ne) >= params.Pcyel)
            bcy(ibin) = bcy(ibin)+1;
        end

        if (Pcplt(ne) >= params.Pcred)
            bcr(ibin) = bcr(ibin)+1;
        end

        % Add this event to the alt-lat dist plot

        if plotmode == 3 % plot mode 3
            if initializing_plot
                % Use scatter function to initiate colorbar, but otherwise use
                % plot function for Pc-sized bubbles
                scatter(NaN,NaN,msiz,mcol,'fill');
                hold on;
            end
            Pcsiz = Pcplt(ne);
            if Pcsiz < params.Pcyel; Pcsiz = params.Pcyel; end
            if Pcsiz > params.Pcred; Pcsiz = params.Pcred; end
            msiz = msiz_Pcyel + (log10(Pcsiz)-params.log10_Pcyel)*msiz_delta;
            % disp(['log10Pc = ' num2str(log10(Pcsiz)) '  msiz = ' num2str(msiz)]);
            ecol = [0 0 0]; % MarkerEdgeColor black
            plot(Amed,Bmed, ...
                'Marker','o', ...
                'MarkerSize',msiz, ...
                'MarkerFaceColor',mcol, ...
                'MarkerEdgeColor',ecol, ...
                'LineStyle','none');
        else % plot modes 0-2
            plot(Amed,Bmed, ...
                'Marker',mrkr, ...
                'MarkerSize',msiz, ...
                'MarkerFaceColor',mcol, ...
                'MarkerEdgeColor',mcol, ...
                'LineStyle','none');
        end

        % Update extrema

        if initializing_plot
            xmn = Amed; xmx = Amed;
            ymn = min(blo(ibin),Bmed);
            ymx = max(bhi(ibin),Bmed);
            % scatter(NaN,NaN,msiz,mcol,'fill');
            hold on;
            initializing_plot = false;
        else
            xmn = min(xmn,Amed); xmx = max(xmx,Amed);
            ymn = min(ymn,min(blo(ibin),Bmed));
            ymx = max(ymx,max(bhi(ibin),Bmed));
        end
        
    end

end

hold off;

if isempty(xrng)
    xrng = plot_range([xmn xmx],0.05);
end

if isempty(yrng)
    ymx = max(abs([ymn ymx]));
    ymx = min(90,1.05*ymx);
    if (ymx > 75); ymx = 90; end
    ymn = -ymx;
    yrng = [ymn ymx];
end

xlim(xrng);
ylim(yrng);
   
box on;

xlabel('Altitude (km)','FontSize',xfsz,'FontWeight',xfwt);
% ylabel('Latitude (deg)','FontSize',yfsz,'FontWeight',yfwt);

if (yrng(2) == 90)
    set(gca,'YTick',(-90:30:90));
end
set(gca,'YTickLabel',[]);

set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

if altlat_time_cbar
    caxis([0 1]);
    cb = colorbar('Location','EastOutSide', ...
        'Ticks',[0 0.5 1],'TickLabels', ...
        {datestr(zrng(1),datefmt), ...
        datestr(mean(zrng),datefmt), ...
        datestr(zrng(2),datefmt)});
    set(cb,'FontSize',cfsz,'FontWeight',cfwt,'LineWidth',alwd);
    if ~isempty(zlabl)
        ylabel(cb,zlabl,'FontSize',zfsz,'FontWeight',zfwt);
    end
end

% Plot latitude distribution histogram in lower-right panel

subplot(3,20,[(21:26) (41:46)]);

if plotmode == 1
    
    % Green, yellow, and red histogram
    
    bcol = params.gcol; 
    % ecol = bcol;
    ecol = [0 0 0];
    barh(bmd,bct,1,'FaceColor',bcol,'EdgeColor',ecol);
    if any(bcy > 0)
        hold on;
        bcol = params.ycol;
        % ecol = bcol;
        barh(bmd,bcy,1,'FaceColor',bcol,'EdgeColor',ecol);
        hold off;
    end
    if any(bcr > 0)
        hold on;
        bcol = params.rcol;
        % ecol = bcol;
        barh(bmd,bcr,1,'FaceColor',bcol,'EdgeColor',ecol);
        hold off;
    end
    
else
    
    % Gray histogram
    
    bcol = [0.8 0.8 0.8];
    ecol = [0 0 0];
    barh(bmd,bct,1,'FaceColor',bcol,'EdgeColor',ecol);
 
end

% Label axes

xlabel('Frequency','FontSize',xfsz,'FontWeight',xfwt);
ylabel('Latitude (deg)','FontSize',yfsz,'FontWeight',yfwt);
set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

% Make y-ticks 30deg for pole-to-pole latitude span

ylim(yrng);
if (yrng(2) == 90)
    set(gca,'YTick',(-90:30:90));
end

% Plot Pc legend

render_legend = false;
title_region = (2:19);

if plotmode == 1
    
    % Green, yellow and red diamond legend
    
    render_legend = true;
    
    subplot(3,20,(16:20));
    title_region = (1:15);
    
    lmsiz = msiz+4;
    
    nleg = 0;

    mcol = params.rcol; ecol = mcol;
    nleg = nleg+1;
    hleg(nleg) = plot(NaN,NaN, ...
        'Marker',mrkr, ...
        'MarkerSize',lmsiz, ...
        'MarkerFaceColor',mcol, ...
        'MarkerEdgeColor',ecol, ...
        'LineStyle','none', ...
        'LineWidth',2);
    tleg{nleg} = ['Pc \geq ' ...
        expstr_to_pltstr(smart_exp_format(params.Pcred,3))];
    
    hold on;

    mcol = params.ycol; ecol = mcol;
    nleg = nleg+1;
    hleg(nleg) = plot(NaN,NaN, ...
        'Marker',mrkr, ...
        'MarkerSize',lmsiz, ...
        'MarkerFaceColor',mcol, ...
        'MarkerEdgeColor',ecol, ...
        'LineStyle','none', ...
        'LineWidth',2);
    tleg{nleg} = [expstr_to_pltstr(smart_exp_format(params.Pcyel,3)) ...
        ' \leq Pc < ' expstr_to_pltstr(smart_exp_format(params.Pcred,3))];
    
    mcol = params.gcol; ecol = mcol;
    nleg = nleg+1;
    hleg(nleg) = plot(NaN,NaN, ...
        'Marker',mrkr, ...
        'MarkerSize',lmsiz, ...
        'MarkerFaceColor',mcol, ...
        'MarkerEdgeColor',ecol, ...
        'LineStyle','none', ...
        'LineWidth',2);
    tleg{nleg} = ['Pc < ' ...
        expstr_to_pltstr(smart_exp_format(params.Pcyel,3))];
    
    hold off;

elseif plotmode == 3
    
    % Bubble legend
    
    render_legend = true;

    subplot(3,20,(17:20));
    title_region = (1:16);
    
    mcol = 'none';
    
    lPcbuf = floor(params.log10_Pcred):-1:ceil(params.log10_Pcyel);
    if (lPcbuf(1) == params.log10_Pcred)
        lPcbuf = lPcbuf(2:end);
    end
    if (lPcbuf(end) == params.log10_Pcyel)
        lPcbuf = lPcbuf(1:end-1);
    end
    
    ecol = [0 0 0];

    nleg = 0;
    
    msiz = msiz_Pcred;
    nleg = nleg+1;
    hleg(nleg) = plot(NaN,NaN, ...
        'Marker','o', ...
        'MarkerSize',msiz, ...
        'MarkerFaceColor',mcol, ...
        'MarkerEdgeColor',ecol, ...
        'LineStyle','none');
    tleg{nleg} = ['Pc \geq ' ...
        expstr_to_pltstr(smart_exp_format(params.Pcred,3))];

    hold on;
    
    for nbuf=1:numel(lPcbuf)
        msiz = msiz_Pcyel + (lPcbuf(nbuf)-params.log10_Pcyel)*msiz_delta;
        nleg = nleg+1;
        hleg(nleg) = plot(NaN,NaN, ...
            'Marker','o', ...
            'MarkerSize',msiz, ...
            'MarkerFaceColor',mcol, ...
            'MarkerEdgeColor',ecol, ...
            'LineStyle','none');
        tleg{nleg} = ['Pc = ' ...
            expstr_to_pltstr(smart_exp_format(10^lPcbuf(nbuf),3))];
    end
    
    msiz = msiz_Pcyel;
    nleg = nleg+1;
    hleg(nleg) = plot(NaN,NaN, ...
        'Marker','o', ...
        'MarkerSize',msiz, ...
        'MarkerFaceColor',mcol, ...
        'MarkerEdgeColor',ecol, ...
        'LineStyle','none');
    % tleg{nleg} = ['Pc \leq ' ...
    %     expstr_to_pltstr(smart_exp_format(params.Pcyel,3))];
    tleg{nleg} = ...
        [expstr_to_pltstr(smart_exp_format(params.Pc_cutoff_accum_risk,3)) ...;
        '-' expstr_to_pltstr(smart_exp_format(params.Pcyel,3))];
    hold off;
    
end

if render_legend
    
    axpos = get(gca,'OuterPosition');
    % lft = axpos(1);
    % bot = axpos(2);
    % wid = axpos(3);
    % hgt = axpos(4);
    axpos(3) = 1-axpos(1);
    axis off;
    
    legend(hleg,tleg,'Location','SouthEastOutSide', ...
        'Position',axpos, ...
       'FontSize',lfsz-1,'FontWeight',lfwt);
    legend('boxoff');
    
end

% Plot title in uppermost panel

subplot(3,20,title_region);

plot(NaN,NaN);
axis off;
Ntitl = numel(titl);
dy = 1/max(6,Ntitl);
y = 1-dy*((1:Ntitl)+0.5);
for n=1:Ntitl
    text(0.5,y(n),titl{n}, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
end

drawnow;

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