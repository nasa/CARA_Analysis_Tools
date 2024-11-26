function PlotMetricVsSDACurves(PlotMatrix,params)

% PlotMetricVsSDACurves - Plot versus Solar Depression Angle (SDA)
% Syntax: PlotMetricVsSDACurves(PlotMatrix,params);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: 
% Plot of a curve versus SDA for visualizing metrics such as number of 
% satellites sunlit (Ni) or brighter than threshold (Nb). 
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Plot ID strings
EvalID = params.EvalID; RunID = params.RunID;

% Set up for plotting
figure(3); clf;
xfsz = 11;
yfsz = xfsz;
xfwt = 'bold';
yfwt = xfwt;
tfsz = xfsz;
tfwt = 'bold';
tfan = 'italic';
afsz = xfsz;
afwt = 'bold';
alwd = 0.5;
lfsz = xfsz;
lfwt = 'bold';
% msiz = 7;

xlabl = 'Solar Depression Angle (deg)';

% Comparison constellations
if isempty(params.Plotting.Comparison)
    Ncomp = 0;
else
    % Number of compared data sets
    Ncomp = numel(params.Plotting.Comparison);
    % Get the comparison data sets
    tabcomp = cell(Ncomp,1);
    tabfile0 = [params.RunID '_SDAPlot.xlsx']; 
    for ncomp=1:Ncomp
        tabpath = fullfile('output',params.Plotting.Comparison{ncomp});
        tabfile = strrep(tabfile0,params.EvalID, ...
                         params.Plotting.Comparison{ncomp});
        tabfull = fullfile(tabpath,tabfile);
        if ~exist(tabfull,'file')
            error(['Comparison plot table not found: ' tabfile]);
        else
            tabcomp{ncomp} = readtable(fullfile(tabpath,tabfile));
        end
    end
end

% Number constellation satellites and orbital shells
Nctot = sum(params.New.Nc);
Nshell = numel(params.New.Nc);

% Title or legend
clear titl slegnd;
if Ncomp == 0
    nt=1; titl{nt} = strrep(EvalID,'_','\_');
    nt=nt+1; titl{nt} = ['Number of Constellation Satellites = ' ...
        num2str(Nctot)];
    if Nshell == 1
        nt=nt+1; titl{nt} = ['Constellation Altitude = ' ...
            num2str(params.New.Altitude_km) ' km'];
    else
        nt=nt+1; titl{nt} = ['Number of Orbital Shells = ' ...
            num2str(Nshell) ', Altitudes = ' ...
            num2str(min(params.New.Altitude_km)) ' to ' ...
            num2str(max(params.New.Altitude_km)) ' km'];
    end
    nt=nt+1; titl{nt} = ' ';
else
    slegnd = cell(Ncomp+1,1);
    nl=1; slegnd{nl} = strrep(EvalID,'_','\_');
    for ncomp=1:Ncomp
        nl=nl+1;
        slegnd{nl} = strrep(params.Plotting.Comparison{ncomp},'_','\_');
    end
end

% Make the plots
Nplt = 2;
for nplt=1:Nplt

    % Y axis points
    if nplt == 1 && params.Plotting.NsunlitVsSDAPlot
        make_plot = true;
        ydat = PlotMatrix.NsunlitPlot;
        splt = 'NiVsSDA';
        ylabl = cell(1,2);
        if params.Evaluation.UniformDist == 2
            ylabl{1} = 'Number Sunlit Above Observer';
            ylabl{2} = '(uniform thin shell approx.)';
        else
            ylabl{1} = 'Number of Sunlit Constellation';
            ylabl{2} = 'Satellites (Global Peak)';
        end
    elseif nplt == 2 && params.Plotting.NbrightVsSDAPlot
        make_plot = true;
        ydat = PlotMatrix.NbrightPlot;
        splt = 'NbVsSDA';
        if params.Evaluation.UniformDist == 2
            ylabl = cell(1,2);
            ylabl{1} = 'Number Brighter than Recommended';
            ylabl{2} = '(uniform thin shell approx.)';
        else
            % ylabl = cell(1,2);
            % ylabl{1} = 'Number of Brighter than Recommended';
            % ylabl{2} = 'Constellation Satellites (Global Peak)';
            ylabl = 'N_b (Global Peak Value)';
        end
        % Augment title for Nb indicator
        if Ncomp == 0
            nt = numel(titl);
            if params.Evaluation.MaxZenith ~= 90
                titl{nt} = ['N_b = number brighter than recommended within ' ...
                    num2str(params.Evaluation.MaxZenith) '{\circ} of zenith'];
            else
                titl{nt} = 'N_b = Number brighter than recommended above observer';
            end
            % Get global peak Nb value as SDA = 18 deg
            NbGP = interp1(PlotMatrix.SolarDepression,ydat(:,1),18);
            nt=nt+1; titl{nt} = ['Global peak for SDA = 18\circ: N_b = ' ...
                smart_exp_format(NbGP,3)];
            nt=nt+1; titl{nt} = ' ';
        end
    else
        make_plot = false;
    end

    % Only make plot if required
    if make_plot

        % X axis points (renewed for every plot on purpose)
        xplt = PlotMatrix.SolarDepression;
        if params.Plotting.CustomPlots
            xrng = [0 54];
        else
            xrng = [min(xplt) max(xplt)];
        end

        % Get the extrema of the comparison constellation curves
        if Ncomp > 0
            ycmin = Inf; ycmax = -Inf;
            for ncomp = 1:Ncomp
                if nplt == 1
                    xcomp = tabcomp{ncomp}.SolarDepression;
                    ycomp = tabcomp{ncomp}.NsunlitPlot;
                elseif nplt == 2
                    xcomp = tabcomp{ncomp}.SolarDepression;
                    ycomp = tabcomp{ncomp}.NbrightPlot;
                end
                ycmin = min(ycmin,min(ycomp));
                ycmax = max(ycmax,max(ycomp));
            end
        end

        % Y axis scaling
        if params.Plotting.PlotMatrixLatitudeBands
            NLatBand = size(ydat,2)-1;
            % Check for four columns, global, low, med and high latitudes
            if NLatBand ~= 0 && NLatBand ~= 3
                error('Invalid number of latitude bands in PlotMatrix');
            end
            % Use first column, which is global latitude band
            yplt = ydat(:,1);
        else
            % Sum over columns, each containing one shell
            yplt = sum(ydat,2);
        end
        ymax = max(yplt);
        ylog = true;
        if Ncomp > 0
            ymax = max(ymax,ycmax);
        end
        ymn0 = 0.1;
        if nplt == 1 
            ymin = max(ymn0,min(yplt(yplt > 0))/2); 
        else
            ymin = ymn0;
        end
        yrng = [ymin max(params.Evaluation.MaxBright)];
        if ymax > yrng(2)
            yrng(2) = ymax*2;
        elseif ymax < yrng(1)
            yrng(1) = ymax/2;
        end

        % Initialized plot
        plot(NaN,NaN); hold on;

        % Plot the new constellation's curve
        Ncurve = 1;
        hlegnd = CurveVsSDAPlot(xplt,yplt,nplt,Nplt,Ncurve,Ncomp,yrng,params);
        
        % Plot the comparison constellation curves
        if Ncomp == 0
            if params.Plotting.PlotMatrixLatitudeBands && NLatBand == 3
                % Allocate legend handles
                hlegnd = repmat(hlegnd,[NLatBand 1]);
                % Plot latitude band curves
                for nLatBand=1:NLatBand
                    Ncurve = 1+nLatBand;
                    hlegnd(Ncurve) = ...
                        CurveVsSDAPlot(xplt,ydat(:,nLatBand+1), ...
                        nplt,Nplt,Ncurve,Ncomp,yrng,params);
                end
                slegnd = cell(size(hlegnd));
                slegnd{1} = 'All Latitudes';
                slegnd{2} = 'Low Latitudes';
                slegnd{3} = 'Medium Latitudes';
                slegnd{4} = 'High Latitudes';
            else
                slegnd = [];
                if Nshell > 1 && size(ydat,2) == Nshell
                    for nshell=1:Nshell
                        Ncurve = 1+nshell;
                        hlegnd(Ncurve) = ...
                            CurveVsSDAPlot(xplt,ydat(:,nshell), ...
                            nplt,Nplt,Ncurve,Ncomp,yrng,params);
                    end
                end
            end
        else
            % Allocate legend handles
            hlegnd = repmat(hlegnd,size(slegnd));
            % Plot comparison curves
            for ncomp = 1:Ncomp
                if nplt == 1 && params.Plotting.NsunlitVsSDAPlot
                    xcomp = tabcomp{ncomp}.SolarDepression;
                    ycomp = tabcomp{ncomp}.NsunlitPlot;
                elseif nplt == 2 && params.Plotting.NbrightVsSDAPlot
                    xcomp = tabcomp{ncomp}.SolarDepression;
                    ycomp = tabcomp{ncomp}.NbrightPlot;
                end
                Ncurve = 1+ncomp;
                hlegnd(Ncurve) = CurveVsSDAPlot(xcomp,ycomp,nplt,Nplt, ...
                    Ncurve,Ncomp,yrng,params);
            end
            % Overplot main curve
            Ncurve = -1;
            CurveVsSDAPlot(xplt,yplt,nplt,Nplt,Ncurve,Ncomp,yrng,params);
        end

        % Redraw axes over shaded regions
        if nplt == Nplt        
            plot(xrng,[yrng(1) yrng(1)],'-k','LineWidth',alwd);
            plot(xrng,[yrng(2) yrng(2)],'-k','LineWidth',alwd);
            plot([xrng(1) xrng(1)],yrng,'-k','LineWidth',alwd);
            plot([xrng(2) xrng(2)],yrng,'-k','LineWidth',alwd);
        end
        hold off;

        % Axes
        xlim(xrng);
        xticks(xrng(1):6:xrng(2));
        if ylog
            set(gca, 'YScale', 'log');
        else
            set(gca, 'YDir', 'reverse');
        end
        ylim(yrng);
        
        % grid on;
        % set(gca,'MinorGridLineStyle','-')
        % set(gca,'GridLineStyle','-')
        
        Ax = gca;
        Ax.XGrid = 'on';        
        Ax.YGrid = 'on';
        Ax.Layer = 'top';
        Ax.GridColor = 0.5*[1 1 1];
        Ax.MinorGridColor = 0.6*[1 1 1];
        Ax.GridAlpha = 0.5;
        Ax.GridLineStyle = '-';
        Ax.MinorGridLineStyle = '-';

        % Labels
        xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
        ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
        set(gca,'LineWidth',alwd,'FontSize',afsz,'FontWeight',afwt);        

        % Title or legend
        if Ncomp == 0
            title(titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
            if ~isempty(slegnd)
                legend(hlegnd,slegnd,'Location','NorthEast', ...
                    'FontSize',lfsz,'FontWeight',lfwt);
            end
        else
            % Special legend modifications for AMOS paper plots
            if params.Plotting.CustomPlots
                for n=1:numel(slegnd)
                    [prts,~] = string_parts(slegnd{n},'\');
                    if strcmpi(prts{1},'Starlink1stGen')
                        slegnd{n} = 'Starlink 1st Generation';
                    elseif strcmpi(prts{1},'Starlink2ndGen')
                        slegnd{n} = 'Starlink 2nd Generation';
                    elseif strcmpi(prts{1},'OneWebPhase1')
                        slegnd{n} = 'OneWeb Phase 1';
                    elseif strcmpi(prts{1},'OneWebPhase2')
                        slegnd{n} = 'OneWeb Phase 2';
                    elseif strcmpi(prts{1},'Iridium2ndGen')
                        slegnd{n} = 'Iridium 2nd Generation';
                    end
                end                
            end
            legend(hlegnd,slegnd,'Location','NorthOutSide',...
                'FontSize',lfsz,'FontWeight',lfwt);
        end

        % Save plot
        drawnow;
        if Ncomp == 0
            pltfile = [RunID '_' splt '.png'];
        else
            pltfile = [RunID '_Ncomp' num2str(Ncomp) '_' splt '.png'];
        end
        saveas(gcf,fullfile(params.Output.outputpath,pltfile));

    end

end

return
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