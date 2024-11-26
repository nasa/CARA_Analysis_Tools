function PlotLightPollutionTradeSpace(PlotInfo,MMT,MMTpars,params)
% PlotLightPollutionTradeSpace
% Syntax: PlotLightPollutionTradeSpace(PlotInfo,MMT,MMTpars,params);
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
% Plot constellation light polution trade space image
%   x = horizontal-axis = Nc
%   y = vertical-axis = altidude
%   z = color-axis = light pollution level (composite score)
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Extract plot info
EvalID = PlotInfo.EvalID; RunID = PlotInfo.RunID;
use_analog_satellites = PlotInfo.use_analog_satellites;

% X axis of image is Nc * FluxScalingFactor
x1 = log10(params.TradeSpace.NcF_axis(1));
x2 = log10(params.TradeSpace.NcF_axis(2));
Nx = params.TradeSpace.NcF_axis(3);
x = linspace(x1,x2,Nx);

% Y axis of image is altitude in km
y1 = log10(params.TradeSpace.hkm_axis(1));
y2 = log10(params.TradeSpace.hkm_axis(2));
Ny = params.TradeSpace.hkm_axis(3);
y = linspace(y1,y2,Ny);

% Image array
z = NaN(Ny,Nx);

% Arrays for scoring
params.Evaluation.MaxBright(params.Evaluation.MaxBright < 0) = NaN;
ndx = ~isnan(params.Evaluation.MaxBright);
MaxBright = params.Evaluation.MaxBright(ndx);
Nev = numel(MaxBright);
MaxBright = reshape(MaxBright,[Nev 1]);
SDAPoints = params.Evaluation.SDAPoints(ndx);
Nabove  = NaN(Nev,1);
Nsunlit = NaN(Nev,1);
Nbright = NaN(Nev,1);

% Set up for analog satellites
if use_analog_satellites    
    MMTpars.previous_output = MMT;
    MMTpars.Fscale = 1;
    MMTpars.Ncomp = 0;
    MMTpars.Plotting = false;
    MMTpars.UseMMTSaveFiles = false;
    MMTpars.StatAnalysis = false;
    MMTpars.verbose = false;
end
MaxZenith = params.Evaluation.MaxZenith;

% Calculate the tradespace image

% Loop over altitudes
for ny=1:Ny

    % Altitude
    altkm = 10^y(ny);

    % Calculate Fb
    if use_analog_satellites    

        MM2 = AnalyzeMMTData(params.Analog.datapath,           ...
                             altkm,                            ...
                             params.Analog.UTbegin,            ...
                             params.Analog.UTend,              ...
                             params.Evaluation.ExtinctionCoef, ...
                             MMTpars);

        Fb = MM2.FracBrighterThanSatCon1Limit;
        FbGrid2 = MM2.FracBrightGrid;

    else

        % Magnitude statistics
        dM = 5*log10(params.New.Altitude_km/params.New.MzenAltkm);        
        Mmd = dM + params.New.Mzen50;
        Mlo = dM + min(params.New.Mzen05,params.New.Mzen95);
        Mhi = dM + max(params.New.Mzen05,params.New.Mzen95);

        % SatCon-1 limit
        Msc = SatCon1Limit(altkm);

        % Estimate FracBrighterThanSatCon1Limit using asymmetric Gaussian with
        % sigmas that match the Mz05 and Mz95 points
        Madj = -5*log10(params.New.Altitude_km/altkm);

        Fb = estimateFb(Msc,Mmd+Madj,0.05,Mlo+Madj,0.95,Mhi+Madj);

        % Calculate fraction of off-zenith mags bright than cutoff
        Nzag = numel(params.Grid.ZenAngGrid);
        Rekm = 6378.137;
        Rnormkm = params.New.Altitude_km;
        akm = Rekm+Rnormkm; akm2 = akm^2;
        FbGrid2 = NaN(size(params.Grid.ZenAngGrid));
        for nzag=1:Nzag
            czag = cos(params.Grid.ZenAngGrid(nzag));
            szag = sin(params.Grid.ZenAngGrid(nzag));
            rho = sqrt(akm2-(Rekm*szag)^2) - Rekm*czag;
            dM = Madj+5*log10(rho/Rnormkm);
            FbGrid2(nzag) = estimateFb(Msc,Mmd+dM,0.05,Mlo+dM,0.95,Mhi+dM);
        end

    end

    % Interpolation grid for Fb
    if params.Grid.FracBrightInterp
        % Interpolate brightness using variable
        ZaInterp = params.Grid.ZenAngGrid;
        FbInterp = FbGrid2;
    else
        % Interpolate brightness using constant Fb = Fbnew
        ZaInterp = [0 pi/2];
        FbInterp = [Fb Fb];
    end

    % Evaluate constellation for Nc = 1
    for n=1:Nev
        % Number of new constellation satellites above a low-latitude
        % observer calculated with the uniform shell approximation
        [Nabove(n),Nsunlit(n),Nbright(n)] = NumExpectedUniformShell( ...
            1,altkm,MaxZenith,ZaInterp,FbInterp,-SDAPoints(n));
    end

    % Nbright = Fb*Nsunlit;
    LightPollutionLevel = max(Nbright./MaxBright);

    for nx=1:Nx

        % Nc * FluxScalingFactor
        NcF = 10.^x(nx);

        % Calculate the light pollution level
        z(ny,nx) = NcF*LightPollutionLevel;

    end

end

% Make plot

figure(4); clf;

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

msiz = 7;

zplt = log10(z);
zrng = [-1 1];
zplt(zplt < zrng(1)) = zrng(1);
zplt(zplt > zrng(2)) = zrng(2);

imagesc(x,y,zplt);    
axis xy;

% Colormap
cmap = flipud(hot(256)); colormap(cmap);

% xlabl = 'Log_{10}[Number of Constellation Satellites]';
xlabl = 'Number of Constellation Satellites';
xticks(log10([10 30 100 300 1e3 3e3 1e4]));
xticklabels({'10','30','100','300','1e3','3e3','1e4'});
hax = gca;
hax.XAxis.MinorTickValues = log10([1e1 2e1 3e1 4e1 5e1 6e1 7e1 8e1 9e1 ...
                                   1e2 2e2 3e2 4e2 5e2 6e2 7e2 8e2 9e2 ...
                                   1e3 2e3 3e3 4e3 5e3 6e3 7e3 8e3 9e3 ...
                                   1e4]);
hax.XAxis.MinorTick = 'On';

% ylabl = 'Log_{10}[Altitude (km)]';
ylabl = 'Altitude (km)';
yticks(log10([300 1e3 3e3 1e4]));
yticklabels({'300','1e3','3e3','1e4'});
hax.YAxis.MinorTickValues = log10([1e1 2e1 3e1 4e1 5e1 6e1 7e1 8e1 9e1 ...
                                   1e2 2e2 3e2 4e2 5e2 6e2 7e2 8e2 9e2 ...
                                   1e3 2e3 3e3 4e3 5e3 6e3 7e3 8e3 9e3 ...
                                   1e4]);
hax.YAxis.MinorTick = 'On';

% Overplot contours
hold on;
gray = 0.7*[1 1 1];
seqcon = [0 0];
contour(x,y,zplt,seqcon, ...
    'LineStyle','-','LineColor',gray,'LineWidth',1,'ShowText','Off');    
seqcon = log10([params.Evaluation.LowToHigh(2) ...
                params.Evaluation.LowToHigh(4)]);
contour(x,y,zplt,seqcon, ...
    'LineStyle','--','LineColor',gray,'LineWidth',1,'ShowText','Off');    
plot(log10(params.New.Nc),log10(params.New.Altitude_km),'o', ...
     'MarkerFaceColor',gray,'MarkerEdgeColor','k','MarkerSize',msiz);
plot(log10(params.New.Nc),log10(params.New.Altitude_km),'+', ...
     'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',msiz);
hold off;

grid on;

xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
set(gca,'LineWidth',alwd,'FontSize',afsz,'FontWeight',afwt);

clear titl;
nt=1; titl{nt} = strrep(EvalID,'_','\_');
nt=nt+1; titl{nt} = ['Light pollution level = ' ...
   smart_exp_format(PlotInfo.LightPollutionLevel,3)];
nt=nt+1; titl{nt} = ['Light pollution risk = ' ...
   PlotInfo.LightPollutionRisk];
k = strfind(PlotInfo.Recommendation,' to mitigate');
recstr = PlotInfo.Recommendation(1:k);
% nt=nt+1; titl{nt} = ['RECOMMENDATION: ' recstr];
nt=nt+1; titl{nt} = ' ';
title(titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);

cb = colorbar;
cb.FontSize = yfsz;

crng = zrng;
caxis(crng);

% cb.Label.String = 'Log_{10}[Light Pollution Level]';
cb.Label.String = 'Light Pollution Level';
set(cb,'YTick',log10([0.1 0.3 1.0 3.0 10]));
set(cb,'YTickLabel',{'0.1' '0.3' '1.0' '3.0' '10'});

pltfile = [RunID '_TradeSpace.png'];
saveas(gcf,fullfile(params.Output.outputpath,pltfile));

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