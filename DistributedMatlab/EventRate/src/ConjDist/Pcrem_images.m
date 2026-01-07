function Pcrem_images(L10PremImage,L10Premmd,L10Premlo,L10Premhi, ...
                      YvalImage,Yval0,PcumImage,RmvrImage,Ylabl,Titl)
% Pcrem_images - Render images of Pcum and Nmvr as a function of L10Prem 
% (x-axis) and a y-axis variable (e.g., Tmis, Fred, or Tcom)
%
% Syntax: Pcrem_images(L10PremImage,L10Premmd,L10Premlo,L10Premhi, ...
%                      YvalImage,Yval0,PcumImage,RmvrImage,Ylabl,Titl)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Render images of Pcum and Nmvr as a function of L10Prem 
% (x-axis) and a y-axis variable (e.g., Tmis, Fred, or Tcom)
%
% =========================================================================
%
% Input:
%
%   L10PremImage - Array of Log10 Pc values to plot                [Np x 1]
% 
%   L10Premmd    - median Pc value at nominal solution
% 
%   L10Premlo    - low percentile Pc value at nominal solution
% 
%   L10Premhi    - high percentile Pc value at nominal solution
% 
%   YvalImage    - Array of y-values to plot                       [Ny x 1]
% 
%   Yval0        - Nominal y-value
% 
%   PcumImage    - Array of Pc values corresponding to each       [Np x Ny]
%                  Prem-Fred combination
%
%   RmvrImage    - Array of RMM values corresponding to each      [Np x Ny]
%                  Prem-Fred combination. Can be set as empty  
%                  to skip RMM rate plotting
% 
%   Ylabl        - y-axis label
% 
%   Titl         - Plot title
% =========================================================================
%
% Output:
%
%   None
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

% Initialize plot

clf;

xfsz = 11;
yfsz = xfsz;

xfwt = 'bold';
yfwt = xfwt;

tfsz = xfsz;
tfwt = 'bold';
tfan = 'italic';

lfsz = xfsz;
lfwt = 'bold';

afsz = xfsz;
afwt = 'bold';
alwd = 1;

% Contour and solution-point colors

contcolr = [0 0 0];
solncolr = 'b';
solnmrkr = 'o';
solnmsiz = 6;

% Yellow-to-red colormap

Nmap = 256;
cmap = zeros(Nmap,3);
cmap(:,1) = 1;
cmap(:,3) = 0;
cmap(:,2) = cmap(:,3)+(1-cmap(:,3)).*((0:Nmap-1)'/(Nmap-1));
cmap = flipud(cmap);
colormap(cmap);

% Determine if plotting RMM rate or not

plotRMMrates = ~isempty(RmvrImage);

% Image: x=L10Prem, y=Yval, z=Pcum

if plotRMMrates
    subplot(3,1,2,'replace');
else
    subplot(3,1,[2 3],'replace');
end

xrng = [L10PremImage(1) L10PremImage(end)];
yrng = [YvalImage(1) YvalImage(end)];

img = log10(PcumImage)';
imagesc(L10PremImage,YvalImage,img);
axis xy;

hold on;

% Overplot contours of log10(Pcum)

mincon = min(-1,floor(L10PremImage(1)));
maxcon = -1;

% Plot 1e-N solid-line contours with labels

seqcon = (mincon:1:maxcon);
contour(L10PremImage,YvalImage,img,seqcon, ...
    'LineStyle','-','LineColor',contcolr,'ShowText','On');

% Plot 3e-N dotted-line contours without labels

seqcon = seqcon + log10(3);
contour(L10PremImage,YvalImage,img,seqcon, ...
    'LineStyle',':','LineColor',contcolr,'ShowText','Off');

% Mark nominal solution

if ~isempty(L10Premmd)
    plot(L10Premmd,Yval0,solnmrkr,'MarkerSize',solnmsiz, ...
        'MarkerFaceColor',solncolr,'MarkerEdgeColor',solncolr);
    plot([L10Premlo L10Premhi],[Yval0 Yval0],'-','Color',solncolr);
end

hold off;

xlim(xrng);
ylim(yrng);

xlabl = 'Log_{10}[Remediation Threshold Pc]';
ylabl = Ylabl;

xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
set(gca,'LineWidth',alwd,'FontSize',afsz,'FontWeight',afwt);

cb = colorbar;
cb.Label.String = 'Log_{10}[Cumulative Pc]';
cb.FontSize = yfsz;

crng = [min(img(:)),max(img(:))];
if crng(1) == crng(2)
    crng(1) = crng(1)-1e-3;
    crng(2) = crng(2)+1e-3;
end
caxis(crng);

% Image: x=L10Prem, y=Yval, z=Nmvr

if plotRMMrates

    subplot(3,1,3,'replace');

    img = log10(RmvrImage)';
    imAlpha = ones(size(img));
    ndx = isinf(img);
    imAlpha(ndx) = 0;        
    imagesc(L10PremImage,YvalImage,img,'AlphaData',imAlpha);
    gscl = (1+max(sum(cmap,2))/3)/2;
    set(gca,'color',gscl*[1 1 1]);
    axis xy;

    hold on;

    % Overplot contours of log10(Nmvr)

    ndx = ~ndx;
    if any(ndx)
        crng = [min(img(ndx)) max(img(ndx))];
        if crng(1) > 0
            crng(1) = 0;
        else
            crng(1) = max(crng(2)-2,crng(1));
        end
    else
        crng = [-1 1];
    end

    mincon = floor(crng(1));
    maxcon = ceil(crng(2));

    seqcon = (mincon:1:maxcon);
    contour(L10PremImage,YvalImage,img,seqcon, ...
        'LineStyle','-','LineColor',contcolr,'ShowText','On');

    seqcon = seqcon + log10(3);
    contour(L10PremImage,YvalImage,img,seqcon, ...
        'LineStyle',':','LineColor',contcolr,'ShowText','Off');


    if ~isempty(L10Premmd)
        plot(L10Premmd,Yval0,solnmrkr,'MarkerSize',solnmsiz, ...
            'MarkerFaceColor',solncolr,'MarkerEdgeColor',solncolr);
        if (L10Premlo ~= L10Premmd)
            plot([L10Premlo L10Premmd],[Yval0 Yval0],'-','Color',solncolr);
        end
        if (L10Premhi ~= L10Premmd)
            plot([L10Premhi L10Premmd],[Yval0 Yval0],'-','Color',solncolr);
        end
    end

    hold off;

    xlim(xrng);
    ylim(yrng);

    xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
    ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
    set(gca,'LineWidth',alwd,'FontSize',afsz,'FontWeight',afwt);

    cb = colorbar;
    cb.Label.String = 'Log_{10}[RMM Rate (yr^{-1})]';
    cb.FontSize = yfsz;

    % caxis([0,max(img(:))]);

    caxis(crng);
    
end

% Title

subplot(3,1,1,'replace');

axpos = get(gca,'OuterPosition');
axpos(2) = axpos(2)-0.025;
axpos(4) = 0.001;
set(gca,'OuterPosition',axpos)
axis off;

title(Titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);

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