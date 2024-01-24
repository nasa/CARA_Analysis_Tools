function GenConjPlanePlot(PcFoster,ConjData)
%
% GenConjPlanePlot - This function generates the Conjunction Plane Plot.
%
% Syntax:            GenConjPlanePlot(PcFoster,ConjData)
%
% Inputs:
%   PcFoster       - 2D probability of collision (Pc) value generated using Foster method
%                   (PcFosterMethod.m)
%   ConjData       - A vector of values produced by PcFosterMethod.m
%                   [Contains: semimajor & semiminor axes (conjunction ellipse),
%                    clock angle, HBR, miss distance, x1sig, 
%                    and RIC sigma values, condition numbers, and approach angle]
%
% Output:
%   []             - Matlab figure (FIG file)
%
% Examples/Validation Cases: 
%
% Other m-files required: None (However, the inputs are obtained by using other m-files)
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Author: Dragan Plakalovic
% E-Mail: Dragan.Plakalovic@omitron-cos.com
% Jul 2015; Last revision: 13-Jul-2015
%
% ----------------- BEGIN CODE -----------------

    % CREATE PLOT
    
    if (isvector(PcFoster))
        PcCircle = PcFoster(1);
    else
        PcCircle = PcFoster;
    end
    
    if (isnan(PcCircle))
        fprintf('PcCircle has value of NaN. Conjunction plane plot cannot be generated.\n')
        return
    end    

    % Conjuction Data
    a          = ConjData(1);
    b          = ConjData(2);
    ClockAngle = ConjData(3);
    HBR        = ConjData(4);
    MissDist   = ConjData(5);
    MinHBSigmaDistance = (MissDist - HBR) / ConjData(6);

    % Plot settings
    AxisLblSz  = 19;
    TitleSz    = 21;
    LineWidth  = 2;
    FontWeight = 'Normal';
    FontName   = 'Calibri';
    
    % Generate basis ellipse (parameterization)
    u  = (0:0.1:2*pi+0.1)';
    x0 = a * cos(u);
    y0 = b * sin(u);
    
    % Rotate ellipse by clock angle
    x = x0.*cosd(ClockAngle) - y0.*sind(ClockAngle); 
    y = x0.*sind(ClockAngle) + y0.*cosd(ClockAngle);
    
    % Construct HBR
    xHBR = HBR * cos(u) + MissDist;
    yHBR = HBR * sin(u);
    
    % Creating figure
    figure('NumberTitle','off','Name','2D Conjunction Plane','Position',[500,250,840,700]);
    
    hold on
    
    % Plot 1, sqrt(2), and 3-sigma ellipses
    plot(x,y,'Color',[1 0 0],'LineWidth',LineWidth);
    plot(sqrt(2)*x, sqrt(2)*y,'Color',[1 0.8 0.2],'LineWidth',LineWidth);
    plot(3*x,3*y,'Color',[0 1 0.45],'LineWidth',LineWidth);
    
    % Plot HBR
    plot(xHBR,yHBR,'b--','LineWidth',LineWidth);
    
    hold off
    
    % Set axes limits to ensure the plot traces don't appear on the edge of the plot area
    UpperXLim = max([max(3*x) max(xHBR)]) + 0.05*(max([max(3*x) max(xHBR)] - min([min(3*x) min(xHBR)])));
    LowerXLim = min([min(3*x) min(xHBR)]) - 0.05*(max([max(3*x) max(xHBR)] - min([min(3*x) min(xHBR)])));
    UpperYLim = max([max(3*y) max(yHBR)]) + 0.05*(max([max(3*y) max(yHBR)] - min([min(3*y) min(yHBR)])));
    LowerYLim = min([min(3*y) min(yHBR)]) - 0.05*(max([max(3*y) max(yHBR)] - min([min(3*y) min(yHBR)])));
    set(gca,'XLim',[LowerXLim UpperXLim],'YLim',[LowerYLim UpperYLim]);
    
    % Labels
    LegendHBR = strcat('Hard Body Region (',num2str(HBR*1000),' m)');
    title('Conjunction Plane','FontSize',TitleSz,'FontWeight',FontWeight,'FontName',FontName);
    xlabel('Miss Distance Direction [km] ','FontSize',AxisLblSz,'FontWeight',FontWeight,'FontName',FontName);
    ylabel('Relative Out of Plane [km]'   ,'FontSize',AxisLblSz,'FontWeight',FontWeight,'FontName',FontName);
    legend('One Sigma Covariance', 'SQRT Two Sigma Covariance', 'Three Sigma Covariance',LegendHBR,'Location','NorthEast');
    text('Units','normalized','Position',[0.01,-0.11],'FontSize',12,'FontName',FontName,'FontWeight','Bold', ...
         'String',sprintf('Pc Value : %05.2e',PcCircle));
    text('Units','normalized','Position',[0.415,-0.11],'FontSize',12,'FontName',FontName,'FontWeight','Bold', ...
         'String',sprintf('Sigma : %04.2f',MinHBSigmaDistance));
    text('Units','normalized','Position',[0.750,-0.11],'FontSize',12,'FontName',FontName,'FontWeight','Bold', ...
         'String',sprintf('Combined HBR [m] : %g',HBR*1000));
    axis equal
    grid on
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  | 07-13-2015 |  Initial Development
%  