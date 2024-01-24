function [Status] = CollisionConsequenceFigures(PrimaryMass,VRel,OutputFolder,OutputID,Catastrophic,NumOfPieces,massVec)
%
% CollisionConsequenceFigures - Generates Figures for representing the
% collision consequence numbers including the number of pieces
% distribution and secondary object mass distribution
%
% Syntax:   [Status] = CollisionConsequenceFigures(PrimaryMass,VRel,OutputFolder,OutputID,Catastrophic,NumOfPieces,massVec)
%
% Inputs:
%   PrimaryMass     - 1X1 Mass of Primary Object (kg)
%   VRel            - 1x1, 3X1, or 1X3 vector of the relative velocity
%                     between the primary and secondary objects (m/s)
%   OutputFolder    - String location of output file directory folder where
%                     figures will be saved (optional)
%                       - If no output location specified, figures will not
%                           be saved but will remain on screen
%                       - If output location specified, figures will be
%                           saved and then automatically closed
%   OutputID        - String for output figure ID (optional) used to
%                     identify specific event
%   Catastrophic    - [NumOfSamplesX1] logical array indicating whether the
%                     sampled collision is catastrophic
%   NumOfPieces     - [NumOfSamplesX1] array of the number of pieces
%                     expected to be generated from a collision for each sample
%   massVec         - [NumOfSamplesX1] array of the secondary object mass
%                     estimates for each individual sample
%
% Outputs:
%   Status          - Success Status of operation
%                       1 - Figures generated Successfully
%                       0 - Error Encountered
%
% Example/Validation Cases:
%
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: 	None
%                           
% Subfunctions: None
% MAT-files required:       None
%
% See also: none
%
% May 2018; Last revision: 22-May-2018
%
% ----------------- BEGIN CODE -----------------
    
try    
    %% Set Constants
    CatastrophicThreshold               = 40000; % Threshold for Catastrophic Collision (Joules/kg)
    CatastrophicPlotGenerationThreshold = 0.05; % Threshold for whether or not to generate a distinct plot for Catastrophic or Non-Catastrophic Debris Clouds
    CatastrophicThresholdMass           = CatastrophicThreshold*2*PrimaryMass/norm(VRel)^2;
    if sum(Catastrophic)/length(Catastrophic)*100 > 0.005
        CollisionProbabilityText = num2str(sum(Catastrophic)/length(Catastrophic)*100,'%.2f');
    else
        CollisionProbabilityText = '<0.005';
    end
    if ~isempty(OutputFolder)
        CloseFigures            = 1; % Toggle Variable for whether to close figures after generation and saving
    else
        CloseFigures            = 0;
    end
    
    %% Close Existing Collision Consequence figure windows if open
    try
        close('Secondary Mass Distribution')
    catch
    end
    try
        close('Non-Catastrophic Debris Piece Distribution')
    catch
    end
    try
        close('Catastrophic Debris Piece Distribution')
    catch
    end
    
    %% Generate Number of Pieces Distribution Plots for Non-Catastrophic Collisions
    CDFcounts = [0:100/length(NumOfPieces):100];
    % Only generate figure if Non-Catastrophic Collision Probability above 5%
    if (1-sum(Catastrophic)/length(Catastrophic))*100 > 5
        idx       = find(~Catastrophic);
        ExpectedNumberofPieces    = round(10^mean(log10(NumOfPieces(idx))));
        h= figure('Position',[75 75 550 700],...
                   'Name','Non-Catastrophic Debris Piece Distribution');
        % generate Histogram of distribution
        subplot(2,1,1);
        % define bin edges expecting to represet a chi square distribution
        edges = [0:median(NumOfPieces(idx))/4:median(NumOfPieces(idx))*6];
        histogram(NumOfPieces(idx),edges)
        title({'Expected Number of Debris Pieces Generated'; 'If a Non-Catastrophic Collision Occurs'})
        xlabel('Number of Pieces')
        ylabel({'Incidence Count ';['(' num2str(length(NumOfPieces),'%.0f') ' Total Trials)']})
        box on
        grid on
        % Add Text Indicating probability of Catastrophic Collision
        yl = ylim;
        xl = xlim;
        text(xl(2)*0.98,yl(2)*0.97,...
            {[CollisionProbabilityText '% Chance of Catastrophic Collision'];
            [num2str(ExpectedNumberofPieces,'%.0f') ' Expected Debris Pieces '];
             '(Non-Catastrophic Collision)'},...
            'HorizontalAlignment','right',...
            'VerticalAlignment','top',...
            'EdgeColor','k',...
            'BackgroundColor',[1 1 1])
        % Generate CDF of distribution
        subplot(2,1,2);
        hold on
            % Plot CDF
            plot([0 reshape(sort(NumOfPieces(idx),'ascend'),1,length(NumOfPieces(idx))) xl(2)],...
                [CDFcounts(1:length(idx)+1) CDFcounts(length(idx)+1)],...
                'LineWidth',2,...
                'Color','b')
        hold off
        title({'Cumulative Distribution Function of';'Expected Non-Catastrophic Debris Pieces'})
        xlabel('Number of Pieces')
        ylabel({'Cumulative Percentage'})
        box on
        grid on
        xlim(xl);
        ylim([0 100]);
        if ~isempty(OutputFolder)
            % Save Figures
            saveas(h,fullfile(OutputFolder,[OutputID '_CollisionConsequence_NonCatastrophicDebrisPieces.fig']));
            saveas(h,fullfile(OutputFolder,[OutputID '_CollisionConsequence_NonCatastrophicDebrisPieces.png']));
        end
        if CloseFigures
            close(h)
        end
    end
    
    %% Generate Number of Pieces Distribution Plots for Catastrophic Collisions
    % Only generate figure if Non-Catastrophic Collision Probability above 5%
    if sum(Catastrophic)/length(Catastrophic)*100 > 5
        h= figure('Position',[75 75 550 700],...
                   'Name','Catastrophic Debris Piece Distribution');
        idx       = find(Catastrophic);
        ExpectedNumberofPieces    = round(10^mean(log10(NumOfPieces(idx))));
        % generate Histogram of distribution
        subplot(2,1,1);
        % define bin edges expecting to represet a chi square distribution
        edges = [min(NumOfPieces(idx)):(median(NumOfPieces(idx))-min(NumOfPieces(idx)))/4:median(NumOfPieces(idx))+(median(NumOfPieces(idx))-min(NumOfPieces(idx)))*6];
        histogram(NumOfPieces(idx),edges)
        title({'Expected Number of Debris Pieces Generated'; 'If a Catastrophic Collision Occurs'})
        xlabel('Number of Pieces')
        ylabel({'Incidence Count ';['(' num2str(length(NumOfPieces),'%.0f') ' Total Trials)']})
        box on
        grid on
        % Add Text Indicating probability of Catastrophic Collision
        yl = ylim;
        xl = xlim;
        text(xl(1)+(xl(2)-xl(1))*0.98,yl(2)*0.97,...
            {[CollisionProbabilityText '% Chance of Catastrophic Collision'];
            [num2str(ExpectedNumberofPieces,'%.0f') ' Expected Debris Pieces '];
             '(Catastrophic Collision)'},...
            'HorizontalAlignment','right',...
            'VerticalAlignment','top',...
            'EdgeColor','k',...
            'BackgroundColor',[1 1 1])
        % Generate CDF of distribution
        subplot(2,1,2);
        hold on
            % Plot CDF
            plot([xl(1) reshape(sort(NumOfPieces(idx),'ascend'),1,length(NumOfPieces(idx)))],...
                [CDFcounts(length(NumOfPieces)-length(idx)+2) CDFcounts(length(NumOfPieces)-length(idx)+2:end)],...
                'LineWidth',2,...
                'Color','b')
        hold off
        title({'Cumulative Distribution Function of';'Expected Catastrophic Debris Pieces'})
        xlabel('Number of Pieces')
        ylabel({'Cumulative Percentage'})
        box on
        grid on
        xlim(xl);
        ylim([0 100]);
        if ~isempty(OutputFolder)
            % Save Figures
            saveas(h,fullfile(OutputFolder,[OutputID '_CollisionConsequence_CatastrophicDebrisPieces.fig']));
            saveas(h,fullfile(OutputFolder,[OutputID '_CollisionConsequence_CatastrophicDebrisPieces.png']));
        end
        if CloseFigures
            close(h)
        end
    end
    
    %% Generate Secondary Object Mass Distribution Plots
    g= figure('Position',[75 75 550 700],...
               'Name','Secondary Mass Distribution');
    % generate Histogram of distribution
    subplot(2,1,1);
    % define bin edges expecting to represet a chi square distribution
    edges = [0:median(massVec)/4:median(massVec)*6];
    histogram(massVec,edges)
    title({'Predicted Secondary Object Mass Distribution'})
    xlabel('Secondary Mass (kg)')
    ylabel({'Incidence Count ';['(' num2str(length(massVec),'%.0f') ' Total Trials)']})
    box on
    grid on
    % Add Text Indicating probability of Catastrophic Collision
    yl = ylim;
    xl = xlim;
    text(xl(2)*0.98,yl(2)*0.97,...
        [CollisionProbabilityText '% Chance of Catastrophic Collision'],...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top',...
        'EdgeColor','k',...
        'BackgroundColor',[1 1 1])
    % Generate CDF of distribution
    subplot(2,1,2);
    hold on
        % Plot CDF
        plot([0 reshape(sort(massVec,'ascend'),1,length(massVec))],...
            [0:100/length(massVec):100],...
            'LineWidth',2,...
            'Color','b')
        % Add Catastrophic Collision Threshold
        if CatastrophicThresholdMass < xl(2)
            plot([CatastrophicThresholdMass CatastrophicThresholdMass],...
                 [0 100],...
                 'LineWidth',2,...
                 'Color','r')
             legend('CDF','Catastrophic Collision Threshold')
        else
            text(xl(2)*0.98,3,...
                {'Catastrophic Collision Threshold';
                 ['(40,000 J/kg): ' num2str(CatastrophicThresholdMass,'%0.2f') ' kg']},...
                'HorizontalAlignment','right',...
                'VerticalAlignment','bottom',...
                'EdgeColor','k',...
                'BackgroundColor',[1 1 1])
        end
    hold off
    title({'Cumulative Distribution Function of';'Secondary Object Mass'})
    xlabel('Secondary Mass (kg)')
    ylabel({'Cumulative Percentage'})
    box on
    grid on
    xlim(xl);
    ylim([0 100]);
    if ~isempty(OutputFolder)
        % Save Figures
        saveas(g,fullfile(OutputFolder,[OutputID '_CollisionConsequence_SecondaryMass.fig']));
        saveas(g,fullfile(OutputFolder,[OutputID '_CollisionConsequence_SecondaryMass.png']));
    end
    if CloseFigures
        close(g)
    end
    Status = 1;
catch
    %% Close Existing Collision Consequence figure windows if open
    try
        close('Secondary Mass Distribution')
    catch
    end
    try
        close('Non-Catastrophic Debris Piece Distribution')
    catch
    end
    try
        close('Catastrophic Debris Piece Distribution')
    catch
    end
    
    warning('An Error occured attempting to generate collision consequence plots')
    Status = 0;
    return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 05-22-2018 | Initial Development
%