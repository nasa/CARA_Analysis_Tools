function [figSaved, fileName] = SaveConjFigure(fh, conjID, plotType, params)
% SaveConjFigure - Saves a figure associated with a conjuntion ID
%
% Syntax: [figSaved, fileName] = SaveConjFigure(fh, conjID, plotType, params);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Saves the figure depicted in the figure handle passed in according to
%   the conjunction ID and parameters passed in. If no conjunction ID is
%   passed in, then the figure is named according to the
%   params.plot_save_tag parameter.
%
%   The params.plot_save_loc defines the location where the figure is
%   saved. If this parameter is empty, then the figure will not be saved.
%   The file name is constructed from the params.plot_save_loc, conjID, and
%   plotType passed in. All figures will be saved as .png files.
%
% =========================================================================
%
% Input:
%
%   fh - Matlab figure handle with the figure to save.
%
%   conjID - Conjunction ID used as part of the filename. This field can be
%            empty.
%
%   plotType - Character field appended to the end of the file name.
%
%   params - Auxilliary input parameter structure with the following
%            fields:
%
%     plot_save_loc - Defines the save location for the plot. If this
%                     field is empty, then no plot is saved.
%
%     plot_save_ext - (Optional) Extension (a.k.a. file type) for the saved
%                     figure. (default = '.png')
%
%     plot_save_tag - Defines a tag which prepends the file name. This
%                     field can be empty, but must be filled in if the
%                     conjID is empty.
%
% =========================================================================
%
% Output:
%
%   figSaved - Boolean indicating whether or not the figure was saved.
%
%   fileName - Name of the output file if figSaved is true.
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % Return right away if no save location is set, assume that we aren't
    % supposed to save file
    figSaved = false;
    fileName = [];
    if isempty(params.plot_save_loc)
        return;
    end
    
    % Check the conjID and plot_save_tag parameters, both cannot be empty
    if isempty(conjID) && isempty(params.plot_save_tag)
        error('No conjunction ID exists and no params.plot_save_tag is defined, cannot save figure');
    end
    
    % Check file save extension
    if ~isfield(params,'plot_save_ext')
        plotSaveExt = '.png';
    else
        plotSaveExt = params.plot_save_ext;
    end
    
    % Construct the file name
    if ~isempty(conjID)
        fileName = conjID;
        if ~isempty(params.plot_save_tag)
            % Also add the plot_save_tag if it is defined
            if strcmpi(params.plot_save_tag(1),'_')
                % Final tag
                fileName = [fileName params.plot_save_tag];
            else
                % Initial tag
                fileName = [params.plot_save_tag '_' fileName];
            end
        end 
    else
        fileName = params.plot_save_tag;
        if startsWith(fileName,'_')
            fileName = ['UndefinedConj' fileName];
        end
    end
    if isempty(plotType)
        error('The plotType parameter must be filled in, cannot save figure');
    end
    fileName = fullfile(params.plot_save_loc,[fileName '_' plotType plotSaveExt]);
    
    % Make sure the figure is fully rendered before saving the file
    drawnow;
    pause(0.1);
    saveas(fh,fileName);
    figSaved = true;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-23 | Initial Development
% L. Baars       | 2023-Jun-26 | Added fileName as an output parameter.
% D. Hall        | 2023-Aug-07 | Added distinction between an initial tag
%                                vs. a final tag.
% L. Baars       | 2025-Aug-22 | Start file names with 'UndefinedConj' if
%                                conjID is empty and plot_save_tag begins
%                                with '_'. Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
