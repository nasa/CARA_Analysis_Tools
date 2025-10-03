function h = figure_set_up(wid_to_hgt_ratio, visible)
% figure_set_up - Creates a new figure with the aspect ratio passed in.
%
% Syntax: h = figure_set_up(wid_to_hgt_ratio, visible);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    wid_to_hgt_ratio - Floating point value depicting the aspect ratio.
%    visible          - (Optional) Boolean which determines if the figure
%                       should be displayed.
%                       Defaults to true.
%
% =========================================================================
%
% Output:
%
%    h - Figure handle
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % Determine the screen to display the figure on (largest screen)
    screenInfo = get(0,'MonitorPositions');
    numScreens = size(screenInfo,1);
    largestScreen = numScreens;
    largestScreenDim = screenInfo(largestScreen,[3 4]);
    for i = largestScreen-1:-1:1
        if screenInfo(i,3)*screenInfo(i,4) > largestScreenDim(1)*largestScreenDim(2)
            largestScreen = i;
            largestScreenDim = screenInfo(i,[3 4]);
        end
    end
    
    % Get the screen size and position from the bottom left corner of the
    % screen
    scrsz = screenInfo(largestScreen,:);
    sclft = scrsz(1) + round(0.05*scrsz(3));
    scbot = scrsz(2) + round(0.05*scrsz(4));

    % Make the max extent of the figure 85% of the smaller of the screen
    % width or screen height. Size the other axis according to the aspect
    % ratio.
    if (scrsz(3) > scrsz(4))
        schgt = 0.85*scrsz(4);
        scwid = schgt*wid_to_hgt_ratio;
    else
        scwid = 0.85*scrsz(3);
        schgt = scwid/wid_to_hgt_ratio;
    end

    % Set the visibility parameter
    if (nargin < 2) || isequal(visible,true)
        visible = 'on';
    elseif isequal(visible,false)
        visible = 'off';
    end

    % Create the figure with some default values standard to all CARA
    % figures
    h = figure('Units','pixels',                             ...
               'Renderer','zbuffer',                         ...
               'Toolbar','figure',                           ...
               'Position', [sclft scbot scwid schgt],        ...
               'PaperPositionMode','auto',                   ...
               'MenuBar', 'figure',                          ...
               'NumberTitle','off',                          ...
               'Visible',visible,                            ...
               'Color', [1 1 1],                             ...
               'InvertHardcopy','off');

    return;       
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-20 | Initial Development
% L. Baars       | 2025-Aug-25 | Updated code for public release. Added
%                                code to search for and display the figure
%                                on the largest screen.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================