function [Dvec] = ThreeCompComponent(comp,dat)
% ThreeCompComponent - Generate a design matrix vector for an OCS model component
% Syntax: Dvec = ThreeCompComponent(comp,dat);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved. 
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------


switch lower(comp)
    
    case {'specsphere'}
        f = ones(size(dat.Phase));
        
    case {'diffsphere'}
        p = dat.Phase*pi/180;
        f = (8/3/pi)*(sin(p)+(pi-p).*cos(p));

    case {'diffpabnormalfacet'}
        p = dat.Phase*pi/180;
        f = 2*(1+cos(p));
        
    otherwise
        error('Invalid CalSphereComponent')
        
end

Dvec = f/4/pi;

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