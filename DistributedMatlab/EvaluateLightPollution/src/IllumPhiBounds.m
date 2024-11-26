function pmps = IllumPhiBounds(t,a,cts,sts,Re)
% IllumPhiBounds - Illuminated azimithal angle.
% pmps = IllumPhiBounds(t,a,cts,sts,Re)
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
% Illuminated azimuthal angle p = phi bounds function, which gives the
% range of phi for which the LOS intersects the sunlit part of the orbital
% shell with semi-major axis = a. The LOS zenith angle, t = theta, can be
% a scalar or array, intended for use in vectorized integrations. The 
% scalar solar zenith angle, ts, is provided as cts = cos(ts) and 
% sts = sin(ts) to avoid repetitive trig function evaluations.
%
% LOS = line of sight
% SMA = semi-major axis
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Cos and sin of LOS zenith angle, t = theta
ct = cos(t); st = sin(t);

% Range along LOS from observer to orbital shell with semi-major axis = a
ba = sqrt(a^2-(Re*st).^2);
rhoa = ba-Re*ct;

% Cos of phi - phisun, the difference between the terminator azimuth, p,
% and the sun azimuth, ps
cpmps = (-1./st/sts).*(ct*cts+(sqrt(rhoa.*(rhoa+2*Re*ct))+Re*cts)./rhoa);    

% Calculate pmps = phi - phisun using the clipped cos limits
cpmps(cpmps < -1) = -1;
cpmps(cpmps >  1) =  1;
pmps = acos(cpmps);

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