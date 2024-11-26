function amass = AirMassRozenberg1966(z,Degrees)
% AirMassRozenberg1966 - Calculate airmass from zenith angle.
% Syntax: amass = AirMassRozenberg1966(z);
%         amass = AirMassRozenberg1966(z,Degrees);
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
% 
% Calculate airmass from zenith angle. Adopted from:
%   Rozenberg, G. V. (1966).  Twilight:  A Study in Atmospheric Optics 
%   (New York: Plenum Press), translated from the Russian by R. B. Rodman, 
%   p. 160.
% 
% =========================================================================
%
% Input:
%
%    z          -   [Scalar or Vector] Zenith angle in radians (default) 
%                   or degrees (optional).
%    Degrees    -   [Boolean] Optional flag indicating if zenith angle is
%                   provided in degrees.
%
% =========================================================================
%
% Output:
%
%   amass       -   [Scalar or Vector] Calculated airmass.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------


if nargin < 2; Degrees = false; end

% Zenith and cosine
if Degrees
    z = z*(pi/180);
end
cosz  = cos(z);

% Do not allow negative cosine values
cosz(cosz < 0) = 0;

% Calculate approximate airmass
amass = 1 ./ ( cosz + 0.025 .* exp(-11 .* cosz) );

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