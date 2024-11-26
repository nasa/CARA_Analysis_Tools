function T = JulianCenturies(JD)
% JulianCenturies - Computes the Julian Century value from a Julian Date
% Syntax: T = JulianCenturies( JD )
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
% 
% Algorithm Source: IAU Standards of Fundamental Astronomy (SOFA)
%
% =========================================================================
%
% Input:
%    JD -  Julian Date  
%
% =========================================================================
%
% Output:
%    T -  Julian Century
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------
    
    T = (JD - 2451545) ./ 36525;

end

% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% N. Sabey       | Aug - 2014 |  Initial Development
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================