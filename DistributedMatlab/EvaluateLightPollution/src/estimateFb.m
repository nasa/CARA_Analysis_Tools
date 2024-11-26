function Fb = estimateFb(Mb,Mmd,Flo,Mlo,Fhi,Mhi)
% estimateFb - Estimate fraction brighter than threshold
% Syntax: Fb = estimateFb(Mb,Mmd,Flo,Mlo,Fhi,Mhi);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate fraction brighter than theshold Mb, assuming an 
% asymmetric Gaussian. 
%
% =========================================================================
%
% Input:
%
%    Mb  -  Magnitude threshold
%    
%    Mmd -  Median magnitude plus adjustment
%
%    Flo -  0.05 (percentile)
%
%    Mlo -  5th percentile magnitude plus adjustment
%
%    Fhi -  0.95 (percentile)
%
%    Mhi -  95th percentile magnitude plus adjustment
% 
% =========================================================================
%
% Output:
%
%    Fb - Fraction brighter than threshold Mb  
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Calculate the Gaussian sigma value either above or below the median
sqrt2 = sqrt(2);
if Mmd >= Mb
    sig = (Mlo-Mmd)/sqrt2/erfinv(-1+2*Flo);
else
    sig = (Mhi-Mmd)/sqrt2/erfinv(-1+2*Fhi);
end

% Calculate the fraction of magnitudes in the asym Gaussian that are
% brighter than Mb
Xb = (Mb-Mmd)/sig/sqrt2;
Fb = 0.5*(1+erf(Xb));

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