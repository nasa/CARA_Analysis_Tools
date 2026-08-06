function CMBEPH = combine_min_eph(STDEPH,MINEPH,StepFracMin,verbose)
% combine_min_eph - Combine a min. rel. dist. ephemeris with a standard
%                   ephemeris.
%
% Syntax: CMBEPH = combine_min_eph(STDEPH, MINEPH, StepFracMin, verbose);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    STDEPH        -   Standard ephemeris structure with fields:
%                        .T - Time vector                   [1xN]
%                        .X - State vectors                 [6xN]
%                        .P - Covariance matrices           [6x6xN]
%
%    MINEPH        -   Minimum relative distance ephemeris structure
%                      with the same fields as STDEPH
%
%    StepFracMin   -   Fractional time-step threshold below which a
%                      standard ephemeris point is excluded for being
%                      too close to a minimum ephemeris point
%
%    verbose       -   Flag to display command window output
%
% =========================================================================
%
% Output:
%
%   CMBEPH        -   Combined ephemeris structure with fields:
%                       .T - Time vector                    [1xM]
%                       .X - State vectors                  [6xM]
%                       .P - Covariance matrices            [6x6xM]
%                       .N - Number of points in combined 
%                            ephemeris
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

NstdTimes = numel(STDEPH.T); % NminTimes = numel(MINEPH.T);

excl_eph_time = false(size(STDEPH.T)); Nexcl = 0;

idx = (STDEPH.T(1) <= MINEPH.T) & (MINEPH.T <= STDEPH.T(NstdTimes));
idx = find(idx); Nidx = numel(idx);

for ii = 1:Nidx
    
    i = idx(ii);
    
    [c,b,coincidence] = find_bracketing_indices(MINEPH.T(i),STDEPH.T);

    if coincidence
        % Minimum perfectly coincides with an ephemeris point
        excl_eph_time(c) = true; Nexcl = Nexcl+1;
    else
        StepFrac = (MINEPH.T(i)-STDEPH.T(c))/(STDEPH.T(b)-STDEPH.T(c));
        if StepFrac < 0 || StepFrac > 1
            error('Invalid time-step fraction calculated');
        end
        if StepFrac < StepFracMin
            excl_eph_time(c) = true;  Nexcl = Nexcl+1;
        end
    end
        
end

if verbose
    disp([' Number of eph times excluded as too close to minima: ' ...
        num2str(Nexcl)]);
end

% Combine the ephemeris tables, excluding eph points as required
ndx = ~excl_eph_time;    
CMBEPH.T  = cat( 2 , MINEPH.T , STDEPH.T(ndx)     );
CMBEPH.X  = cat( 2 , MINEPH.X , STDEPH.X(:,ndx)   );
CMBEPH.P  = cat( 3 , MINEPH.P , STDEPH.P(:,:,ndx) );

% Sort the combined ephemeris
[CMBEPH.T,ndx] = sort(CMBEPH.T);
CMBEPH.X = CMBEPH.X(:,ndx);
CMBEPH.P = CMBEPH.P(:,:,ndx);

CMBEPH.N = numel(CMBEPH.T);

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================