function [dTCA,X1CA,X2CA] = FindNearbyCA(X1,X2,MotionMode,RelTol)
% FindNearbyCA - Find the close-approach (CA) point from an input primary
%                and secondary position/velocity inertial state.
%
% Syntax: [TCA,X1CA,X2CA] = FindCA(X1,X2,MotionMode,RelTol)
%
%    X1      - Primary object's pos/vel state vector in ECI coordinates
%              (6x1) [m & m/s] or [km & km/s] 
%    X2      - Secondary object's pos/vel state vector in ECI coordinates
%              (6x1) [m & m/s] or [km & km/s] 
%    MotionMode - 'LINEAR' (currently implemented)
%               - 'TWOBODY' (possible future addition)
%    RelTol  - Tolerance for finding CA (reserved for TWOBODY MotionMode)
%
% Outputs:
%    dTCA    - Offset time to CA (s)
%    X1CA    - Primary  ECI state at CA (6x1) [m & m/s] or [km & km/s] 
%    X2CA    - Seconary ECI state at CA (6x1) [m & m/s] or [km & km/s] 
%
% Example/Validation Cases:
%    Conjunctions with unacceptably large TCA offsets:
%     25544_conj_44437_20190727_160507_20190727_160507 (dTCA = 497.068 s)
%     43042_conj_43043_20190505_061624_20190428_061624 (dTCA = 28.329 s)
%     25544_conj_44437_20190803_160507_20190727_160507 (dTCA = 27.198 s)
%     26998_conj_81790_20180707_061819_20180707_061819 (dTCA = 0.173 s)
%    Conjunctions with negligbly small TCA offsets:
%     40115_conj_24925_20181003_073833_20180928_144418 (dTCA = 9.91E-06 s)
%     43689_conj_31620_20190725_050922_20190719_152452 (dTCA = 1.92E-08 s)
%     
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Last revision: 2018-SEP-04
%
% ----------------- BEGIN CODE -----------------

    % Set up defaults
    
    Nargin = nargin;

    if Nargin < 3 || isempty(MotionMode)
        MotionMode = 'LINEAR';
    else
        MotionMode = upper(MotionMode);
    end
    
    % Initialize motion mode
    
    switch MotionMode
        case 'LINEAR'
            % Linear motion mode
            motion = 1;
        case 'TWOBODY'
            % 2-body motion mode
            motion = 2;
            % Get relative tolerance for convergence
            if Nargin < 4 || isempty(RelTol)
                RelTol = 1E-12; %#ok<NASGU>
            end
        otherwise
            error('Invalid motion mode');
    end

    % Handle different motion modes

    if motion == 1 % Linear motion

        % Primary and secondary cartesian positions and velocities
        r1 = X1(1:3);
        v1 = X1(4:6);
        r2 = X2(1:3);
        v2 = X2(4:6);

        % Relative velocity
        v = v2-v1;
        vmag2 = v'*v;

        % Handle zero relative velocity case
        if vmag2 == 0
            % No TCA offset can be calculated for zero relative velocity
            dTCA = NaN;
            % Primary and secondary CA positions remain unchanged for the
            % case of linear motion with zero relative velocity
            r1CA = r1;
            r2CA = r2;
        else
            % TCA offset for linear relative motion
            dTCA = -((r2-r1)'*v)/vmag2;
            % Primary and secondary positions at the linear-motion TCA
            r1CA = r1+dTCA*v1;
            r2CA = r2+dTCA*v2;
        end

        % Return final linear motion states
        
        X1CA = [r1CA; v1];
        X2CA = [r2CA; v2];
        
    elseif motion == 2 % Two-body motion
        
        error('TWOBODY motion mode not implemented');
        
    end
    
    return;

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2019-AUG-07 | Initial Development
% D. Hall        | 2019-SEP-04 | Expanded comments and description
%

end