function [r,v,P] = StateCovInterp(T, time, Pos, Vel, Cov)
%
% StateCovInterp - Interpolate ECI-frame position and velocity vectors
% and (optionally) the associated ECI-frame covariance matrix.
% Given a time T in days, along with 5 bracketing ephemeris times in days,
% the associated 5 ECI-frame position/velocity vectors each in a 
% 5x3 array, and the 5 associated covariance matrices in
% a nxnx5 array, interpolate the pos and vel vectors at time T using
% Lagrange 5-point interpolation, and the covariance at time T using the
% two-body STM blending method.
%
% Syntax: [r,v,P] = StateCovInterp(T, time, Pos, Vel, Cov)
%
% Inputs:
%    T    -  Time to interpolate to in days, which must be between the first
%            and last epoch times
%    time -  Times in days for the 5 bracketing times (5x1 or 1x5 vector)
%            Assumed to be sorted into monotonically increasing order.
%    Pos  -  Position vectors in ECI coordinates in km (5x3 matrix)
%    Vel  -  Velocity vectors in ECI coordinates in km/s (5x3 matrix)
%    Cov  -  Covariances for five states in (km,s) units (6x6x5 matrix)
%            NOTE: Covariance interpolation is optional, and can be
%            prevented by leaving Cov out of the input variable list or
%            inputing it as an empty set []. Pos/vel vector interpolation
%            is always performed.
%
%    NOTE: This function uses km and km/s units for pos/vel/cov variables. 
%
% Outputs:
%    r -  Interpolated position at T in km        (1x3 vector)
%    v -  Interpolated velocity at T in km/s      (1x3 vector)
%    P -  Interpolated covariance in (km,s) units (6x6 matrix)
%
% Other m-files required:
%   LagrangeStateInterp.m
%   TBSTMCovInterp.m
% Subfunctions: None
% MAT-files required: None
% See also: None
%
% ----------------- BEGIN CODE -----------------

%% Check for sufficient input

Nargin = nargin;
if (Nargin < 4)
    error('Too few input variables');
elseif (Nargin > 5)
    error('Too many input variables');
end

% Ensure correct time array dimensions
sztime = size(time);
if min(sztime) ~= 1 || max(sztime) ~= 5
    error('Incorrect time input variable dimensions');
end

% Ensure correct pos/vel array dimensions
if ~isequal(size(Pos), [5 3]) || ...
   ~isequal(size(Vel), [5 3])
    error('Incorrect Pos or Vel input variable dimensions');
end

% Check for empty covariance, meaning no covariance interpolation is
% required
if Nargin == 4
    Cov = [];
end
isemptyCov = isempty(Cov);

% Check for correct non-empty covariance array dimensions
if ~isemptyCov && ~isequal(size(Cov), [6 6 5])
    error('Incorrect Cov input variable dimensions');
end

% Check if input Covarainces 3X3 or Off-Diagonals=0
if isemptyCov
    UseLagrangeCovInterpolation = false;
else
    if sum(sum(Cov(4:6,4:6,1))) == 0
        % Zero veolcity covariances require Lagrange cov. interpolation
        UseLagrangeCovInterpolation = true;
    else
        % Non-zero veolcity covariances allow two-body state transition
        % matrix (TBSTM) method of cov. interpolation
        UseLagrangeCovInterpolation = false;
    end
end

%% Check to see if T coincides with an ephemeris point. 

% The tolerance 100*eps('double') was chosen here to guard against slight
% differences between the ASW solution for T and O/O-propagated states to T.
% The two aren't typically exactly the same, so this limit allows for
% times to be considered "equal" when they are within some reasonable
% tolerance. For conjunctions between 15 - 20 km/s relative velocity,
% the difference between these values results in ~ 0.03 - 0.04 mm.
dt = T-time;
[dtmin,i] = min(abs(dt));
if dtmin < 100*eps('double')
    % T coincides with time(i)
    r = Pos(i,:);
    v = Vel(i,:);
    if isemptyCov
        P = [];
    else
        P = Cov(:,:,i);
    end    
    return;
end

% Ensure the T is bounded by the input times, otherwise return null values
if T < time(1) || T > time(5)
    r = [];
    v = [];
    P = [];
    return;
end

if UseLagrangeCovInterpolation
    
    warning('Covariance Nodes presented are insufficiently populated for 2-Body Interpolation (required 6X6 fully populated covariance), covariance interpolation will be done using Lagrangian methods (may result in NPD covariances)')
    % Interpolate pos/vel/Cov state using 5-point Lagrange method
    [r,v,P] = LagrangeInterp(T, time, Pos, Vel, Cov);
    
else
    
    % Interpolate pos/vel state using 5-point Lagrange method
    [r,v] = LagrangeStateInterp(T, time, Pos, Vel);
    
    % Check if covariance really needs to be interpolated
    if isemptyCov || (nargout < 3)
        % No covariance interpolation required
        P = [];
        return
    end
    
    % Find the two ephemeris times that bracket time T
    % i = max(find(dt > 0));
    i = find(dt > 0,1,'last');
    if isempty(i) || (i == 5)
        error('Failed to find bracketing times');
    end
    
    % Preceding time quantities
    T1 = time(i);
    X1 = [Pos(i,:) Vel(i,:)]';
    P1 = Cov(:,:,i);
    
    % Following time quantities
    i = i+1;
    T2 = time(i);
    X2 = [Pos(i,:) Vel(i,:)]';
    P2 = Cov(:,:,i);
    
    % Interpolate the covariance using two-body STM blending method
    X = [r v]';
    P = TBSTMCovInterp(T,X,T1,X1,P1,T2,X2,P2);
    
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2020-MAR-31 | Initial Development
% T. Lechtenberg | 2022-Jan-06 | Addition of Catch for 3X3 Covariances to
%                                use lagrange interpolation
% D. Hall        | 2024-Jan-29 | Fixed bug handling empty covariances,
%                                which implies no covariance interpolation
%                                is required. Also fixed bug in the
%                                "UseLagrangeCovInterpolation" branch,
%                                which did not have "Cov" in the input
%                                parameter list to the "LagrangeInterp"
%                                function.