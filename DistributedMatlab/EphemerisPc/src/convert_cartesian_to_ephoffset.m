function [EOt,out] = convert_cartesian_to_ephoffset(t,Xt,MeanXt,eph,verbose)
% convert_cartesian_to_ephoffset - Convert cartesian state to ephemeris-
%                                  offset (EO) state.
%
% Syntax: [EOt, out] = convert_cartesian_to_ephoffset(t, Xt, MeanXt, eph);
%         [EOt, out] = convert_cartesian_to_ephoffset(t, Xt, MeanXt, eph, verbose);
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
%    t          -   Time (eph.T units, e.g., days)
%
%    Xt         -   Cartesian state at time t               [6x1]
%                   (eph.X units, e.g., km and km/s)
%
%    MeanXt     -   Mean trajectory cartesian state         [6x1]
%                   (same units as Xt)
%
%    eph        -   Ephemeris structure with fields:
%                     .T             - Ephemeris times      [1xN]
%                     .X             - Mean states          [6xN]
%                     .N             - Number of points
%                     .SecPerTimeUnit - Seconds per time
%                                      unit (optional,
%                                      default = 86400)
%
%    verbose    -   Verbosity level: 0, 1, 2, or 3
%                   (optional, default = 1)
%
% =========================================================================
%
% Output:
%
%   EOt         -   Eph. offset state transformed from Xt    [6x1]
%                   (km and km/s)
%
%   out         -   Auxiliary output
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initializations
Nargin = nargin;
if Nargin < 5; verbose = []; end
if isempty(verbose); verbose = 1; end % Levels 0, 1, 2, or 3

% Persistent parameters
persistent ExtrapPars fms_options

% Define the ephemeris extrapolation parameters
if isempty(ExtrapPars)
    % Issue a warning for extrapolations beyond first number of eph.
    % points, but issue no errors because eph. offset coordinates can
    % sometimes require state extrapolations
    ExtrapPars = [Inf Inf];
end

% Get eph. time units
if ~isfield(eph,'SecPerTimeUnit') || isempty(eph.SecPerTimeUnit)
    SecPerTimeUnit = 86400; % Assume time units of days
else
    SecPerTimeUnit = eph.SecPerTimeUnit;
end

% Extract pos and vel
rt = Xt(1:3); vt = Xt(4:6);

% Extract mean pos and vel
rbart = MeanXt(1:3); vbart = MeanXt(4:6);

% Squared mean vel
% vbart2 = vbart' * vbart;

% Calculate the time offset
if isequal(Xt,MeanXt)
    
    % If field point is equal to the mean-traj point, then there is no time
    % offset
    tau = t; rbartau = rbart; vbartau = vbart; converged = 1;
    
else
    
    % Initial estimate of offset time
    t0 = t; r0 = rbart; v0 = vbart;
    
    % Begin and end times for the ephemeris
    TephBeg = eph.T(1); TephEnd = eph.T(eph.N);
    
    % Iteratively refine offset time
    iterating = true; converged = 0; adtmin = Inf;
    iter = 0; itermax = 10; drtol = 1e-10; dr2best = Inf;
    while iterating
        % Estimate offset time assuming rectilinear motion from current
        % estimate (r0,v0)
        iter = iter+1;
        dt1 = ( v0' * (rt-r0) ) / ( v0' * v0 ); % Increment in seconds
        adt1 = abs(dt1);
        t1 = t0 + dt1/SecPerTimeUnit;
        % Interpolated ephemeris state at offset time, t1, allowing state
        % extrapolation at the edges
        [r1,v1] = interpState(t1,eph,[],ExtrapPars); r1 = r1'; v1 = v1';
        if isempty(r1)
            iterating = false; converged = -1;
        else
            if verbose > 2
                % Set the extrapolation flag
                if t1 < TephBeg
                    Extrap = -1;
                elseif t1 > TephEnd
                    Extrap = 1;
                else
                    Extrap = 0;
                end
                % Report results
                disp([' iter = ' num2str(iter) ...
                    ' dt = ' smart_exp_format(dt1,4) ...
                    ' dr = ' smart_exp_format(norm(r1-r0),4) ...
                    ' t = ' smart_exp_format(t1,8) ...
                    ' Extrap = ' num2str(Extrap)]);
            end
            % Save the best estimate for tau
            dr = r1-rt;
            dr2 = dr' * dr;
            if dr2 < dr2best
                dr2best = dr2; t1best = t1; r1best = r1; v1best = v1;
            end
            % Check for convergence
            if norm(r1-r0) <= drtol || adt1 > adtmin 
                iterating = false; converged = 1;
            elseif iter >= itermax 
                iterating = false;
            else
                t0 = t1; r0 = r1; v0 = v1; adtmin = min(adtmin,adt1);
            end
        end
    end
    
    % Process converged results
    if converged >= 0
        
        % Adopt the best offset time
        tau = t1best; rbartau = r1best; vbartau = v1best;
        
        % % Test tau solution with Matlab fminsearch function
        if verbose > 3 % || t > 1
            t0 = t;
            % Set options for fminsearch function
            if isempty(fms_options)
                fms_options = optimset('fminsearch');
                fms_options.Display = 'iter';
                fms_options.TolX = 1e-18 / 86400;
                fms_options.TolFun = 1e-18;
            end
            % Use Matlab function to find refined offset time
            drfun = @(tt) sqrt(sum((rt-interpState(tt,eph,[],ExtrapPars)').^2));
            t1 = fminsearch(drfun,t0,fms_options);
            [r1,v1] = interpState(t1,eph,[],ExtrapPars); r1 = r1'; v1 = v1';
            [t1best t0 t1best-t0]*SecPerTimeUnit
            [t1best t1 t1best-t1]*SecPerTimeUnit
            [r1best r1 r1best-r1]
            [v1best v1 v1best-v1]
            keyboard;
        end
        
    end
    
end

out.converged = converged;
if converged < 0
    
    if verbose > 0
        warning('Fatally unconverged ephemeris offset time');
    end
    EOt = NaN(6,1); out = [];
    
else
    
    if converged == 0 && verbose > 1
        warning('Unconverged ephemeris offset time');
    end
    
    % Calculate the along-traj. vector (like V of VNB)
    Xhat_tau = vbartau'/norm(vbartau);
    % Calculate first traj.-normal vector (like N of VNB)
    Yhat_tau = cross(rbartau,vbartau); Yhat_tau = Yhat_tau'/norm(Yhat_tau);
    % Calculate second traj.normal vector (like B of VNB)
    Zhat_tau = cross(Xhat_tau,Yhat_tau);
    
    % Calculate the ephemeris-offset "position" state
    dt = (tau-t)*SecPerTimeUnit;
    dr = rt-rbartau;
    
    % Eph. offset velocity state
    dv = vt-vbart;
    
    % Eph. offset full state vector
    EOt = [norm(vbart)*dt ...
           Yhat_tau*dr ...
           Zhat_tau*dr ...
           Xhat_tau*dv ...
           Yhat_tau*dv ...
           Zhat_tau*dv]';

    % Other output
    out.tau = tau;
    out.Xbar_tau = [rbartau; vbartau];
    out.Xhat_tau = Xhat_tau;
    out.Yhat_tau = Yhat_tau;
    out.Zhat_tau = Zhat_tau;
      
end

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