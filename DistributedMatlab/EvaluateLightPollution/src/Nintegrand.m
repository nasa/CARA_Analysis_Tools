function I = Nintegrand(t,ca,sa,Rs,Re,NeedIllumination,HW20function)
% Nintegrand - Calculate integrand
% Syntax: I = Nintegrand(t,ca,sa,Rs,Re,NeedIllumination,HW20function);
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
% Calculate the integrand for the number of satellites above an observatory
% using the uniform-shell apporoximation, normalized to a constellation
% number of Nc = 1.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

if nargin < 7
    HW20function = true;
end

% Integrating this function over zenith angle, theta or t, yields the
% average number of satellites: Na for total above thetamax, Ni for the
% number illuminated by the Sun (at altitude angle a), Ni.

% Also see Hainaut and Williams (2020)

% Constants
Rs2 = Rs^2; Re2 = Re^2;

% Sine and cosine of zenith angle
st = sin(t); % st2 = st.^2;
ct = cos(t); % ct2 = ct.^2;

% Range function
bt = sqrt( Rs2 - Re2*st.^2 ); % b(theta) function

% Range from obs to sat
rho = bt - Re*ct;    
    
% Azimuthal angle for solar illumination

if NeedIllumination
    
    % Calculate the az angle range restricted to illuminated region    
    p0 = NaN(size(t));
    cp0top = sa*(rho.*ct + Re) + sqrt( rho.*(2*ct*Re+rho)  );
    cp0bot = rho.*st*ca;
    cp0 = -cp0top./cp0bot;
    
    % All range of p illuminated
    all_lit = cp0 <= -1;
    p0(all_lit) = pi;
    
    % No range of p illuminated
    all_drk = cp0 >= 1;
    p0(all_drk) = 0;
    
    % Some range of p illuminated
    ndx = ~(all_lit | all_drk);
    p0(ndx) = acos(cp0(ndx));

else
    
    % Az angle extends from -pi <= p <= pi if illimination is not required
    p0 = pi*ones(size(t));
    
end

% Derivative dfsdt

if HW20function

    % Derivative of HW20 eq 10 (ee LightPollution.nb Mathematica file)
    % dfsdt = ( ( 1 - (ct*Re) ./ (Rs*sqrt(1-Re2*st2/Rs2)) ) ...
    %          .* sin( t - asin(Re*st/Rs) ) ) / 2;
    dfsdt = ( ( 1 - (ct*Re) ./ bt ) ...
             .* sin( t - asin(Re*st/Rs) ) ) / 2;
    
    % Integrand function
    I = p0 .* dfsdt / pi;

else
         
    % From CLPUnifShell.nb
    I = (p0 .* st .* rho.^2 ./ bt) / (2*pi*Rs);
         
end

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