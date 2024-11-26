function I = NvNiNbThetaIntegrand(tht,a,si,clO,slO,cbO,sbO,cbS,sbS,UniformDist)
% NvNiNbThetaIntegrand - Calculate integrand for Nv, Ni, or Nb
% Syntax: I = NvNiNbThetaIntegrand(tht,a,si,clO,slO,cbO,sbO,cbS,sbS,UniformDist);
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
% using the first-order integrated Kessler density approaximation.
%
% Calculation uses coordinate system with Sun mostly along -x axis, so that
%   LocalTimeObs = 0  radians corresponds to midnight
% and 
%   LocalTimeObs = pi radians corresponds to noon
%
% =========================================================================
%
% Input:
% tht = Array for zenith angle of LOS (i.e., axial angle)
% a = Constellation orbital shell semi-major axis (SMA)
% si = Constellation orbital shell sin(inclination)
% (clO,slO) = Cosine & sine of lO = local time of observer
% (cbO,sbO) = Cosine & sine of bO = geocentric latitude of observer
% (cbS,sbS) = Cosine & sine of bS = geocentric latitude of sub-solar point
%             (leave these empty for Nv integrand for which solar 
%              illumination is not required)
% UniformDist = Flag indicating uniform-shell distribution vs Kessler
%               distribution of constellation satellites
%
% =========================================================================
%
% Output:
% I = Integrand function array (same dimension as tht)
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Input array dimension
sizetht = size(tht);

% Initialize output
I = zeros(sizetht);

% Earth radius
Re  = 6378.137;

% Constants
twopi = 2*pi;
    
% Calculate bound in azimuth (phi) angle due to solar illumination
if isempty(cbS) % This indicates solar illumination is not required
    
    % Full phi range when no solar illumination required
    pbeg = zeros(sizetht);
    pend = twopi + pbeg;

else

    % Unit vector to Sun (mostly along -x axis)
    shat = [-cbS; 0; sbS];
    
    % North, east and radial unit vectors
    nhatO = [-clO*sbO; -slO*sbO; cbO];
    ehatO = [-slO; clO; 0];
    rhatO = [clO*cbO; slO*cbO; sbO];
    
    % Azimuth angle of solar unit vector
    pS = atan2(ehatO'*shat,nhatO'*shat);
    
    % Zenith angle of solar unit vector
    ctS = rhatO' * shat;
    stS = sqrt(1-ctS^2);
    
    % Sunlit azimuth bounds
    pdel = IllumPhiBounds(tht,a,ctS,stS,Re);
    pbeg = pS - pdel;
    pend = pS + pdel;
    
end

% Only integrate for zenith angles with finite azimuth (phi) integration
% intervals, i.e., with phi(end) > phi(begin). The remainder have integrals
% equal to zero.
azint = find(pend > pbeg); Nvzint = numel(azint);

if Nvzint > 0
    
    % Select the subset of LOS zenith angles to integrate numerically
    theta = tht(azint); Dtheta = size(theta);
    
    % Sine and cosine of theta arrays
    st = sin(theta); ct = cos(theta);
    
    % Earth position parameters
    % rhatO = [clO*cbO; slO*cbO; sbO];
    % rvecO = Re*rhatO;

    % Constellation SMA squared, etc.
    a2 = a^2; si2 = si^2; Re2 = Re^2;

    % Range function array, i.e., b(a,tht) function
    st2 = st.^2;
    bat = sqrt( a2 - Re2*st2 ); 

    % Range array from observer to constellation orbital shell
    rhoat = bat - Re*ct;
    
    if UniformDist
        
        % Calculate the uniform distribution integrand
        delphi = pend(azint)-pbeg(azint);
        I(azint) = ( st .* rhoat.^2 ./ bat ) .* delphi / (4*a*pi);
        
    else 
        
        % Calculate the Kessler integtrand for the points within bounds
    
        % Calculate phi bounds due to inclination bounds on Kessler dist
        pSouth1 = NaN(Dtheta); pSouth2 = pSouth1; pNorth1 = pSouth1; pNorth2 = pSouth1;

        % ECI z component of observer position
        zO = Re*sbO;

        % Calculate three quantities required for the calculation
        q0 = a*si;
        q1 = (zO+sbO*(rhoat.*ct));
        q2 = cbO*(rhoat.*st);

        % Calculate the azimuth, p, bounds due to inclination limits
        %  expressed using the inequality cpSouth < cos(p) < cpNorth

        % Bounds with cos(p) > cpSouth
        cpSouth = -(q0+q1)./q2;
        ndx = abs(cpSouth) < 1; % Only accept real bounds
        pSouth2(ndx) = acos(cpSouth(ndx));
        pSouth1(ndx) = -pSouth2(ndx);

        % Bounds with cos(phi) < cpNorth
        cpNorth =  (q0-q1)./q2;
        ndx = abs(cpNorth) < 1; % Only accept real bounds
        pNorth1(ndx) = acos(cpNorth(ndx));
        pNorth2(ndx) = -pNorth1(ndx);   

        % Allocate integraion array
        intphi = zeros(size(theta));
        
        % Warning message ID to suppress
        % MSGID = 'MATLAB:integral:NonFiniteValue';

        % for nn=1:Nvzint
        for nn=Nvzint:-1:1

            % Current azimuthal angle (i.e., theta) index
            n = azint(nn);

            % Calculate the  waypoints to use during the phi integration
            PhiWayPoints = [pNorth1(nn),pNorth2(nn),pSouth1(nn),pSouth2(nn)];
            idx = ~isnan(PhiWayPoints);
            PhiWayPoints = PhiWayPoints(idx);
            if any(PhiWayPoints)
                PhiWayPoints = [PhiWayPoints-twopi PhiWayPoints PhiWayPoints+twopi]; %#ok<AGROW>
                idx = pbeg(n) < PhiWayPoints & PhiWayPoints < pend(n);
                PhiWayPoints = unique(PhiWayPoints(idx));
            end

            % Define the integrand for the phi integral
            cbOst = cbO*st(nn);
            sbOct = sbO*ct(nn);
            intfun = @(ppp)real(1./sqrt( ...
                si2 - ((zO+rhoat(nn)*(cbOst*cos(ppp)+sbOct))/a).^2 )); 
            
            % Do the phi integral in intervals, putting singularities at
            % at the interval's beginning and/or ending integration limits
            PhiBounds = [pbeg(n) PhiWayPoints pend(n)];
            M = numel(PhiBounds)-1;
            for m=1:M
                % Integration bounds for this interval
                pa = PhiBounds(m);
                pb = PhiBounds(m+1);
                % Integrand at midpoint
                pm = (pa+pb)/2;
                intmid = intfun(pm);
                % If midpoint is zero, then integral is also zero
                if intmid > 0
                    % Calculate integral with default tolerances
                    % warning('off');
                    intphinm = integral(intfun,pa,pb);
                    % If default tolerances yield bad result, try
                    % using no absolute tolerance as a mitigation
                    if isinf(intphinm) || isnan(intphinm)
                        intphinm = integral(intfun,pa,pb,'AbsTol',Inf);
                        
                        % % If that fails, try adjusting the bounds slightly.
                        % % This happens when an LOS rarely traverses right
                        % % near an inclination singularity point
                        % if isinf(intphinm) || isnan(intphinm)
                        %     % Acceptable adjustment values (smaller is
                        %     % better)
                        %     pEps = 10.^(-6:0.5:-2); NpEps = numel(pEps);
                        %     % Check values 
                        %     iab = intfun([pa pb]);
                        %     adjusta = isinf(iab(1));
                        %     adjustb = isinf(iab(2));
                        %     adjusting = adjusta | adjustb;
                        %     % Adjust until endpoints are both finite, or
                        %     % the acceptable adjustment values expire
                        %     npEps = 1; ppa = pa; ppb = pb; dp = pb-pa;
                        %     while adjusting
                        %         dpEps = dp*pEps(npEps);
                        %         if adjusta; ppa = pa+dpEps; end
                        %         if adjustb; ppb = pb-dpEps; end
                        %         iab = intfun([ppa ppb]);
                        %         adjusta = isinf(iab(1));
                        %         adjustb = isinf(iab(2));
                        %         npEps = npEps+1;
                        %         adjusting = (adjusta | adjustb) & (npEps <= NpEps);
                        %     end
                        %     % Calculate adjust estimate if appropriate
                        %     if ~adjusta && ~adjustb
                        %         intphinm = integral(intfun,ppa,ppb, ...
                        %             'AbsTol',Inf) * (pb-pa)/(ppb-ppa);
                        %         disp(['Adjustment worked: ' num2str(pEps(npEps-1))]);
                        %     else
                        %         disp('Adjustment failed');
                        %     end
                        % else
                        %     disp('Mitigation worked');
                        % end
                        
                        % Turn warnings back on, and report failures
                        % warning('on');
                        % if isinf(intphinm) || isnan(intphinm)
                        %     warning('MATLAB:NvNiNbThetaIntegrand:NonFiniteValue', ...
                        %         'Infinite or NaN encountered despite mitigation/adjustment efforts');
                        % end
                        
                    end
                    % Add this intervals contribution to the integral
                    intphi(nn) = intphi(nn)+intphinm;
                end
            end

        end

        % Calculate the Kessler integtrand for the points within bounds
        I(azint) = ( st .* rhoat.^2 ./ bat .* intphi ) ...
                   / (2*a*pi^2);
        
    end
    
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