function I = NvNiNbThetaPhiIntegrand(tht,phi,a,si2,clO,slO,cbO,sbO,cbS,sbS,ZGrid,FGrid,UniformDist)

% NvNiNbThetaPhiIntegrand - Calculate integrand for Nv, Ni, or Nb.
% Syntax: I = NvNiNbThetaPhiIntegrand(tht,phi,a,si2,clO,slO,cbO,sbO,cbS,sbS,ZGrid,FGrid,UniformDist);
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
% using the first-order integrated Kessler density approaximation
%
% =========================================================================
%
% Input:
% phi = Array for azimuth of line of sight (LOS) in radians
% tht = Array for zenith angle of LOS (i.e., axial angle) in radians
% a = Constellation orbital shell semi-major axis (SMA) in km
% si2 = Constellation orbital shell sin(inclination) squared
% (clO,slO) = Cosine & sine of lO = local time of observer
% (cbO,sbO) = Cosine & sine of bO = geocentric latitude of observer
% (cbS,sbS) = Cosine & sine of bS = geocentric latitude of sub-solar point
%             Leave these arrays empty if solar illumination factor is not
%             required (i.e., for Nv calculations)
% (ZGrid,FGrid) = Zenith angle and brighter-than-recommended fraction
%                 interpolation grids. 
%             Leave these arrays empty if brighter-than-recommended factor
%             is not required (i.e., for Nv or Ni calculations)
% UniformDist = Flag indicating uniform-shell distribution vs Kessler
%               distribution of constellation satellites
%
% =========================================================================
%
% Output:
% I = Integrand function for Nv, Ni or Nb (theta,phi) integration
%
% Calculation uses coordinate system with Sun mostly along -x axis, so that
%   LocalTimeObs = 0  radians corresponds to midnight
% and 
%   LocalTimeObs = pi radians corresponds to noon
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Convert input (phi,tht) arrays into row vectors
OrigSize = size(phi);
Nptvec = numel(phi);
Dptvec = [1 Nptvec];
phi = reshape(phi,Dptvec);
tht = reshape(tht,Dptvec);

% Initialize integrand output
I = zeros(Dptvec);

% Ground-based observer position parameters
Re  = 6378.137; Re2 = Re^2;
rhatO = [clO*cbO; slO*cbO; sbO];
rvecO = Re*rhatO;

% Constellation shell squared semi-major axis
a2 = a^2; 

% Sine and cosine of LOS azimuthal angle arrays
sp = sin(phi); cp = cos(phi);

% Sine and cosine of LOS zenith angle arrays
st = sin(tht); ct = cos(tht);

% Range function array, i.e., bat = b(a,tht) function
st2 = st.^2;
bat = sqrt( a2 - Re2*st2 ); 

% Range array from observer to constellation orbital shell, rho(a,tht)
rhoat = bat - Re*ct;

% North and each unit vectors used for azimunth and zenith (phi,theta)
% LOS angles; the set (nhatO,ehatO,rhatO) forms a right-handed triad of
% unit vectors
nhatO = [-clO*sbO; -slO*sbO; cbO];
ehatO = [-slO; clO; 0];

% Position vectors of LOS intersections with constellation orbital shell,
% rvecapt = rvec(a,phi,theta)
strep = repmat(st,[3 1]);
rvecapt = repmat(rvecO,Dptvec)  + repmat(rhoat,[3 1]) .* ( ...
          repmat(nhatO,Dptvec) .* repmat(cp,[3 1]).*strep + ...
          repmat(ehatO,Dptvec) .* repmat(sp,[3 1]).*strep + ...
          repmat(rhatO,Dptvec) .* repmat(ct,[3 1]) );

% Find LOS intersections that are within  inclination bounds
if UniformDist
    % All within inclination bounds for uniform shell distribution
    inb = true(Dptvec);
else
    % Latitudes of LOS intersections with constellation orbital shell
    % sbapth = sin(beta(a,phi,theta)) 
    sbapt2 = rvecapt(3,:).^2 / a2;
    % Find LOS intersections that are within the inclination bounds
    inb = sbapt2 < si2;
end

% Number of shell-LOS intersection points within inclination bounds
Ninb = sum(inb);

% Process points that are in bounds
if Ninb ~= 0

    % Calculate the integrand for Nv - the total number of satellites above
    % the observer
    if UniformDist
        % Calculate the uniform distribution integrand
        I(inb) = ( st(inb) .* rhoat(inb).^2 ./ bat(inb) ) / (4*a*pi);
    else 
        % Calculate the Kessler integtrand for the points within bounds
        I(inb) = ( st(inb) .* rhoat(inb).^2 ./ bat(inb)  ...
                   ./ sqrt(si2-sbapt2(inb)) ) / (2*a*pi^2);
    end

    % Calculate and apply the solar illumination factors, if required
    if ~isempty(cbS)

        % Unit vector to Sun, mostly along -x axis because the midnight
        % plane (from which local time) is along the +x axis by definition
        shat = [-cbS; 0; sbS];

        % Close approach ranges from LOS shell intersections towards Sun
        shatrep = repmat(shat,[1 Ninb]);
        rhopca = -sum( rvecapt(:,inb) .* shatrep , 1 );

        % Close approach points
        rvecca = rvecapt(:,inb) + shatrep.*repmat(rhopca,[3 1]);

        % Illumination is zero when CA points are towards Sun and 
        % within Earth sphere
        r2ca = sum( rvecca .* rvecca , 1 );    
        Illum = rhopca < 0 | r2ca >= Re2;

        % Inlcude illumination factors in integrand array
        I(inb) = I(inb) .* Illum;   

    end
    
    % Multiply by brighter-than-recommended fraction of satellites, if
    % required, interpolated from grid of input values
    if ~isempty(ZGrid)
        I(inb) = I(inb) .* interp1(ZGrid,FGrid,tht(inb));
    end
    
end

% Reshape output integrand array to the original size
I = reshape(I,OrigSize);

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