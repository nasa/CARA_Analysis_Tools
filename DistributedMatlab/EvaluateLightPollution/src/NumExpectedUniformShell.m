function [Nacon,Nicon,Nbcon] = NumExpectedUniformShell(Nc,hkm,zenmaxdeg,zaint,fbint,asundeg,HW20,verbose)
% NumExpectedUniformShell - Estimate satellites above using uniform shell approx
% Syntax: [Nacon,Nicon,Nbcon] = NumExpectedUniformShell(Nc,hkm,zenmaxdeg,zaint,fbint,asundeg,HW20,verbose);
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
% Use uniform shell approx to estimate the number of constellation
% satellites above an Earth-based observer. Also estimate the normalized
% flux from those satellites
%
% =========================================================================
%
% Input:
%  Nc = Number of constellation satellites
%  h = Constellation altitude
%  zenmaxdeg = Max observation zenith angle
%  asundeg = Solar elevation (i.e., altitude) angle, which is the negative
%            of the solar depression angle
%  zaint = Zenith angle interpolation array (rad) 
%          for brighter-than-recommended fraction
%  fbint = Fraction array for brighter-than-recommended fraction
%  HW20 = Flag to use Hainaut & Williams (2020) integrand
%         expression, which produces equivalent results as
%         Hall (2021) integrand
%  verbose = Verbose flag
%
% =========================================================================
%
% Output:
%  Nacon = Number of constellation satellites above a low-lat. observer
%          within specified angle of the zenith direction
%  Nicon = Number of those satellites that are also illumiated by the Sun
%  Nbcon = Number of those satellites that are also brighter than
%          recommended
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Initializations
Nargin = nargin;

na = 0;
na=na+1; if Nargin < na; Nc = 1600; end
na=na+1; if Nargin < na; hkm = 550; end
na=na+1; if Nargin < na; zenmaxdeg = 60; end
na=na+1; if Nargin < na; asundeg = -18; end
na=na+1; if Nargin < na; zaint = [0 pi/2]; end
na=na+1; if Nargin < na; fbint = [1 0]; end
na=na+1; if Nargin < na; HW20 = false; end
na=na+1; if Nargin < na; verbose = false; end

% Solar altitude angle (e.g., end of Ast. twilight = 18 deg)
a = asundeg * pi/180; ca = cos(a); sa = sin(a);

% Use t = theta for zenith angle
tmax = zenmaxdeg * pi/180;
ctmax = cos(tmax); stmax = sin(tmax);

% Earth radius
Re = 6378.137;

% Constellation altitude and radius
h = hkm;
Rs = Re+h;

% Check if asun < 0
if (a > 0)
    error('Estimation for a > 0 not tested');
end

% Calculate the total number of satellites above zenith cutoff 
NeedIllumination = false;   
Naint = @(t)Nintegrand(t,ca,sa,Rs,Re,NeedIllumination);
Nacon = Nc*integral(Naint,0,tmax);

% Analytical form by Hall (2021)
% Nacon1 = (Nc*(Rs - Re*stmax^2 - ctmax*sqrt(Rs^2 - Re^2*stmax^2)))/(2*Rs);
btmax = sqrt(Rs^2 - Re^2*stmax^2);
Nacon1 = (Nc/2/Rs) * (Rs - Re*stmax^2 - ctmax*btmax);

% Analytical form from HW20
fncz = asin(Re*stmax/Rs);
fncz = cos(tmax-fncz);
fncz = 0.5*(1-fncz);
Nacon2 = Nc * fncz;

if verbose
    Nfig = 15;
    disp(['Nacon = ' smart_exp_format(Nacon ,Nfig) ' ' ...
                     smart_exp_format(Nacon1,Nfig) ' ' ...
                     smart_exp_format(Nacon2,Nfig)]);
end

% Calculate the number of illuminated satellites above zenith cutoff 
NeedIllumination = true;   
Niint = @(t)Nintegrand(t,ca,sa,Rs,Re,NeedIllumination,HW20);
Nicon = Nc*integral(Niint,0,tmax);

if verbose
    disp(['Nicon = ' smart_exp_format(Nicon,Nfig)]);
end

Nbint = @(t)Nintegrand(t,ca,sa,Rs,Re,NeedIllumination,HW20).*interp1(zaint,fbint,t);
Nbcon = Nc*integral(Nbint,0,tmax);

if verbose
    disp(['Nbcon = ' smart_exp_format(Nbcon,Nfig)]);
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