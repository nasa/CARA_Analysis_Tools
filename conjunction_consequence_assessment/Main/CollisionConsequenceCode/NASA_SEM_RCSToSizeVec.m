function x = NASA_SEM_RCSToSizeVec(z)
%
% This function uses the NASA Size Estimation Model (SEM) to convert
% a vector of normalized RCS values into normalized characteristic size
% values.
%
% INPUT:
%
%   z = Vector of normalized RCS values, z = RCS/lambda^2, which is
%   a dimensionless quantity and not in dB.
%
% OUTPUT:
%
%   x = Vector of normalized size estimates, x = diameter/lambda.
%
% REFERENCE:
%
%   C. L. Stokely et al., "Haystack and HAX Radar Measurements
%   of the Orbital Debris Environment; 2003", JSC-62815, Nov 2006.
%

% Initialize output

x = NaN(size(z));

% Limits of the Mie region in z

zmie1 = 0.001220;
zmie2 = 2.835;

% Find subsets of the input z vector within the optical,
% Rayleigh and Mie regimes

optical = (z > zmie2);

rayleigh = (z < zmie1);

mie = ~(optical | rayleigh);
    
% Process all z values in optical regime

if any(optical)
    x(optical) = sqrt(4*z(optical)/pi);
end

% Process all z values in Mie regime

if any(mie)
    
    % Interpolate from table of the Mie region values
    
    % Table 1 of JSC-62815

    xztab = [0.10997    0.001220  ; ...
             0.11685    0.001735  ; ...
             0.12444    0.002468  ; ...
             0.13302    0.003511  ; ...
             0.14256    0.004993  ; ...
             0.15256    0.007102  ; ...
             0.16220    0.01010   ; ...
             0.17138    0.01437   ; ...
             0.18039    0.02044   ; ...
             0.18982    0.02907   ; ...
             0.20014    0.04135   ; ...
             0.21237    0.05881   ; ...
             0.22902    0.08365   ; ...
             0.25574    0.1190    ; ...
             0.30537    0.1692    ; ...
             0.42028    0.2407    ; ...
             0.56287    0.3424    ; ...
             0.71108    0.4870    ; ...
             0.86714    0.6927    ; ...
             1.0529     0.9852    ; ...
             1.2790     1.401     ; ...
             1.5661     1.993     ; ...
             1.8975     2.835    ];
    
    % Estimate size using logarithmic interpolation
    
    x(mie) = exp(interp1(log(xztab(:,2)),log(xztab(:,1)),log(z(mie)))); 
             
end

% Process all z values in Rayleigh regime

if any(rayleigh)
    x(rayleigh) = (4*z(rayleigh)/(9*pi^5)).^(1/6);
end

% Check for undefined values

if any(isnan(x))
    error('Some size values remain undefined');
end

end