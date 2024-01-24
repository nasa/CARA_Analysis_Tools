function period = orbit_period(r0vec,v0vec,GM)

% Calculate orbital period

% If GM value is not input, assume m and s units
if nargin < 3
    % GM = Grav. constant times central mass, the EGM-96 value in m^3/s^2
    GM = 3.986004418e14;
end

% Calculate the period assuming two-body motion
r0     = sqrt(sum(r0vec.*r0vec));
v02    = sum(v0vec.*v0vec);
beta   = 2*GM/r0 - v02;
period = (2*pi) * GM * beta^(-1.5);

return;
end