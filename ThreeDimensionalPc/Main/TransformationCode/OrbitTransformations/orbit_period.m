function period = orbit_period(r0vec,v0vec,GM)

% Calculate orbital period

r0     = sqrt(sum(r0vec.*r0vec));
v02    = sum(v0vec.*v0vec);
beta   = 2*GM/r0 - v02;

period = (2*pi) * GM * beta^(-1.5);

return;
end