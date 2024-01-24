function [sig] = RelAtmoDensUnc(altitude_km)

% Rough estimate of the relative atmospheric density uncertainty as a
% function of altitude (from S.Casali's FORTRAN program cdm_pc_cc.f).

x = altitude_km;
x(x <  200) = 200;
x(x > 1000) = 1000;

sig = 0.9644444444E-01 + x .* (-0.9070150220E-03  ...
                       + x .* ( 0.4856002331E-05  ...
                       + x .* (-0.6813649313E-08  ...
                       + x .*   0.2884615385E-11)));

return;
end