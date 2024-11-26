function [ rsun ] = SunPos(JDut1)
% SunPos - Computes the sun-earth vector using the analytical formulation
% shown in the Astronomical Almanac (1992), producing a J2000 MEME vector
% with an accuracy of approximately 0.01 degree
% Syntax: rsun = SunPos(JDut1);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%    JDut1 - Julian Date (UT1)
%
% =========================================================================
%
% Output:
%    rsun - Sun-Earth vector in J2000 MEME
%
% =========================================================================
%
% References:
%   Algorithm from Fundamentals of Astrodynamics and Applications by 
%   David A. Vallado
% 
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------


    % AU to km conversion factor
    AU2KM = 149597870.700;

    % Computing Julian Centuries (UT1)
    Tut1 = JulianCenturies(JDut1);
    
    % Mean longitude of the sun (degs)
    lambdaMSun = 280.4606184 + 36000.77005361 .* Tut1;
    
    % ASSUMPTION: TDB approx. equals Tut1 - (SHOULD BE TT approx equals Tut1?)
    Ttdb = Tut1;
    
    % Mean anomaly of the sun
    Msun = 357.5277233 + 35999.05034 .* Ttdb;
    
    % Ecliptic longitude of the sun (the ecliptic latitude never exceeds
    % 0.000333 degrees and is usually assumed to be 0)
    lambdaEcliptic = lambdaMSun + 1.914666471 .* sind(Msun) + 0.019994643 .* sind(2 .* Msun);
    
    % Distance from the earth to the sun in AU
    r = 1.000140612 - 0.016708617 .* cosd(Msun) - 0.000139589 .* cosd(2 .* Msun);
    
    % Obliquity of the ecliptic
    eps = 23.439291 - 0.0130042 .* Ttdb;
    
    % Computing J2000 MEME vector components and converting to km
    rsun(:,1) = r .* cosd(lambdaEcliptic) .* AU2KM;
    rsun(:,2) = r .* cosd(eps) .* sind(lambdaEcliptic) .* AU2KM;
    rsun(:,3) = r .* sind(eps) .* sind(lambdaEcliptic) .* AU2KM;
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY --------------------
% Developer      |    Date     |     Description
%----------------------------------------------------
% N. Sabey       | Sept - 2014 |  Initial Development
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================