function [J] = J_CART2EQ_Analytic(State)
%
% J_CART2EQ_Analytic - Analytically computes the Jacobian matrix using partial derivatives.
%                      This Jacobian matrix is used to either convert 
%                      Cartesian ECI covariances into Equinoctial element set.
%
% Syntax:            [J] = Jacobian_CART2EQ(State)
%
% Inputs:
%    State  - State vector in Cartesian ECI coordinates.         Units: km and km/s
%             (1x6 row vector)      
%
% Outputs:
%    J - Jacobian of Transformation Matrix (6x6 matrix)
    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu  = 3.986004418e5;

   
    % Convert ECI state to Keplerian elements (Radians, km)
    [KEP] = Cart2Kep(State,'Eccentric','Rad');
    a     = KEP(1);
    e     = KEP(2);
    i     = KEP(3);
    Omega = KEP(4);
    w     = KEP(5);
    EA    = KEP(6);
    
    % Compute mean anomaly from eccentric anomaly
    % M = EA - e * sin(EA);

    % Equinoctial element set in terms of Keplerian element set
    [EQ] = Kep2Eq(KEP);
    Af  = EQ(1); %e * cos(w + Omega);
    Ag  = EQ(2); %e * sin(w + Omega);
    % L = EQ(3); %mod(Omega + w + M, 2*pi);  % Mean Longitude [Units: rad]
    n   = EQ(4); %sqrt(mu/a^3);              % Mean Motion    [Units: 1/s]
    chi = EQ(5); %tan(i/2) * sin(Omega);
    psi = EQ(6); %tan(i/2) * cos(Omega);
    
    % Get Retrograde Factor
    fr = 1;
    if abs(i) > pi/2 && abs(i)<3*pi/2
        fr = -1;
    end
    
    % Jacobian matrix between Keplerian and Equinoctial elements
    % (Matrix od partial derivatives, dKEP/dEQ)
    J_KEP_EQ      = zeros(6,6);
    J_KEP_EQ(1,4) = (-2/(3*n))*nthroot(mu/n^2,3);
    J_KEP_EQ(2,1) = Af / sqrt(Af^2 + Ag^2);
    J_KEP_EQ(2,2) = Ag / sqrt(Af^2 + Ag^2);
    J_KEP_EQ(3,5) = 2*fr*chi / ((1 + chi^2 + psi^2)*sqrt(chi^2 + psi^2));
    J_KEP_EQ(3,6) = 2*fr*psi / ((1 + chi^2 + psi^2)*sqrt(chi^2 + psi^2));
    J_KEP_EQ(4,5) =  psi / (chi^2 + psi^2);
    J_KEP_EQ(4,6) = -chi / (chi^2 + psi^2);
    J_KEP_EQ(5,1) = -Ag / (Af^2 + Ag^2);
    J_KEP_EQ(5,2) =  Af / (Af^2 + Ag^2);
    J_KEP_EQ(5,5) = -fr*psi / (chi^2 + psi^2);
    J_KEP_EQ(5,6) =  fr*chi / (chi^2 + psi^2);
    J_KEP_EQ(6,1) =  Ag / (Af^2 + Ag^2);
    J_KEP_EQ(6,2) = -Af / (Af^2 + Ag^2);
    J_KEP_EQ(6,3) = 1; 

    % Convert Keplerian elements to Perifocal coordinates via eccentric anomaly
    r     = a * (1 - e*cos(EA));
    EAdot = sqrt(mu/a) / r;
    xp    = a * (cos(EA) - e);
    yp    = a * sin(EA) * sqrt(1-e^2);
    zp    = 0;
    xdotp = -a * sin(EA) * EAdot;
    ydotp =  a * sqrt(1-e^2) * cos(EA) * EAdot;
    zdotp = 0;
    
    % Rotation matrix from Perifocal to ECI Cartesian coordinates 
    PQW2ECI = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i), ...
               -1*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i), ...
               sin(Omega)*sin(i); ...
               sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i), ...
               -1*sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i), ...
               -1*cos(Omega)*sin(i); ...
               sin(w)*sin(i), cos(w)*sin(i), cos(i)];
           
    % Compute auxiliary partial derivatives necessary for assembling
    % Jacobian matrix, dECI/dEQ, between Cartesian ECI and Keplerian elements
    
    % Partial derivative of xp w.r.t. a,e,M
    dxp_da = xp/a;
    dxp_de = -a - (yp^2)/(r*(1-e^2));
    dxp_dM = xdotp / n;
    
    % Partial derivative of yp w.r.t. a,e,M
    dyp_da = yp/a;
    dyp_de = (xp*yp)/(r*(1-e^2));
    dyp_dM = ydotp / n;
    
    % Partial derivative of zp w.r.t. a,e,M
    dzp_da = 0;
    dzp_de = 0;
    dzp_dM = 0;
    
    % Partial derivative of xdotp w.r.t. a,e,M
    dxdotp_da = -xdotp/(2*a);
    dxdotp_de = xdotp*(a/r)^2*(2*(xp/a) + e/(1-e^2)*(yp/a)^2);
    dxdotp_dM = -n * xp * (a/r)^3;
    
    % Partial derivative of ydotp w.r.t. a,e,M
    dydotp_da = -ydotp/(2*a);
    dydotp_de = n/(sqrt(1-e^2))*(a/r)^2*(xp^2/r - yp^2/(a*(1-e^2)));
    dydotp_dM = -n * yp * (a/r)^3;
    
    % Partial derivative of zdotp w.r.t. a,e,M
    dzdotp_da = 0;
    dzdotp_de = 0;
    dzdotp_dM = 0;
    
    % Partial of PQW2ECI w.r.t. Omega (RAAN)
    dPQW2ECI_dOmega = [-PQW2ECI(2,1), -PQW2ECI(2,2), -PQW2ECI(2,3) ;
                        PQW2ECI(1,1),  PQW2ECI(1,2),  PQW2ECI(1,3) ;
                             0      ,       0      ,        0     ];
                       
    % Partial of PQW2ECI w.r.t. w (argument of perigee)
    dPQW2ECI_dw     = [ PQW2ECI(1,2), -PQW2ECI(1,1),        0      ;
                        PQW2ECI(2,2), -PQW2ECI(2,1),        0      ;
                        PQW2ECI(3,2), -PQW2ECI(3,1),        0     ];
                    
    % Partial of PQW2ECI w.r.t. i (inclination)
    dPQW2ECI_di     = [ sin(Omega)*sin(w)*sin(i),  sin(Omega)*cos(w)*sin(i),  sin(Omega)*cos(i) ;
                       -cos(Omega)*sin(w)*sin(i), -cos(Omega)*cos(w)*sin(i), -cos(Omega)*cos(i) ;
                             sin(w)*cos(i)      ,        cos(w)*cos(i)     ,       -sin(i)     ];
                         
    % Jacobian matrix between Cartesian ECI and Equinoctial elements
    J_CART_KEP = zeros(6,6);
    J_CART_KEP(1:3,1) = PQW2ECI         * [dxp_da; dyp_da; dzp_da];
    J_CART_KEP(4:6,1) = PQW2ECI         * [dxdotp_da; dydotp_da; dzdotp_da];
    J_CART_KEP(1:3,2) = PQW2ECI         * [dxp_de; dyp_de; dzp_de];
    J_CART_KEP(4:6,2) = PQW2ECI         * [dxdotp_de; dydotp_de; dzdotp_de];
    J_CART_KEP(1:3,3) = dPQW2ECI_di     * [xp; yp; zp];
    J_CART_KEP(4:6,3) = dPQW2ECI_di     * [xdotp; ydotp; zdotp];
    J_CART_KEP(1:3,4) = dPQW2ECI_dOmega * [xp; yp; zp];
    J_CART_KEP(4:6,4) = dPQW2ECI_dOmega * [xdotp; ydotp; zdotp];
    J_CART_KEP(1:3,5) = dPQW2ECI_dw     * [xp; yp; zp];
    J_CART_KEP(4:6,5) = dPQW2ECI_dw     * [xdotp; ydotp; zdotp];
    J_CART_KEP(1:3,6) = PQW2ECI         * [dxp_dM; dyp_dM; dzp_dM];
    J_CART_KEP(4:6,6) = PQW2ECI         * [dxdotp_dM; dydotp_dM; dzdotp_dM];
    
    % Convert VCM Equinoctial covariance to Cartesian ECI
    J = inv(J_KEP_EQ) * inv(J_CART_KEP);

end