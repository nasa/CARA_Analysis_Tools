function Nxint = NvNiNbIntegration(zenmax,phibeg,phiend,                ...
                    a,si,clO,slO,cbO,sbO,cbS,sbS,                       ...
                    Zint,Fint,UniformDist,                              ...
                    I2AbsTol,I2RelTol,I2MaxFunEvals)

% NvNiNbIntegration - Integration function for Nv, Ni, or Nb.
% Syntax: Nxint = NvNiNbIntegration(zenmax,phibeg,phiend,               ...
%                   a,si,clO,slO,cbO,sbO,cbS,sbS,                       ...
%                   Zint,Fint,UniformDist,                              ...
%                   I2AbsTol,I2RelTol,I2MaxFunEvals)
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Turn warning off for the integration process
warning_off = true; warning('off');

% Handle quad2d only mode
if I2MaxFunEvals < 0
    
    % This special mode forces the use of quad2d
    Nxint = Nv;
    
else

    % Integrand function
    if isempty(Zint)
        % Integrand for Nv or Ni integrals, with empty Fb interpolation grids.
        Nxfun = @(ttt) ...
            NvNiNbThetaIntegrand(ttt,a,si,clO,slO,cbO,sbO,cbS,sbS,UniformDist);
    else
        % Integrand for Nb integrals, with populated Fb interpolation grids.
        % Multiply by brighter-than-recommended fraction of satellites, Fb, 
        % interpolated from the grid of input values.
        Nxfun = @(ttt) ...
            interp1(Zint,Fint,ttt) .* ...
            NvNiNbThetaIntegrand(ttt,a,si,clO,slO,cbO,sbO,cbS,sbS,UniformDist);    
    end

    % Calculate integral using specified tolerances
    Nxint = integral(Nxfun,0,zenmax,'AbsTol',I2AbsTol,'RelTol',I2RelTol);
    
end

% Use quad2d as backup integrator
if isnan(Nxint) || isinf(Nxint)
    
    % Define the backup integrand function
    si2 = si^2;
    Nxfun = @(ttt,ppp)NvNiNbThetaPhiIntegrand(ttt,ppp,                  ...
                        a,si2,clO,slO,cbO,sbO,cbS,sbS,                  ...
                        Zint,Fint,UniformDist);

    % Calculate the backup integral using specified parameters
    Nxint = quad2d(Nxfun,0,zenmax,phibeg,phiend,                        ...
                    'AbsTol',I2AbsTol,                                  ...
                    'RelTol',I2RelTol,                                  ...
                    'MaxFunEvals',max(2000,abs(I2MaxFunEvals)));
                
    if isnan(Nxint) || isinf(Nxint)
        error('Integration failed.');
    end
    
end

if warning_off; warning('on'); end

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