function [Nx,out] = NumExpectedKessler(SDADeg,ObsLat,SunLat,CalcNv,CalcNi,CalcNb,ZaInterp,FbInterp,KnownLT,params)
% NumExpectedKessler - Calculate satellites above observer
% Syntax: [Nx,out] = NumExpectedKessler(SDADeg,ObsLat,SunLat,CalcNv,CalcNi,CalcNb,ZaInterp,FbInterp,KnownLT,params);
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
% Calculate the maximum statistically expected number of satellites above
% ground-based observers for a given solar depression angle (SDA)
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Initialize output
out = [];

% Check for null calculation
if ~(CalcNv || CalcNi || CalcNb)
    warning('No (Nv,Ni,Nb) calculations being performed');
end

% Constants
Re = 6378.137;
EarthObliquityDeg = 23.5;
DegPerRad = 180/pi;
twopi = 2*pi;
cbOeps = 30*eps;

% 2D integration accuracy tolerance parameters
IntMode = params.Evaluation.NumExpectedKesslerIntMode;
if numel(IntMode) == 3
    % Integration parameters specified explicitly
    I2AbsTol0     = IntMode(1); 
    I2RelTol      = IntMode(2);
    I2MaxFunEvals = IntMode(3);
elseif numel(IntMode) == 1
    AbsIntMode = abs(IntMode);
    if AbsIntMode == 1
        % Normal accuracy mode
        I2AbsTol0 = 1e-6; 
        I2RelTol = 5e-5;
        I2MaxFunEvals = 10000;
    elseif AbsIntMode == 2
        % High accuracy mode
        I2AbsTol0 = 1e-8;
        I2RelTol = 1e-6;
        I2MaxFunEvals = 20000;
    else
        error('Invalid integration mode');
    end
else
    error('Invalid integration mode');
end

% Special mode to force the use of quad2d integrations
if IntMode < 0
    I2MaxFunEvals = -I2MaxFunEvals;
end

% Uniform vs nonuniform (i.e., Kessler) distribution used for each
% orbital shell of constellation satellites
UniformDist = params.Evaluation.UniformDist;

% Observer latitudes to evaluate (degrees)
ObsLat = unique(ObsLat); NObsLat = numel(ObsLat);
if any(abs(ObsLat) > 90)
    error('Observer location latitude(s) out of range');
end

% Sub-solar latitudes to evaluate (degrees)
SunLat = unique(SunLat); NSunLat = numel(SunLat);
if any(abs(SunLat) > EarthObliquityDeg)
    error('Sub-solar point latitude(s) out of range');
end

% Observer and sub-solar latitudes and trig functions
bObs = ObsLat/DegPerRad; cbObs = cos(bObs); sbObs = sin(bObs);
bSun = SunLat/DegPerRad; cbSun = cos(bSun); sbSun = sin(bSun);

% Maximum zenith angle considered for evaluation
zenmax = params.Evaluation.MaxZenith/DegPerRad;
czmax = cos(zenmax); szmax2 = 1-czmax^2; szmax = sqrt(szmax2);

% Cosine and sine of zenith angle (i.e., theta) of observer-to-sun vector, 
% ctS = cos(thetaSun) = cos(SDA+90deg)
if SDADeg < 0 || SDADeg > 90
    error('Invalid solar depression angle');
end
ctS = cos((SDADeg+90)/DegPerRad); stS = sqrt(1-ctS^2);    

% Number of shells and constellation satellites per shell
Nshell  = numel(params.New.Nc);

% Total number of latitude grid points
NumLatGrid = NObsLat*NSunLat;
LatGridDim = [NObsLat,NSunLat];

% Number of shells and constellation satellites per shell Nlat x Nshell
Ncshell = reshape(params.New.Nc,[1 Nshell]);
NcRep = repmat(Ncshell,[NumLatGrid 1]);

% Calculate the semi-major axis values for all orbital shells of the new
% constellation
SMA = reshape(Re+params.New.Altitude_km,[1 Nshell]);
SMARep = repmat(SMA,[NumLatGrid 1]);

% Thin shell approximation for Nv at maximum zenith
bthin = sqrt( SMA.^2 - Re^2 * szmax2 );
fthin = Re*szmax2 + czmax*bthin;
Nvthin0 = (1 - fthin./SMA) / 2;
% Nvthin = Ncshell .* Nvthin0;

% Prepare for Ni or Nb calculation by calculating sunlit azimuth angle
% factor, phiithin, as given by eq (12) of Hall (2021 AMOS).
% No sunlit satellites exist within constellation shells
% that have phiithin = 0, so this save CPU time by preventing the
% calculation of Ni and Nb in those cases
if CalcNi || CalcNb
    CalcNi_or_CalcNb = true;
    SDARad = SDADeg/DegPerRad;
    cSDA = cos(SDARad); sSDA = sin(SDARad);
    rho = bthin - Re*czmax;
    gtop = sSDA*(rho*czmax+Re)-sqrt(rho.*(2*Re*czmax+rho));
    gbot = rho*szmax*cSDA;
    gam = gtop ./ gbot;
    gam(gam < -1) = -1;
    gam(gam >  1) =  1;
    phiithin = acos(gam);
else
    CalcNi_or_CalcNb = false;
    phiithin = [];
end

% Absolute tolerance depends on thin shell approximation
I2AbsTol = max(1e-10,I2AbsTol0*Nvthin0);
I2AbsRep = repmat(I2AbsTol,[NumLatGrid 1]);

if ~UniformDist
    % Inclinations required for nonuniform Kessler distribution
    incl = params.New.Inclination_deg;
    % indx = incl > 90;
    % incl(indx) = 180-incl(indx);
    incl = reshape(incl/DegPerRad,[1 Nshell]);
    siRep = repmat(sin(incl),[NumLatGrid 1]);    
else
    siRep = NaN(NumLatGrid,Nshell);
end

% Calculate arrays of subscripts (NEED TO VECTORIZE!!!!)
nObsArr  = NaN(1,NumLatGrid); nSunArr  = NaN(1,NumLatGrid);
bObsArr  = NaN(1,NumLatGrid); bSunArr  = NaN(1,NumLatGrid);
cbObsArr = NaN(1,NumLatGrid); sbObsArr = NaN(1,NumLatGrid); 
cbSunArr = NaN(1,NumLatGrid); sbSunArr = NaN(1,NumLatGrid); 
for i=1:NumLatGrid
    % Subscripts of Obs and Sun lattitude arrays
    [nO,nS] = ind2sub(LatGridDim,i);
    nObsArr(i)  = nO;         nSunArr(i) = nS;
    bObsArr(i)  = ObsLat(nO); bSunArr(i) = SunLat(nS);
    cbObsArr(i) = cbObs(nO); sbObsArr(i) = sbObs(nO);
    cbSunArr(i) = cbSun(nS); sbSunArr(i) = sbSun(nS);
end

% Initialize arrays of local time (LT) and number of all, sunlit, and
% brighter-than-recommended satellites (Nv,Ni,Nb)
LT = NaN(1,NumLatGrid); 
% Nv = LT;
% Ni = LT;
% Nb = LT;
NvShell = NaN(NumLatGrid,Nshell);
NiShell = NaN(NumLatGrid,Nshell);
NbShell = NaN(NumLatGrid,Nshell);

% Allocate calculation buffers
NvBuf = NaN(NumLatGrid,1);
NiBuf = NaN(NumLatGrid,1);
NbBuf = NaN(NumLatGrid,1);

% Loop over the orbital shells
for nshell = 1:Nshell
    
    % Illuminated azimuth half-width
    if CalcNi_or_CalcNb
        phiithin_shell = phiithin(nshell);
    else
        phiithin_shell = [];
    end

    % Calculate arrays of values in parallel, just in case NumLatGrid is large
    parfor i=1:NumLatGrid

        % Cosine and sine values of observer and sub-solar latitudes
        cbO = cbObsArr(i); sbO = sbObsArr(i); bODeg = bObsArr(i);
        cbS = cbSunArr(i); sbS = sbSunArr(i); bSDeg = bSunArr(i);

        % Calculate the cosine of the local time of the observer that
        % corresponds to the input SDA (LT = lambdaObs = lO)
        if abs(cbO) > cbOeps
            % Cos(LT) for non-polar points
            clO = (sbO*sbS-ctS)/(cbO*cbS);
            % If local time angle is known to exist, the account for possible
            % round-off errors
            if KnownLT
                clO = max(-1,min(1,clO));
            end
        else
            % Effective Cos(LT) for polar points 
            if KnownLT
                clO = 0; % Set pole LT to 6am = Pi/2
            else
                if bODeg > 0
                    % N pole
                    SDAPole = -bSDeg;
                else
                    % S pole
                    SDAPole = bSDeg;
                end
                if SDAPole == SDADeg
                    clO = 0; % Set pole LT to 6am = Pi/2
                else
                    clO = Inf; % Set pole LT non-real value
                end
            end
        end
        
        % If |cos(LT)| <= 1 then LT angle is real
        if -1 <= clO && clO <= 1

            % Sine of LT angle and associated sine value
            LT(i) = acos(clO);
            slO = sqrt(1-clO^2); % Use dawn-side value which produces same output as dusk

            % Loop over the orbital shells
            % for nshell = 1:Nshell

                % Number of constellation satellites in this shell
                Nc = NcRep(i,nshell);

                % Semi-major axis of this shell
                a = SMARep(i,nshell);

                % Sine(inclination)
                si = siRep(i,nshell);

                % Calculate Nv integral, if required
                if CalcNv

                    % Azimuth limits for Nv integration
                    phibeg = 0; phiend = twopi;
                    % No solar illumination parameters required for integration
                    cbSint = []; sbSint = []; 
                    % No brighter-than-recommended fraction grid required
                    Zint = []; Fint = [];

                    % Calculate Nv as a 2D integral over (theta,phi)
                    Nvint = NvNiNbIntegration(zenmax,phibeg,phiend,     ...
                        a,si,clO,slO,cbO,sbO,cbSint,sbSint,             ...
                        Zint,Fint,UniformDist,                          ...
                        I2AbsRep(i),I2RelTol,I2MaxFunEvals);

                    % Issue error if both integration methods fail
                    if isnan(Nvint)
                        error('Failed Nv solid angle integration');
                    else
                        NvBuf(i) = Nc * Nvint;
                    end

                end

                % Set up to calculate Ni or Nb
                if CalcNi_or_CalcNb && phiithin_shell > 0
                    % Integrate over only sunlit azimuth
                    % angles, i.e., over the azimuthal angle limits
                    % bounding the region where the intersection points of
                    % of the LOS and orbital shell are sunlit
                    shat = [-cbS; 0; sbS];
                    nhat0 = [-clO*sbO; -slO*sbO; cbO];
                    ehat0 = [-slO; clO; 0];
                    pS = atan2(ehat0'*shat,nhat0'*shat);
                    pbeg = @(tt)pS-IllumPhiBounds(tt,a,ctS,stS,Re);
                    pend = @(tt)pS+IllumPhiBounds(tt,a,ctS,stS,Re);
                else
                    pbeg = 0; pend = twopi;
                end

                % Calculate Ni integral, if required
                if CalcNi

                    % If no azimuthal angles are sunlit, then Ni also zero
                    if phiithin_shell == 0
                        Niint = 0;
                    else

                        % Calculate the Ni integral

                        % No brighter-than-recommended fraction grid required
                        Zint = []; Fint = [];

                        % Solar illumination required for Ni integration
                        cbSval = cbS; sbSval = sbS;

                        % Calculate Nb as a 2D integral over (theta,phi)                    
                        Niint = NvNiNbIntegration(zenmax,pbeg,pend,     ...
                            a,si,clO,slO,cbO,sbO,cbSval,sbSval,         ...
                            Zint,Fint,UniformDist,                      ...
                            I2AbsRep(i),I2RelTol,I2MaxFunEvals);

                    end

                    % Issue error if both integration methods fail
                    if isnan(Niint)
                        error('Failed Ni solid angle integration');
                    else
                        NiBuf(i) = Nc * Niint;
                    end

                end

                % Calculate Nb, if required
                if CalcNb

                    % If Nithin is zero, then Nb is also zero
                    if phiithin_shell == 0
                        Nbint = 0;
                    else

                        % Brighter-than-recommended fraction grid required
                        Zint = ZaInterp; Fint = FbInterp(nshell,:);

                        % Solar illumination required for Ni integration
                        cbSval = cbS; sbSval = sbS;

                        Nbint = NvNiNbIntegration(zenmax,pbeg,pend,     ...
                            a,si,clO,slO,cbO,sbO,cbSval,sbSval,         ...
                            Zint,Fint,UniformDist,                      ...
                            I2AbsRep(i),I2RelTol,I2MaxFunEvals);

                    end

                    % Issue error if both integration methods fail
                    if isnan(Nbint)
                        error('Failed Nb solid angle integration');
                    else
                        NbBuf(i) = Nc * Nbint;
                    end

                end

            % end
            
        end

    end

    if CalcNv; NvShell(:,nshell) = NvBuf; end
    if CalcNi; NiShell(:,nshell) = NiBuf; end
    if CalcNb; NbShell(:,nshell) = NbBuf; end

end

% Reshape output arrays. Recall: LatGridDim = [NObsLat,NSunLat]

out.ObsLat = ObsLat; out.SunLat = SunLat;
out.LocalTime = reshape(LT,LatGridDim)*DegPerRad;

if CalcNv
    out.NvShell = NaN(NObsLat,NSunLat,Nshell);
    for nshell=1:Nshell
        out.NvShell(:,:,nshell) = reshape(NvShell(:,nshell),LatGridDim);
    end
    out.Nv = sum(out.NvShell,3);
    Nx = out.Nv;
end
if CalcNi
    out.NiShell = NaN(NObsLat,NSunLat,Nshell);
    for nshell=1:Nshell
        out.NiShell(:,:,nshell) = reshape(NiShell(:,nshell),LatGridDim);
    end
    out.Ni = sum(out.NiShell,3);
    Nx = out.Ni;
end
if CalcNb
    out.NbShell = NaN(NObsLat,NSunLat,Nshell);
    for nshell=1:Nshell
        out.NbShell(:,:,nshell) = reshape(NbShell(:,nshell),LatGridDim);
    end
    out.Nb = sum(out.NbShell,3);    
    Nx = out.Nb;
end

% Only output numbers for Nx
Nx(isnan(Nx)) = 0;

% Reshape as required
if NumLatGrid > 1
    if numel(SunLat) == 1
        Nx = reshape(Nx,size(ObsLat));
    elseif numel(ObsLat) == 1
        Nx = reshape(Nx,size(SunLat));
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