function [xmnma,ymnma,xmxma,ymxma,converged,nbisect,x,y,imnma,imxma] = ...
    refine_bounded_extrema(fun,x1,x2,Ninitial,Nbisectmax,extrema_types, ...
                           TolX,TolY,endpoints,verbose,check_inputs)
%==========================================================================
%
% Find and refine all of the the extrema of the function "y = fun(x)" 
% within the interval (x1,x2) and initially divided up into Ninital points.
%
%==========================================================================
%
% This function has two basic usage modes:
%
% MODE 1: The bounds (x1,x2) are input as scalars, and the function first
% calculates Ninitial (x,y) points spanning that range inclusively, and
% subsequently refines all of the extrema found in that initial set using a
% bisection method.
%
% MODE 2: Upon input, Ninitial = [], x1 holds a precaculated x vector, and
% x2 holds the correponding precalculated y vector, meaning that
% x2 = fun(x1).  In this mode "Ninitial" must be input as an empty set.
% The function refines all of the extrema found in that initial (x,y) set
% using a bisection method.
%
%==========================================================================
%
% Input:
%
%   fun         = Anonymous function used to find and refine extrema.
%
%   (x1,x2)     = MODE 1: Scalars holding the x bounds.
%                 MODE 2: Vectors holding the Ninitial (x,y) values, which
%                         have been sorted by increasing x value.
%
%   Ninitial     = MODE 1: Number of initial points spanning (x1,x2)
%                          inclusively.
%                  MODE 2: Empty set.
%
%   Nbisectmax  = Maximum number of bisections allowed for refinement.
%
%   extrema_types = 1 for only minima, 2 for only maxima, and 3 for both.
%
%   TolX        = Tolerance for x-value convergence. If scalar, then this
%                 is the absolute tolerance level. If 2-vector, then these
%                 are the absolute and relative tolerance levels.
%
%   TolY        = Tolerance for y-value convergence. If scalar, then this
%                 is the absolute tolerance level. If 2-vector, then these
%                 are the absolute and relative tolerance levels.
%
%   endpoints   = Flag indicating that extrema at endpoints should also 
%                 be refined. Note, this function does no shifting, so the
%                 endpoint extrema will only be refined within the original
%                 bounds inclusively.
%
%   verbose     = Verbosity
%
%   check_input = Flag to suppress default and error checking of the above
%                 input parameters.  If check_input = false, then all
%                 inputs are assumed to be populated correctly.
%
%==========================================================================
%
% Output:
%
%   xmnma        = Vector holding the x-values of all found/refined minima.
%
%   ymnma        = Vector holding the y-values of all found/refined minima.
%
%   xmxma        = Vector holding the x-values of all found/refined maxima.
%
%   ymxma        = Vector holding the y-values of all found/refined maxima.
%
%   converged   = Flag indicating convergence to within (TolX,TolY) with
%                 fewer than Nbisectmax bisections.
%
%   nbisect     = Number of bisection iterations.
%
%   (x,y)       = Vectors holding (x,y) values holding all of the points
%                 used in the bisection refinement.
%
%   imnma       = Vector of (x,y) vector indices for the minima, so that
%                 xmnma(n) = x(imnma(n)) and ymnma(n) = y(imnma(n)).
%
%   imxma       = Vector of (x,y) vector indices for the maxima, so that
%                 xmxma(n) = x(imxma(n)) and ymxma(n) = y(imxma(n)).
%
%==========================================================================
%
% Author:
%
%   Doyle T. Hall, Omitron Inc.
%
% =========================================================================

% Initializations

xmnma = []; ymnma = []; imnma = [];
xmxma = []; ymxma = []; imxma = [];

nbisect = 0; converged = false;

% Check run mode and set defaults for inputs

mode_1 = ~isempty(Ninitial);

Nargin = nargin;

narg = 11;
if (Nargin < narg) || isempty(check_inputs); check_inputs = true; end

if check_inputs

    % Verbosity (0 => no reports, 1 => few reports to user, 2 => many reports)

    narg = narg-1;
    if (Nargin < narg) || isempty(verbose); verbose = 0; end
    
    % Default is to not include endpoints

    narg = narg-1;
    if (Nargin < narg) || isempty(endpoints); endpoints = false; end
    
    % Default tolerances, absolute and relative
    
    narg = narg-1;
    if (Nargin < narg) || isempty(TolY); TolY = 1e-6; end

    narg = narg-1;
    if (Nargin < narg) || isempty(TolX); TolX = 1e-6; end
        
    % Default extrema types
    
    narg = narg-1;
    if (Nargin < narg) || isempty(extrema_types); extrema_types = 3; end
    
    if (extrema_types < 1) || (extrema_types > 3)
        error('Illegal value for extrema_types: 1 for only minima, 2 for only maxima, and 3 for both.');
    end
    
    % Default max number of bisections
    
    narg = narg-1;
    if (Nargin < narg) || isempty(Nbisectmax); Nbisectmax = 100; end
    
    % Default number of points 
    
    if mode_1
        
        % Mode where (x1,x2) specify the bounding endpoints, and Ninitial
        % specifies the number of initial inclusively spanning points.
        
        if ~isscalar(Ninitial) || ~isnumeric(Ninitial)
            error('MODE 1: Illegal value for Ninitial');
        elseif (Ninitial < 5)
            warning('MODE 1: Ninitial too small (must be 5 or larger); setting Ninitial = 5.')
            Ninitial = 5;
        end
        
        if ~isscalar(x1) || ~isnumeric(x1) || isinf(x1) || isnan(x1)
            error('MODE 1: Illegal value for x1 bound.');
        end

        if ~isscalar(x2) || ~isnumeric(x2) || isinf(x2) || isnan(x2)
            error('MODE 1: Illegal value for x2 bound.');
        end
        
    else
        
        % Mode where (x1,x2) specify the initial, pre-calculated (x,y)
        % vectors
        
        Nx = numel(x1);
        Ny = numel(x2);
        
        if (Nx ~= Ny)
            error('MODE 2: Input (x,y) vectors have unequal dimensions');
        elseif (Nx < 5)
            warning('MODE 2: Input (x,y) vectors have fewer than 5 points; recommend at least 5 for stable operation.');
        end
        
    end
    
end

% If only one tolerance supplied, use it as the absolute tolerance 
if numel(TolX) == 1; TolX = [TolX,NaN]; end
if numel(TolY) == 1; TolY = [TolY,NaN]; end

% NaN values indicate that tolerances are not to be considered
AbsTolXConsidered = ~isnan(TolX(1));
RelTolXConsidered = ~isnan(TolX(2));
AbsTolYConsidered = ~isnan(TolY(1));
RelTolYConsidered = ~isnan(TolY(2));

NoTolXConsidered = ~AbsTolXConsidered && ~RelTolXConsidered;
NoTolYConsidered = ~AbsTolYConsidered && ~RelTolYConsidered;

if NoTolXConsidered && NoTolYConsidered
    % At least one tolerance, X or Y, must be considered
    error('Tolerance in X and/or Y must be specified');
else
    AllTolXConsidered = AbsTolXConsidered & RelTolXConsidered;
    AllTolYConsidered = AbsTolYConsidered & RelTolYConsidered;
end

% Set up initial (x,y) grid.

if mode_1
    
    % Calculate the initial (x,y) values for mode 1

    if (verbose > 0)
        disp(' ');
        disp(['Calculating ' num2str(Ninitial) ' initial points.']);
    end
    
    % Min and max bounds
    
    minx = min(x1,x2);
    maxx = max(x1,x2);
    
    % Check for finite interval

    if minx == maxx
        error('MODE 1: (x1,x2) bounds cannot be equal.');
    end
    
    % Calculate initial x vector

    x = linspace(minx,maxx,Ninitial);

    % Calculate the y values, assuming vectorized "fun".

    y = fun(x);
    
else
    
    % Use the input initial (x,y) values for mode 2
    
    x = x1;
    y = x2;
    
    if (verbose > 0)
        disp(' ');
        disp(['Using ' num2str(numel(x)) ' initial input points.']);
    end
    
    % Check for sorted input

    if check_inputs
        
        diffx = diff(x);
        
        if any(diffx == 0)
            error('MODE 2: Input x vector cannot have redundant values.');
        elseif min(diffx) <= 0
            warning('MODE 2: Input x vector not sorted - sorting.');
            [x,nsrt] = sort(x);
            y = y(nsrt);
        end
        
    end
    
    % x1 = x(1);
    % x2 = x(end);
    
end

% Start the interative process to find and refine the extrema

still_refining = true;

while still_refining
    
    % Find all interior extrema
    
    switch extrema_types
        case 1 % only minima
            [~,~,ymnma,imnma] = extrema(y,endpoints);
            xmnma = x(imnma);
        case 2 % only maxima
            [ymxma,imxma,~,~] = extrema(y,endpoints);
            xmxma = x(imxma);
        case 3 % both minima and maxima
            [ymxma,imxma,ymnma,imnma] = extrema(y,endpoints);
            xmnma = x(imnma);
            xmxma = x(imxma);
    end
    
    % Combined extrema
    
    % xexma = [xmnma xmxma];
    % yexma = [ymnma ymxma];
    iexma = [imnma imxma];
    Nexma = numel(iexma);

    % If there are any extrema, see if they need refining
    
    if (Nexma == 0)
        
        still_refining = false;
        converged = true;
        
        if verbose
            disp('WARNING: No interior extrema found');
        end
        
    else
        
        % Current number of points

        N = numel(x);
        Nm1 = N-1;
    
        % Define new points bracketing the extrema using bisection
        
        nnew = 0;
        xnew = NaN(1,2*Nexma);
        DelXeps = false;
        
        % if nbisect == Nbisectmax
        %     keyboard;
        % end
        
        for nnmd=1:Nexma
            
            % Mid-point for bisection
            nmd = max(2,min(Nm1,iexma(nnmd)));
            
            % Low and high x-value points
            nlo = nmd-1; DelXlo = x(nmd)-x(nlo);
            nhi = nmd+1; DelXhi = x(nhi)-x(nmd);
            
            % Bisect low-X segement of required
            DelX = DelXlo;
            if (DelX <= max(eps(x(nmd)),eps(x(nlo)))); DelXeps = true; end
            if NoTolXConsidered
                % For no x tolerences, keep the low and high x spacing
                % within a factor of two
                if DelXeps
                    DelXabovetol = false;
                else
                    DelXabovetol = DelXlo > 2*DelXhi;
                end
            else
                if AllTolXConsidered
                    DelXabovetol = (DelX > TolX(1)) && ...
                                   (DelX > TolX(2)*max(abs([x(nmd) x(nlo)])));
                else
                    if AbsTolXConsidered
                        DelXabovetol = (DelX > TolX(1));
                    else
                        DelXabovetol = (DelX > TolX(2)*max(abs([x(nmd) x(nlo)])));
                    end
                end
            end
            if NoTolYConsidered
                DelYabovetol = false;
            else
                DelY = abs(y(nmd)-y(nlo));
                if AllTolYConsidered
                    DelYabovetol = (DelY > TolY(1)) && ...
                                   (DelY > TolY(2)*max(abs([y(nmd) y(nlo)])));
                else
                    if AbsTolYConsidered
                        DelYabovetol = (DelY > TolY(1));
                    else
                        DelYabovetol = (DelY > TolY(2)*max(abs([y(nmd) y(nlo)])));
                    end
                end
            end
            if (DelXabovetol || DelYabovetol)
                nnew = nnew+1;
                xnew(nnew) = (x(nmd)+x(nlo))/2;
            end
            
            % Bisect high-X segement of required
            DelX = DelXhi;
            if (DelX <= max(eps(x(nmd)),eps(x(nhi)))); DelXeps = true; end
            if NoTolXConsidered
                % For no x tolerences, keep the low and high x spacing
                % within a factor of two
                if DelXeps
                    DelXabovetol = false;
                else
                    DelXabovetol = DelXhi > 2*DelXlo;
                end
            else
                if AllTolXConsidered
                    DelXabovetol = (DelX > TolX(1)) && ...
                                   (DelX > TolX(2)*max(abs([x(nhi) x(nmd)])));
                else
                    if AbsTolXConsidered
                        DelXabovetol = (DelX > TolX(1));
                    else
                        DelXabovetol = (DelX > TolX(2)*max(abs([x(nhi) x(nmd)])));
                    end
                end
            end
            if NoTolYConsidered
                DelYabovetol = false;
            else
                DelY = abs(y(nhi)-y(nmd));
                if AllTolYConsidered
                    DelYabovetol = (DelY > TolY(1)) && ...
                                   (DelY > TolY(2)*max(abs([y(nhi) y(nmd)])));
                else
                    if AbsTolYConsidered
                        DelYabovetol = (DelY > TolY(1));
                    else
                        DelYabovetol = (DelY > TolY(2)*max(abs([y(nhi) y(nmd)])));
                    end
                end
            end
            if (DelXabovetol || DelYabovetol)
                nnew = nnew+1;
                xnew(nnew) = (x(nmd)+x(nhi))/2;
            end
            
        end
        
        % Check for convergence, which occurs when there are no new points
        % produced via bisections
    
        if (nnew == 0)
            converged = true;
            still_refining = false;
        elseif DelXeps
            converged = false;
            still_refining = false;
            if (verbose > 0)
                disp(['WARNING: Bisected delta-x comparable to epsilon,' ...
                    ' possibly due to noisy/discontinuous objective function']);
            end
        end
        
        % If still refining, then set up for the next bisection iteration
           
        if still_refining

            if verbose > 0
                disp(['Bisection iteration ' num2str(nbisect+1) ...
                    ' requires ' num2str(nnew) ' new points']);
            end
            
            xnew = unique(xnew(1:nnew));
            ynew = fun(xnew);
            
            xnew = cat(2,x,xnew);
            ynew = cat(2,y,ynew);
            
            [x,nsrt] = sort(xnew);
            y = ynew(nsrt);

            % Check if the number of bisections is too large
            
            nbisect = nbisect+1;
            if (nbisect > Nbisectmax)
                % Discontinue refining
                still_refining = false;
                if (verbose > 0)
                    disp('WARNING: Maximum number of bisections exceeded');
                end
                % Recalculate the extrema, since new points have been added
                switch extrema_types
                    case 1 % only minima
                        [~,~,ymnma,imnma] = extrema(y,endpoints);
                        xmnma = x(imnma);
                    case 2 % only maxima
                        [ymxma,imxma,~,~] = extrema(y,endpoints);
                        xmxma = x(imxma);
                    case 3 % both minima and maxima
                        [ymxma,imxma,ymnma,imnma] = extrema(y,endpoints);
                        xmnma = x(imnma);
                        xmxma = x(imxma);
                end
            end
            
        end
    
    end
    
end

% Report results

if verbose > 0
    if converged
        disp(['Convergence achieved after ' num2str(nbisect) ' iterations']);
    else
        disp(['WARNING: No convergence achieved after ' num2str(nbisect) ' iterations']);
    end
end

return;
end