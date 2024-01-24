function [xmin,ymin,converged,nbisect,nshift,xbuf,ybuf] = shift_bisect_minsearch( ...
    fun,x0,y0,Nbisectmax,Nshiftmax,tolX,tolY,verbose)

% Find the minimum of the function "y = fun(x)" given the initial search
% interval (x1,x2).
%
% Because this function uses an algorithm that performs both shifting and
% bisection in search for the minimum, (x1,x2) does not necessarily have to
% bound the minimum.
%
% This function finds one and only one minimum, so if there are multiple
% minima over (x1,x2), then this function will only find one at most.
%

% Set defaults for inputs
Nargin = nargin;

% Initial number of grid points (should be 5 or larger)
N0 = numel(x0);
if (N0 < 5)
    error('shift_bisect_minsearch requires five or more points');
end
if N0 ~= numel(y0)
    error('shift_bisect_minsearch needs same number of (x,y) grid points');
end
if any(diff(x0) < 0)
    error('shift_bisect_minsearch needs monotonically increasing grid points');
end

% Max number of bisections
na = 4;
if (Nargin < na) || isempty(Nbisectmax)
    Nbisectmax = 100;
end

% Max number of shifting iterations
na = na+1;
if (Nargin < na) || isempty(Nshiftmax)
    Nshiftmax = 100;
end

% Tolerance for x convergence (can be NaN to only use y convergence)
na = na+1;
if (Nargin < na) || isempty(tolX)
    tolX = 1e-4;
end

% Tolerance for y convergence (must be positive)
na = na+1;
if (Nargin < na) || isempty(tolY)
    tolY = 1e-4;
end
if isnan(tolY) || tolY < 0
    error('Invalid tolY parameter');
end

% Verbosity (0 => no reports, 1 => few reports to user, 2 => many reports)
na = na+1;
if (Nargin < na) || isempty(verbose)
    verbose = 1;
end

% Check to see if output buffering is required
buffering = (nargout > 5);

% Find minimum value over initial grid
if max(y0) ~= min(y0)
    % Find min point
    [~, nmin] = min(y0);
else
    % Use center point if all initial y values are equal
    nmin = round((Ninitial-1)/2+1);
end

% Initialize the 5-point bisection search buffer

if (nmin < 3)
    % Minimum is near or at left boundary
    n1 = 1;
    n2 = 5;
elseif (nmin > N0-2)
    % Minimum is near or at right boundary
    n1 = N0-4;
    n2 = N0;
else
    % Minimum point in center somewhere
    n1 = nmin-2;
    n2 = nmin+2;
end

x = x0(n1:n2);      % Holds 5 x values
y = y0(n1:n2);      % Holds 5 y values
c = false(size(x));  % Holds 5 logical values (true => y needs to be calculated)

% Initialize the buffers to hold new variables

xnew = zeros(size(x));
ynew = xnew;
cnew = true(size(x));

% Perform iterative search

if (verbose > 0)
    if buffering
        bufstr = ' (with output buffering)';
    else
        bufstr = '';
    end
    disp(['Performing shifting + bisection search for minimum' bufstr]);
end

% Iteration and convergence flags
still_iterating  = true; converged = false;

nshift = 0;  % Counts shifts to left or right
nbisect = 0; % Counts bisections

if buffering
    xbuf = x0;
    ybuf = y0;
end

% Iterate until converged or stopped
while still_iterating

    % Calculate any required y values
    for n=1:5
        if c(n)
            y(n) = fun(x(n));
            c(n) = false;
            if buffering
                xbuf = cat(2,xbuf,x(n));
                ybuf = cat(2,ybuf,y(n));
            end
        end
    end
    
    % Find min among the 5 calculated values
    [ymin, nmin] = min(y);
    if (ymin == max(y))
        % Use center point if all y values are equal
        nmin = 3;
    end
    xmin = x(nmin);
    
    if (verbose > 1)
        disp(['  (x,y) = ' num2str([x(nmin) y(nmin)])]);
    end    
    
    % Either shift if min is on edge, or bisect if min is interior
    
    if (nmin == 1)

        % First x point is the lowest y value (perform shift to left)
        
        if (verbose > 1)
            disp([' Shifting left, nshift = ' num2str(nshift)]);
        end    

        xnew(2:5) = x(1:4);
        ynew(2:5) = y(1:4);
        cnew(2:5) = false;

        dx        = x(2)-x(1);
        xnew(1)   = x(1)-dx;
        cnew(1)   = true;
        
        nshift    = nshift+1;

    elseif (nmin == 5)

        % Last x point is the highest y value (perform shift to right)

        if (verbose > 1)
            disp([' Shifting right, nshift = ' num2str(nshift)]);
        end    
        
        xnew(1:4) = x(2:5);
        ynew(1:4) = y(2:5);
        cnew(1:4) = false;

        dx        = x(5)-x(4);
        xnew(5)   = x(5)+dx;
        cnew(5)   = true;
        
        nshift    = nshift+1;

    else
        
        % One of the interior points has the highest y value, so check for
        % convergence, and set up to perform next bisection if required

        nlo = nmin-1; nhi = nmin+1;
        
        % Check for convergence
        
        dxlo = abs(x(nmin)-x(nlo));
        dxhi = abs(x(nmin)-x(nhi));
        dx   = max([dxlo dxhi]);
        
        if verbose > 2
            disp([' dx = ' num2str(dx) ' dx/tolX = ' num2str(dx/tolX)]);
        end
        
        if (dx < tolX) || (dx < 10*eps(x(nmin)))
            
            converged = true;
            if (verbose > 0)
                disp(['Convergence achieved because dx < tolX or dx < epsX (nbisect = ' ...
                    num2str(nbisect) ' nshift = ' num2str(nshift) ')']);
            end 
            
        else
            
            dylo = y(nlo)-y(nmin);
            dyhi = y(nhi)-y(nmin);
            dy   = max([dylo dyhi]);
            
            if verbose > 2 && ~isnan(tolY)
                disp([' dy = ' num2str(dy) ' dy/tolY = ' num2str(dy/tolY)]);
            end
            
            if (dy < tolY)
                converged = true;
                if (verbose > 0)
                    disp(['Convergence achieved because dy < tolY (nbisect = ' ...
                        num2str(nbisect) ' nshift = ' num2str(nshift) ')']);
                end                
            end
            
        end

        % If not yet converged then set up to perform next bisection
        
        if ~converged

            xnew(1) = x(nlo);
            ynew(1) = y(nlo);
            cnew(1) = false;

            xnew(2) = (x(nlo)+x(nmin))/2;
            cnew(2) = true;

            xnew(3) = x(nmin);
            ynew(3) = y(nmin);
            cnew(3) = false;

            xnew(4) = (x(nhi)+x(nmin))/2;
            cnew(4) = true;

            xnew(5) = x(nhi);
            ynew(5) = y(nhi);
            cnew(5) = false;

            nbisect = nbisect+1;
            
        end
        
    end

    % Check for NaN values, or too many shift or bisection iterations
    
    if any(isnan(y))
        
        % NaN values detected
        if (verbose > 0)
            disp('NaN y-values detected; returning unconverged result');
        end
        still_iterating = false;
  
    elseif (nbisect > Nbisectmax)

        % Too many bisections
        if (verbose > 0)
            disp('Maximum number of bisections exceeded; returning unconverged result');
        end
        still_iterating = false;

    elseif (nshift > Nshiftmax)

        % Too many shifts
        if (verbose > 0)
            disp('Maximum number of shifts exceeded; returning unconverged result');
        end
        still_iterating = false;

    end
    
    % Check for convergence        
    if converged
        % Stop iterating if converged
        still_iterating = false;
    else
        % Copy the new buffers to set up for next iteration
        x = xnew;
        y = ynew;
        c = cnew;
    end
    
end

% Sort the buffered values
if buffering
    [xbuf,srt] = sort(xbuf);
    ybuf = ybuf(srt);
end

return;
end