function [xmn, xmx] = plot_range(xvec_input, padfactor, zerofactor)

% Calculate sensible plot ranges for an input data vector

% Exclude NaN and +/-Inf data points

ndx = ~isnan(xvec_input) & ~isinf(xvec_input);

if any(ndx)
    
    % Default zero and pad factors

    if (nargin < 3); zerofactor = 0.0; end
    if (nargin < 2); padfactor  = 0.1; end
    
    % Epsilon

    xeps = eps;
    
    % If single padfactor supplied, double it to apply to both side

    if (numel(padfactor) == 1)
        xpad = [padfactor padfactor];
    else
        xpad = padfactor;
    end
    
    % Ensure positive pad factors

    xpad(1) = max(xeps,xpad(1));
    xpad(2) = max(xeps,xpad(2));
    
    % Good data and extrema
    
    xvec = xvec_input(ndx);
    xmn0 = min(xvec);
    xmx0 = max(xvec);

    if (abs(xmx0-xmn0) < xeps)
        if (xmn0 == 0)
            xmn = -1;
            xmx = 1;
        else
            xmn1 = (1-xpad(1))*xmn0;
            xmx1 = (1+xpad(2))*xmx0;
            xmn = min([xmn1 xmx1]);
            xmx = max([xmn1 xmx1]);
        end
    else
        xdel = xmx0-xmn0;
        xmn  = xmn0-xdel*xpad(1);
        xmx  = xmx0+xdel*xpad(2);
    end
    
    % Make lower limit zero if required

    if (zerofactor > xeps)
        if ((xmn > 0) && (xmn <= xmx*zerofactor)); xmn=0; end
    end

else
    
    % All NaN or Inf data
    
    xmn = 0;
    xmx = 1;
    
end

% Output one or two arguments

if nargout == 1
    xmn = [xmn,xmx];
end

return