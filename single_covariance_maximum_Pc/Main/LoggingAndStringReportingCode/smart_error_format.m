function [xstr,estr,pstr] = smart_error_format(x,e,esigfig)

% Create formats for a quantity and its uncertainty, with a reasonable
% number of significant figures for each.  More specifically, create
% output strings that express the quantity, x, and error, e, 
% with enough digits to show 2 significant figures on the error, 
% and the corresponding number of significant figures on the quantity
% itself.

% Special handling for nonpositive uncertainties, etc.

if (e <= 0) || isnan(x) || isnan(e) || isinf(x) || isinf(e)
    xstr = num2str(x);
    estr = num2str(e);
    pstr = [xstr ' \pm ' estr]; 
    return;
end

if nargin < 3
    % Default to 2 significant figures for errors
    esigfig = 2;
end
esigfig = max(1,round(esigfig));

% Calculate the power-of-ten exponent of x

a = abs(x);

if (x ~= 0)
    p = floor(log10(a));
else
    p = floor(log10(e));
end

% Calculate the number of decimal places to the right of the decimal
% point to use in the formatting

fnorm = 10^p;

xdec = max(1,fix(-log10(e/fnorm)+esigfig));
edec = esigfig-1;

% xdec = fix(-log10(e/fnorm)+2);
% edec = 1;

xfmt = ['%-0.' num2str(xdec) 'e'];
xstr0 = strrep(num2str(x,xfmt),'E','e'); % Use lower-case e if exp format

efmt = ['%-0.' num2str(edec) 'e'];
estr0 = strrep(num2str(e,efmt),'E','e'); % Use lower-case e if exp format

% Ensure that there isn't a more compact way to express the two numbers

x0 = str2double(xstr0);
[xstr1,~,~] = smart_exp_format(x0,xdec+1,true);

xstr2 = num2str(x0);

l0 = length(xstr0); l1 = length(xstr1); l2 = length(xstr2);

if l2 <= max(l0,l1)
    xstr = xstr2;
elseif l1 < l0
    xstr = xstr1;
else
    xstr = xstr0;
end

e0 = str2double(estr0);
[estr1,~,~] = smart_exp_format(e0,edec+1,true);

estr2 = num2str(e0);

l0 = length(estr0); l1 = length(estr1); l2 = length(estr2);

if l2 <= max(l0,l1)
    estr = estr2;
elseif l1 < l0
    estr = estr1;
else
    estr = estr0;
end

% Handle the special case where neither xstr or estr use exponential
% notation, and the xstr needs to have zeros added for a sensible match

[xstr] = pad_x_to_match_e(xstr,estr);

% Create a string appropriate for a plot if required

if (nargout > 2)
    
    k = strfind(xstr,'e');
    if isempty(k)
        pstr = [plot_notation(xstr,'') ' \pm ' plot_notation(estr,'')];
    else
        xxstr = xstr(1:k-1); 
        pp = str2double(xstr(k+1:end));
        ee = str2double(estr)/10^pp;
        eefmt = ['%-0.' num2str(edec+1) 'g'];
        eestr = num2str(ee,eefmt);
        kk = strfind(eestr,'e');
        if isempty(kk)
            [xxstr] = pad_x_to_match_e(xxstr,eestr);
            pstr = [xxstr ' \pm ' eestr];
            if (pp ~= 0)
                pstr = plot_notation(['(' pstr ')E' num2str(p)],' ');
            end
        else
            pstr = [plot_notation(xstr,'') ' \pm ' plot_notation(estr,'')];
        end
    end
    
end

return;
end

% =========================================================================

function xout = plot_notation(xstr,pad)

% Convert to TeX style plot notation

if (nargin < 1) || isempty(pad); pad=''; end

k = strfind(lower(xstr),'e');

if isempty(k)
    xout = xstr;
else
    xout = [xstr(1:k-1) pad '\times' pad '10^{' xstr(k+1:end) '}'];
end

return;
end

% =========================================================================

function [xstr] = pad_x_to_match_e(xstr,estr)

% Handle the special case where neither xstr or estr use exponential
% notation, and the xstr needs to have zeros added for a sensible match

e_xstr = strfind(xstr,'e');
e_estr = strfind(estr,'e');

if isempty(e_xstr) && isempty(e_estr)
    
    % No exponential in either xstr or estr
    
    p_estr = strfind(estr,'.');
    
    if ~isempty(p_estr)
        
        % Decimal point in estr
        
        p_xstr = strfind(xstr,'.');
        
        % If no decimal point in xstr, then it must be an integer,
        % so add '.0' to then end
        
        if isempty(p_xstr)
            xstr = [xstr '.0'];
            p_xstr = strfind(xstr,'.');
        end
        
        % The number of numerals in estr to the right of the decimal point
        
        Ne = length(estr)-p_estr;
        Nx = length(xstr)-p_xstr;
        
        % Add zeros as required
        
        if (Ne > Nx)
            xstr = [xstr repmat('0',[1 Ne-Nx])];
        end
        
    end
    
end

return;
end