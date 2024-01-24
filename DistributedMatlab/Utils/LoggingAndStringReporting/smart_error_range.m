function [xmdstr,xlostr,xhistr] = smart_error_range(xmd,xlo,xhi)

% Format strings smartly for an estimated quantity (xmd) and a bracketing
% range (xlo < xmd < xhi)

xdel   = max(abs([xmd-xlo xmd-xhi]));

xmdstr = smart_error_format(xmd,xdel);
xlostr = smart_error_format(xlo,xdel);
xhistr = smart_error_format(xhi,xdel);

xmdstr = smart_exp_format(str2double(xmdstr));
xlostr = smart_exp_format(str2double(xlostr));
xhistr = smart_exp_format(str2double(xhistr));

return;
end