function y_lpc = makelpc(x,f_s)
% makelpc(x,f_s) makes a 20ms AR(12) analysis of x
% and returns av matrix of lpc coefficients.

window_size = 20e-3; % 20ms
inc = f_s * window_size; % number of frames

t = [inc inc 0];
y_lpc = lpcauto(x,12,t);

end
