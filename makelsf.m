function [X_lsf,Y_lsf,n,frame_length,X_lpc,F_x] = makelsf(x,y,f_s,p)
% makelsf(x,y,f_s,p) makes a 20ms AR(12) analysis of x
% and returns av matrix of LSF coefficients.


window_size = 20e-3; % 20ms
inc = f_s * window_size; % samples per frame


% Source vector
N_loop = floor(length(x)/inc);
F_x = zeros(N_loop,inc);
for i=1:N_loop
   F_x(i,:) = x(1+inc*(i-1):inc*i); 
end

if N_loop < length(x)/inc
   F_x(N_loop+1,:) = zeros(inc,1);
   F_x(N_loop+1,1:(length(x)-inc*N_loop)) = x(inc*N_loop+1:end); 
end


% Target vector
N_loop = floor(length(y)/inc);
F_y = zeros(N_loop,inc);
for i=1:N_loop
   F_y(i,:) = y(1+inc*(i-1):inc*i); 
end

if N_loop < length(y)/inc
   F_y(N_loop+1,:) = zeros(inc,1);
   F_y(N_loop+1,1:(length(y)-inc*N_loop)) = y(inc*N_loop+1:end); 
end


% Dimensions
[n,frame_length] = size(F_x); % SIZE(F_x) NEQ SIZE(F_y)

% Make AR(p) coefficients
X_lpc = lpc(F_x',p);
Y_lpc = lpc(F_y',p);

% Transformation LPC --> LSF
X_lsf = zeros(n,p);
for i=1:n
    X_lsf(i,:) = poly2lsf(X_lpc(i,:));
end

Y_lsf = zeros(n,p);
for i=1:n
    Y_lsf(i,:) = poly2lsf(Y_lpc(i,:));
end

end
