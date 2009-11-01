function [X_lpc,E,F_x] = split_lpc(x,f_s,p)


window_size = 10e-3; % 10ms
inc = f_s * window_size; % samples per frame

% Source vector
N_x = length(x);
N_loop = floor(N_x/inc);
F_x = zeros(N_loop,inc);
for i=1:N_loop
   F_x(i,:) = x(1+inc*(i-1):inc*i); 
end

if N_loop < N_x/inc
   F_x(N_loop+1,:) = zeros(inc,1);
   F_x(N_loop+1,1:(N_x-inc*N_loop)) = x(inc*N_loop+1:end); 
end

% Make AR(p) coefficients
[X_lpc,E] = lpc(F_x',p);

end
