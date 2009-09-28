%% SCRIPT FILE
% Terje Gundersen 15.09.2009
% Last updated 22.09.2009
close all;
clear all;

%% Read file
[x,f_s]=wavread('data/000358_JF.wav');


%% Split file in 20ms segments
window_size = 20e-3; % 20ms
inc = f_s * window_size; % samples per frame
N_loop = floor(length(x)/inc);

F = [];
for i=1:N_loop
   F(i,:) = x(1+inc*(i-1):inc*i); 
end

if N_loop < length(x)/inc
   test(N_loop+1,:) = x(inc*N_loop:end); 
end

% Dimensions
[m,n] = size(F);


%% Make AR(12) coefficients
lpc_order = 12;
X_lpc = lpc(F',lpc_order);


% Subtract estimated signal from real signal
error = zeros(m,n);
for i=1:m
    error(i,:) = filter([1 X_lpc(i,2:end)], 1, F(i,:));
end


%% Transformation LPC --> RC
X_rc = zeros(m,lpc_order);
for i=1:m
    X_rc(i,:) = poly2rc(X_lpc(i,:));
end
% status = lpcconv('ar','rf',X_lpc,Y_lpc,12);


%% EM algorithm
N_class = 5;
P = zeros(m,N_class);

gm_obj = gmdistribution.fit(X_rc,N_class,'CovType','diagonal'); % EM-alg
P = posterior(gm_obj,X_rc); % Posterior probability



%% Conversion Function
