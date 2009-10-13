%% TRAINING
% Terje Gundersen 15.09.2009
close all;
clear all;

%% Read file
[x,f_s]=wavread('data/t01s000228.wav');
[y,fs_y]=wavread('data/t03s000228.wav');


%% Split file in 20ms segments
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
[n,frame_length] = size(F_x);


%% Make AR(p) coefficients
p = 12; % LPC order
X_lpc = lpc(F_x',p);
Y_lpc = lpc(F_y',p);


%% Transformation LPC --> RC
X_lsf = zeros(n,p);
for i=1:n
    X_lsf(i,:) = poly2lsf(X_lpc(i,:));
end
Y_lsf = zeros(n,p);
for i=1:n
    Y_lsf(i,:) = poly2lsf(Y_lpc(i,:));
end


%% EM algorithm
m = 8; % Number of mixture models
gm_obj = gmdistribution.fit(X_lsf,m,'CovType','diagonal'); % EM-alg


%% Matrices
P = posterior(gm_obj,X_lsf); % Posterior probability

% Convert the vector Sigma into a diagonal matrix and invert it.
sigma_diag = zeros(m,p);
for i=1:m
	sigma_diag(i,:) = 1./gm_obj.Sigma(1,:,i);
end

% Compute V and Gamma for each p
V = zeros(m,p);
Gamma = zeros(m,p);
% for i=1:n
for k=1:p
	[V(:,k),Gamma(:,k)] = param(k,m,P,X_lsf,Y_lsf,gm_obj,sigma_diag); 
end
% end

save 'variables';
