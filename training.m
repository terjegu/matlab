%% TRAINING
% Terje Gundersen 29.10.2009
close all;
clear all;

%% Read file
[y,fs]=wavread('data/t01s000228.wav');
[x,fs_y]=wavread('data/t03s000228.wav');

[Y,X] = lpcdtw(y,x,fs); % returns time aligned lpc coefficients

%% Convert LPC to LSF
[fn,fl] = size(X);
p = fl-1;
X_lsf = zeros(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X(i,:));
end

Y_lsf = zeros(fn,p);
for i=1:fn
    Y_lsf(i,:) = poly2lsf(Y(i,:));
end

%% Load GMM
load 'gmm';

%% Compute unknown parameters
P = posterior(gm_obj,X_lsf); % Posterior probability

% Convert the vector Sigma into a diagonal matrix and invert it.
sigma_diag = zeros(m,p);
for i=1:m
	sigma_diag(i,:) = 1./gm_obj.Sigma(1,:,i);
end

% Compute V and Gamma for each p
V = zeros(m,p);
Gamma = zeros(m,p);

for k=1:p
	[V(:,k),Gamma(:,k)] = param(k,m,P,X_lsf,Y_lsf,gm_obj,sigma_diag); 
end

save('variables','V','Gamma','P','sigma_diag');
