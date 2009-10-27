%% TRAINING
% Terje Gundersen 15.09.2009
close all;
clear all;

%% Read file
[x,f_s]=wavread('data/t01s000228.wav');
[y,fs_y]=wavread('data/t03s000228.wav');


%% Split file in 20ms segments and convert to LSF vectors
p = 16; % LPC order
[X_lsf,Y_lsf,n,frame_length] = makelsf(x,y,f_s,p);


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

for k=1:p
	[V(:,k),Gamma(:,k)] = param(k,m,P,X_lsf,Y_lsf,gm_obj,sigma_diag); 
end

save 'variables';
