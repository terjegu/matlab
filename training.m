%% TRAINING
% Terje Gundersen 29.10.2009
close all;
clear all;

%% Load GMM
load 'gmm128';
load 'wavfiles';
N = 641;          % Number of training vectors

%% Compute V and Gamma
p = gm_obj.NDimensions;
m = gm_obj.NComponents;
P = posterior(gm_obj,X_lsf(1:N,:)); % Posterior probability

% Convert the vector Sigma into a diagonal matrix and invert it.
sigma_diag = zeros(m,p);
for i=1:m
	sigma_diag(i,:) = 1./sqrt(gm_obj.Sigma(1,:,i));
end

% Compute V and Gamma for each p
V = zeros(m,p);
Gamma = zeros(m,p);

for k=1:p
	[V(:,k),Gamma(:,k)] = param(k,P,X_lsf(1:N,:),Y_lsf(1:N,:),gm_obj,sigma_diag); 
end

%% Save Data
save('variables128_test','V','Gamma','sigma_diag');