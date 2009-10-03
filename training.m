%% SCRIPT FILE
% Terje Gundersen 15.09.2009
close all;
clear all;

%% Read file
[x,f_s]=wavread('data/000358_JF.wav');


%% Split file in 20ms segments
window_size = 20e-3; % 20ms
inc = f_s * window_size; % samples per frame
N_loop = floor(length(x)/inc);

F = zeros(N_loop,inc);
for i=1:N_loop
   F(i,:) = x(1+inc*(i-1):inc*i); 
end

if N_loop < length(x)/inc
   F(N_loop+1,:) = x(inc*N_loop+1:end); 
end

% Dimensions
[n,p_lpc] = size(F);


%% Make AR(p) coefficients
p = 12; % LPC order
X_lpc = lpc(F',p);


%% Transformation LPC --> RC
X_rc = zeros(n,p);
for i=1:n
    X_rc(i,:) = poly2rc(X_lpc(i,:));
end
% status = lpcconv('ar','rf',X_lpc,Y_lpc,12);


%% EM algorithm
m = 5; % Number of mixture models
gm_obj = gmdistribution.fit(X_rc,m,'CovType','diagonal'); % EM-alg


%% Matrices
P = posterior(gm_obj,X_rc); % Posterior probability

% Convert the vector Sigma into a diagonal matrix and invert it.
sigma = zeros(p,p,m);
for i=1:m
    for j=1:p
        sigma(j,j,i) = 1./gm_obj.Sigma(1,j,i);
    end
end

% Calculate the matrix D = P(C|x) * (x-mu)^T * Sigma^-1
D = zeros(n,m*p);
for i=1:n
    for j=1:m
        D(i,1+(j-1)*p:j*p) =...
            P(i,j)*(X_rc(i,:)-gm_obj.mu(j,:)) * sigma(:,:,j);
    end
end


%% Conversion Function
% param = inv([P';D']*[P D])*[P';D']*y
% V = param(1:m,:);
% Gamma = param(m:end,:);