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
   test(N_loop+1,:) = x(inc*N_loop:end); 
end

% Dimensions
[m,n] = size(F);


%% Make AR(12) coefficients
lpc_order = 12;
X_lpc = lpc(F',lpc_order);


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


%% Matrices
P = posterior(gm_obj,X_rc); % Posterior probability

% Convert the vector Sigma into a diagonal matrix and invert it.
sigma = zeros(lpc_order,lpc_order,N_class);
for i=1:N_class
    for j=1:lpc_order
        sigma(j,j,i) = 1./gm_obj.Sigma(1,j,i);
    end
end

% Calculate the matrix D = P(C|x) * (x-mu)^T * Sigma^-1
D = zeros(m,N_class*lpc_order);
for i=1:m
    for j=1:N_class
        D(i,1+(j-1)*lpc_order:j*lpc_order) =...
            P(i,j)*(X_rc(i,:)-gm_obj.mu(j,:)) * sigma(:,:,j);
    end
end


%% Conversion Function
% param = inv([P';D']*[P D])*[P';D']*y
% V = param(1:N_class,:);
% Gamma = param(N_class:end,:);