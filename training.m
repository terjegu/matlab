%% SCRIPT FILE
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
X_rc = zeros(n,p);
for i=1:n
    X_rc(i,:) = poly2rc(X_lpc(i,:));
end
Y_rc = zeros(n,p);
for i=1:n
    Y_rc(i,:) = poly2rc(Y_lpc(i,:));
end


%% EM algorithm
m = 8; % Number of mixture models
gm_obj = gmdistribution.fit(X_rc,m,'CovType','diagonal'); % EM-alg


%% Matrices
P = posterior(gm_obj,X_rc); % Posterior probability

% Convert the vector Sigma into a diagonal matrix and invert it.
sigma_inv = zeros(p,p,m);
for i=1:m
    for j=1:p
        sigma_inv(j,j,i) = 1./gm_obj.Sigma(1,j,i);
    end
end

% Calculate the matrix D = P(C|x) * (x-mu)^T * Sigma^-1
D = zeros(n,p*m);
for i=1:n
    for j=1:m
        D(i,1+(j-1)*p:j*p) =...
            P(i,j)*(X_rc(i,:)-gm_obj.mu(j,:)) * sigma_inv(:,:,j);
    end
end


%% Conversion variables
param = ([P';D']*[P D])\[P';D']*Y_rc;
V = param(1:m,:);
Gamma = param((m+1):end,:);


%% Conversion function
X_conv = zeros(n,p);

for i=1:n
    for j=1:m
        var = Gamma(1+(j-1)*p:j*p,:)*sigma_inv(:,:,j)*...
            (X_rc(i,:)-gm_obj.mu(j,:))';
        brackets = var' + V(j,:);
        X_conv(i,:) = X_conv(i,:)+P(i,j)*brackets;
    end
end

%% test
% test = zeros(p,p,m);
% for i=1:p
%     for j=1:m
%         test(:,:,j) = Gamma(1+(j-1)*p:j*p,:)*sigma_inv(:,:,j);
%     end
% end
%% reconstruct
% 
% X_lpc_new = zeros(n,p+1);
% for i=1:n
%     X_lpc_new(i,:) = rc2poly(X_conv(i,:));
% end
% 
% error = zeros(n,frame_length);
% for i=1:n
%     error(i,:) = filter([1 X_lpc(i,2:end)], 1, F_x(i,:));
% end
% 
% X_rest = zeros(n,frame_length);
% for i=1:n
%     X_rest(i,:) = filter(1, [1 X_lpc_new(i,2:end)], error(i,:));
% end
% 
% X_final = [];
% for i=1:n
%    X_final = [X_final; X_rest(i,:)']; 
% end
% % 
% 
% wavwrite(X_final,f_s,'data/test.wav')