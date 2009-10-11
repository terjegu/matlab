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

% Compute V and Gamma for each p and apply transfer function
X_conv = zeros(n,p);
for i=1:n
    for k=1:p
        [V,Gamma] = param(k,n,m,P,X_lsf,Y_lsf,gm_obj,sigma_diag); 
        X_conv(i,k) = sum(P(i,:).*(Gamma(:).*(X_lsf(i,k)-gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:))');
    end
end


%% normalise
% for i=1:n
%     for j=1:p
%         if X_conv(i,j)>pi
%             X_conv(i,j)=pi*0.999;
%         elseif X_conv(i,j)<0
%             X_conv(i,j)=0;
%         end
%     end
% end


%% reconstruct
% LSF to LPC
X_lpc_new = zeros(n,p+1);
for i=1:n
    X_lpc_new(i,:) = lsf2poly(X_conv(i,:));
end
% 
% Extract error signal from x
error = zeros(n,frame_length);
for i=1:n
    error(i,:) = filter([1 X_lpc(i,2:end)], 1, F_x(i,:));
end

% inverse filter with error to get y
Y = zeros(n,frame_length);
for i=1:n
    Y(i,:) = filter(1, [1 X_lpc_new(i,2:end)], error(i,:));
end

X_final = [];
for i=1:n
   X_final = [X_final; Y(i,:)']; 
end

% X_final = X_final./max(X_final);
% for i=1:length(X_final)
%     if X_final(i) > 0.5
%         X_final(i) = 0.5;
%     elseif X_final(i) < -0.5
%         X_final(i) = -0.5;
%     end
% end
% close all;
figure(1)
plot(X_final);
figure(2)
plot(x,'r');
wavwrite(X_final,f_s,'data/test.wav')