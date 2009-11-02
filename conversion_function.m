%% CONVERSION FUNCTION
% Terje Gundersen 13.10.2009
close all;
clear all;
load 'variables';
load 'gmm';

%% Read file
[y,fs]=wavread('data/t01s000228.wav');
[x,fs_y]=wavread('data/t03s000228.wav');

[Y,X,K1,index,len] = lpcdtw(y,x,fs); % returns time aligned lp coefficients

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

%% Conversion function
X_conv = zeros(fn,p);
for i=1:fn
    for k=1:p
        X_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_lsf(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end


%% Reconstruct 
X_lpc_new = zeros(fn,p+1);              % LSF to LPC
for i=1:fn
    X_lpc_new(i,:) = lsf2poly(X_conv(i,:));
end

overlap=K1(1,2)-K1(2,1)+1;              % end frame1- start frame2 + 1
X_s = split(x,len,floor(overlap/2));    % Vector to matrix
X_new = X_s(index,:);                   % DTW resulting matrix

e2 = lpcfilt(X_new,X_lpc_new);          % error signal
X2 = lpcifilt2(e2,X_lpc_new);           % reconstructed matrix
x2 = concat(X2,len,floor(overlap/2));	% matrix to vector

%% Plot and write to file
% NFFT = pow2(nextpow2(length(x)));
% f = fs/2*linspace(0,1,NFFT/2+1);
% F_x = abs(fft(x,131072));
% F_x2 = abs(fft(x2,131072));
% F_y = abs(fft(y,131072));
% 
% figure(1)
% subplot(211);
% plot(x2);
% title('Converted, time domain');
% subplot(212);
% plot(f,F_x2(1:NFFT/2+1));
% title('Frequency domain');
% 
% figure(2)
% subplot(211);
% plot(y,'r');
% title('Target, time domain');
% subplot(212);
% plot(f,F_y(1:NFFT/2+1),'r');
% title('Frequency domain');
% 
% figure(3)
% subplot(211);
% plot(x,'g');
% title('Source, time domain');
% subplot(212);
% plot(f,F_x(1:NFFT/2+1),'g')
% title('Frequency domain');

wavwrite(x2,fs,'data/test.wav')