%% CONVERSION FUNCTION
% Terje Gundersen 13.10.2009
close all;
clear all;

%% Read file
[x,f_s]=wavread('data/t01s000228.wav');
[y,fs_y]=wavread('data/t03s000228.wav');


%% Split file in 20ms segments
p = 12; % LPC order
[X_lsf,Y_lsf,n,frame_length,X_lpc,F_x] = makelsf(x,y,f_s,p);


%% Conversion function
load 'variables';

X_conv = zeros(n,p);
for i=1:n
    for k=1:p
        X_conv(i,k) = sum(P(i,:).*(Gamma(:,k).*(X_lsf(i,k)-...
            gm_obj.mu(:,k)).*sigma_diag(:,k)+V(:,k))');
    end
end


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

figure(1)
plot(X_final);
title('Converted');
figure(2)
plot(y,'r');
title('Target');
figure(3)
plot(x,'g');
title('Source');
% wavwrite(X_final,f_s,'data/test.wav')