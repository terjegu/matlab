%% SCRIPT FILE
% Terje Gundersen 15.09.2009
close all;
clear all;


[x,f_s,WMODE,FIDX]=readwav('data/000358_JF.wav');

%% LPC
lpc_order = 12;
window_size = 20e-3; % 20ms
inc = f_s * window_size; % number of frames

F = enframe(x,inc); % Splits the singla into frames of 20ms
[m,n] = size(F);

% Make AR(12) coefficients
X_lpc = zeros(m,lpc_order + 1);
for i=1:m
    X_lpc(i,:) = lpc(F(i,:),lpc_order);
end

% Estimate signal from LPC
x_est = filter([0 -X_lpc(1,2:end)], 1, F(1,:));

% Subtract estimated signal from original signal
error = F(1,:) - x_est;


figure(1)
plot(F(1,:),'r')
hold on;
plot(x_est)
plot(abs(error),'g')
title('Error signal before transformation');

%% Transformation LPC
Y_lpc = zeros(m,lpc_order);
% status = lpcconv('rr','rf',X_lpc,Y_lpc,12);


%% Transform features


%% Transform LPC back


%% x_reset
x_rest = filter([0 -Y_lpc(1,2:end)], 1, F(1,:));

F_out = error + x_rest;

% concat
