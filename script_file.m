%% SCRIPT FILE
% Terje Gundersen 15.09.2009
% Last updated 22.09.2009
close all;
clear all;

%% Read file
[x,f_s]=wavread('data/000358_JF.wav');


%% Split file in 20ms segments
window_size = 20e-3; % 20ms
inc = f_s * window_size; % samples per frame
N_loop = floor(length(x)/inc);

F = [];
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


% Subtract estimated signal from real signal
error = zeros(m,n);
for i=1:m
    error(i,:) = filter([1 X_lpc(i,2:end)], 1, F(i,:));
end


%% Transformation LPC
X_rc = zeros(m,lpc_order);
for i=1:m
    X_rc(i,:) = poly2rc(X_lpc(i,:));
end
% status = lpcconv('ar','rf',X_lpc,Y_lpc,12);


%% Transform features


%% Transform LPC back
Y_lpc = zeros(m,lpc_order+1);
for i=1:m
    Y_lpc(i,:) = rc2poly(X_rc(i,:));
end


%% Inverse filter
x_rest = zeros(m,n);
for i=1:m
    x_rest(i,:) = filter(1, [1 Y_lpc(i,2:end)], error(i,:));
end


%% Concatenate matrix to one vector
y = [];
for i=1:m
   y = [y x_rest(i,:)]; 
end
y = y';


%% Write file
wavwrite(y,f_s,'data/000358_JF_out.wav')