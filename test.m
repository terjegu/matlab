close all;
clear all;

[x,fs]=wavread('data/t03s000228.wav');  % target
y=wavread('data/t01s000228.wav');       % target
window_size = 10e-3;                    % 20ms
len = floor(fs*window_size);            % samples per frame


[X_lpc,Y_lpc] = lpcdtw(x,y,fs);


X_s = split(x,len);                  % Vector to matrix
e = lpcfilt(X_s,X_lpc);              % error signal
X2 = lpcifilt2(e,Y_lpc);             % reconstructed matrix
temp = X2';
x2 = temp(:);                        % matrix to vector

figure(1)
subplot(2,1,1)
plot(x2)
title('Converted')
subplot(2,1,2)
plot(x)
title('Source')

% wavwrite(x2,fs,'data/test.wav')