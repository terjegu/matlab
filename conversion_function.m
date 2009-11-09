%% CONVERSION FUNCTION
% Terje Gundersen 13.10.2009
close all;
clear all;
load 'variables';
load 'gmm';

%% Read file
[x,fs]=wavread('data/t03s000228.wav');

window_size = 10e-3;            % 10ms
len = floor(fs*window_size);	% samples per frame
anal = round(1.5*len);          % samples per analysis frame (overlapping)
p = 16;                         % LPC order (Fs/1000)

X_lpc= lpcauto(x,p,[len anal]); % LPC matrix

%% Convert LPC to LSF
fn = length(X_lpc);
X_lsf = zeros(fn,p);
for i=1:fn
    X_lsf(i,:) = poly2lsf(X_lpc(i,:));
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
X_lpc_conv = zeros(fn,p+1);              % LSF to LPC
for i=1:fn
    X_lpc_conv(i,:) = lsf2poly(X_conv(i,:));
end

overlap = anal-len;                       % end frame1 - start frame2 + 1
X_s = split(x,len,floor(overlap/2));    % Vector to matrix
e2 = lpcfilt(X_s,X_lpc);                % error signal
X2 = lpcifilt2(e2,X_lpc_conv);          % reconstructed matrix
x2 = concat(X2,len,floor(overlap/2));	% matrix to vector

%% Write to file
wavwrite(x2,fs,'data/test_conv.wav')


%% Plot complete signal
[y,fs_y]=wavread('data/t01s000228.wav');	% target

NFFT = pow2(nextpow2(length(x)));
f = fs/2*linspace(0,1,NFFT/2+1);
F_x = abs(fft(x,131072));
F_x2 = abs(fft(x2,131072));
F_y = abs(fft(y,131072));

% Converted
figure(1)
subplot(211);
plot(x2);
title('Converted, time domain');
subplot(212);
plot(f,F_x2(1:NFFT/2+1));
title('Frequency domain');

% Target
figure(2)
subplot(211);
plot(y,'r');
title('Target, time domain');
subplot(212);
plot(f,F_y(1:NFFT/2+1),'r');
title('Frequency domain');

% Source
figure(3)
subplot(211);
plot(x,'g');
title('Source, time domain');
subplot(212);
plot(f,F_x(1:NFFT/2+1),'g')
title('Frequency domain');

%% Plot one lpc frame
[~,Y_lpc] = lpcdtw(x,y,fs);

frame_num = 4;
N = length(X2(frame_num,:));
NFFT = pow2(nextpow2(N));
t = (1:N)/fs*1000;
frame = (N*frame_num+1:N*(frame_num+1));

% Converted
[X2_freqz,f] = freqz(X_lpc_conv(frame_num,:),1,NFFT,fs);
figure(1)
subplot(211);
plot(t,x2(frame));
title('Converted, time domain');
xlabel('t [ms]');
subplot(212);
plot(f/1000,abs(X2_freqz));
title('Frequency domain');
xlabel('f [kHz]');

% Target
[Y_freqz,f] = freqz(Y_lpc(frame_num,:),1,NFFT,fs);
figure(2)
subplot(211);
plot(t,y(frame),'r'); % y before dtw, not correct
title('Target, time domain');
xlabel('t [ms]');
subplot(212);
plot(f/1000,abs(Y_freqz));
title('Frequency domain');
xlabel('f [kHz]');

% Source
[X_freqz,f] = freqz(X_lpc(frame_num,:),1,NFFT,fs);
figure(3)
subplot(211);
plot(t,x(frame),'g');
title('Source, time domain');
xlabel('t [ms]');
subplot(212);
plot(f/1000,abs(X_freqz));
title('Frequency domain');
xlabel('f [kHz]');
