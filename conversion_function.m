%% CONVERSION FUNCTION
% Terje Gundersen 13.10.2009
close all;
clear all;
load 'variables64';
load 'gmm64';

%% Read file
[x,fs]=wavread('data/t03s000228.wav');

window_size = 10e-3;            % 20ms
len = floor(fs*window_size);	% samples per frame
anal = round(1*len);          % samples per analysis frame (overlapping)
p = 16;                         % LPC order (Fs/1000)

X_lpc = lpcauto(x,p,[len anal]); % LPC matrix

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
X_lpc_conv = zeros(fn,p+1);             % LSF to LPC
for i=1:fn
%     X_conv(X_conv<0)=0;
%     X_conv(X_conv>=pi)=pi*0.999;
    X_lpc_conv(i,:) = lsf2poly(X_conv(i,:));
end

% overlap = anal-len;                   % end frame1 - start frame2 + 1,
% floor(overlap/2)

% X_s = split_overlap(x,len,anal-len);  % Vector to matrix
X_s = split(x,len);                     % Vector to matrix
e = lpcfilt(X_s,X_lpc);                 % error signal
X2 = lpcifilt2(e,X_lpc_conv);           % reconstructed matrix
temp = X2';
x2 = temp(:);                           % matrix to vector
% x2 = concat_overlap(X2,len,anal-len); % matrix to vector

%% Write to file
% x2(x2>=1) = 0.999;                    % Prevent clipping
% x2(x2<-1) = -1;
% 
% wavwrite(x2,fs,'data/test_conv.wav')

%% Plot complete signal
y=wavread('data/t01s000228.wav');   % target
% temp = e';
% e2 = temp(:);                       % error

t_x = (1:length(x))/fs;
t_x2 = (1:length(x2))/fs;
t_y = (1:length(y))/fs;
% t_e2 = (1:length(e2))/fs;

NFFT = pow2(nextpow2(length(x)));
f = fs/2/1000*linspace(0,1,NFFT/2+1);
F_x = log10(abs(fft(x,NFFT)));
F_x2 = log10(abs(fft(x2,NFFT)));
F_y = log10(abs(fft(y,NFFT)));
% F_e2 = log10(abs(fft(e2,NFFT)));

% Converted
figure(1)
subplot(311);
plot(t_x,x,'g');
title('Source, time domain');
subplot(312);
plot(t_y,y,'r');
title('Target, time domain');
subplot(313);
plot(t_x2,x2);
title('Converted, time domain');
xlabel('t [s]');

% % Target
% figure(2)
% subplot(311);
% plot(f,F_x(1:NFFT/2+1),'g')
% title('Source');
% ylabel('dB');
% subplot(312);
% plot(f,F_y(1:NFFT/2+1),'r');
% title('Target');
% ylabel('dB');
% subplot(313);
% plot(f,F_x2(1:NFFT/2+1));
% title('Converted');
% xlabel('f [kHz]');
% ylabel('dB');


% % Error
% figure(4)
% subplot(211);
% plot(t_e2,e2,'k');
% title('Error, time domain');
% xlabel('t [s]');
% subplot(212);
% plot(f,F_e2(1:NFFT/2+1),'k');
% title('Frequency domain');
% xlabel('f [kHz]');
% ylabel('dB');

%% Plot one lpc frame
[y,fs_y]=wavread('data/t01s000228.wav'); % target
[~,Y_lpc,index] = lpcdtw(x,y,fs);
% Y_s = split(y,len,0);       % Vector to matrix
% Y_s = Y_s(index);
% temp = Y_s';
% y2 = temp(:);                   % matrix to vector

frame_num = 120;
N = length(X2(frame_num,:));
NFFT = pow2(nextpow2(N));
t = (1:N)/fs*1000;
frame = (N*frame_num+1:N*(frame_num+1));

[X_freqz,f_x] = freqz(1,X_lpc(frame_num,:),NFFT,fs);
[Y_freqz,f_y] = freqz(1,Y_lpc(frame_num,:),NFFT,fs);
[X2_freqz,f_x2] = freqz(1,X_lpc_conv(frame_num,:),NFFT,fs);

% 
% figure(5)
% subplot(311);
% plot(t,x(frame),'g');
% title('Source, time domain');
% xlabel('t [ms]');
% subplot(312);
% plot(t,y(frame),'r'); % y before dtw, not correct
% title('Target, time domain');
% xlabel('t [ms]');
% subplot(313);
% plot(t,x2(frame));
% title('Converted, time domain');
% xlabel('t [ms]');


% subplot(313);
% plot(t,e(frame));
% title('Exitation');
% xlabel('f [kHz]');


figure(6)
subplot(311);
plot(f_x/1000,log10(abs(X_freqz)),'g');
title('Source');
ylabel('dB');
subplot(312);
plot(f_y/1000,log10(abs(Y_freqz)),'r');
title('Target');
ylabel('dB');
subplot(313);
plot(f_x2/1000,log10(abs(X2_freqz)));
title('Converted');
xlabel('f [kHz]');
ylabel('dB');

% 
% figure(7)
% hold on;
% plot(f_x/1000,log10(abs(X_freqz)),'g');
% plot(f_y/1000,log10(abs(Y_freqz)),'r');
% plot(f_x2/1000,log10(abs(X2_freqz)));
% title('Source');
% legend('Source','Target','Converted');
% xlabel('f [kHz]');
% ylabel('dB');
% hold off;

%% Compare GMM
% save('x2_64','X2_freqz');

% figure(8)
% load('x2_16');
% subplot(311);
% plot(f_x2/1000,log10(abs(X2_freqz)),'g');
% title('16 GMM');
% % xlabel('f [kHz]');
% ylabel('dB');
% load('x2_32');
% subplot(312);
% plot(f_x2/1000,log10(abs(X2_freqz)),'r');
% title('32 GMM');
% % xlabel('f [kHz]');
% ylabel('dB');
% load('x2_64');
% subplot(313);
% plot(f_x2/1000,log10(abs(X2_freqz)));
% title('64 GMM');
% xlabel('f [kHz]');
% ylabel('dB');
