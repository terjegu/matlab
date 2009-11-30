%% TEST
close all;
clear all;

%%
[pm_x,~] = textread('data/t03s000228.txt','%f%f');
[x,fs] = wavread('data/t03s000228.wav');  % source
y = wavread('data/t01s000228.wav');       % target
[pm_y,~] = textread('data/t01s000228.txt','%f%f');
pm_x = pm_x*fs;
pm_y = pm_y*fs;
window_size = 10e-3;                    % 10ms
p = 16;

%% Find t = [len anal skip]
N_pmx = length(pm_x);
len_x = window_size*fs;
% len_x = pm_x(1);
for i = 2:N_pmx
   temp = pm_x(i)-pm_x(i-1);
   if temp<len_x
       len_x = floor(temp);
   end
end

skip_x = zeros(N_pmx,1);
anal_x = zeros(N_pmx,1);
anal_x(1) = pm_x(2);
skip_x(1) = min(0,round(pm_x(1)-anal_x(1)/2));
for i = 2:N_pmx-1
    anal_x(i) = round(pm_x(i+1)-pm_x(i-1));
    skip_x(i) = round(pm_x(i)-anal_x(i)/2-(i-1)*len_x);
end
anal_x(N_pmx) = length(x)-pm_x(N_pmx);
skip_x(N_pmx) = round(pm_x(N_pmx)-anal_x(N_pmx)/2-(N_pmx-1)*len_x);
len_x = len_x*ones(N_pmx,1);


% skip_x = zeros(N_pmx,1);
% anal_x = 2*len_x;
% skip_x(1) = round(pm_x(1)-anal_x/2);
% for i = 2:N_pmx-1
%     skip_x(i) = round(pm_x(i)-i*anal_x/2); % pm_x(i)-anal_x/2-(i-1)*len_x
% end
% skip_x(N_pmx) = round(pm_x(N_pmx)-pm_x(N_pmx-1)-anal_x);
% len_x = len_x*ones(N_pmx,1);
% anal_x = anal_x*ones(N_pmx,1);

t_x = [len_x anal_x skip_x];
X_lpc = lpcauto(x,p,t_x);

%% Find t = [len anal skip]
N_pmy = length(pm_y);
len_y = window_size*fs;
% len_y = pm_y(1);
for i = 2:N_pmy
   temp = pm_y(i)-pm_y(i-1);
   if temp<len_y
       len_y = floor(temp);
   end
end

skip_y = zeros(N_pmy,1);
anal_y = zeros(N_pmy,1);
anal_y(1) = pm_y(2);
skip_y(1) = min(0,round(pm_y(1)-anal_y(1)/2));
for i = 2:N_pmy-1
    anal_y(i) = round(pm_y(i+1)-pm_y(i-1));
    skip_y(i) = round(pm_y(i)-anal_y(i)/2-(i-1)*len_y);
end
anal_y(N_pmy) = length(y)-pm_y(N_pmy);
skip_y(N_pmy) = round(pm_y(N_pmy)-anal_y(N_pmy)/2-(N_pmy-1)*len_y);

len_y = len_y*ones(N_pmy,1);
t = [len_y anal_y skip_y];
Y_lpc = lpcauto(y,p,t);

%% DTW
SM = distitar(X_lpc,Y_lpc); % Construct the scores matrix 
SM = SM./(max(SM(:))+0.1);  % Normalise

[p,q,~] = dp2(1-SM);    % Dynamic programming to find the lowest-cost path

m = length(X_lpc);
index = zeros(m,1);     % Matching indecies
for i = 1:m
    index(i) = q(find(p >= i,1));
end
Y_lpc = Y_lpc(index,:); % Update Y_lpc

%%
% [X_lpc,Y_lpc,index] = lpcdtw(x,y,fs);
% Y_lpc2 = lpcauto(y,p,len);
e1 = lpcifilt(x,X_lpc,pm_x-floor(window_size*fs/2));     % Exitation
x1 = lpcfilt2(e1,Y_lpc,pm_x-floor(window_size*fs/2));    % Synthesis
% X_s1 = split_overlap(x,len,anal_x-len);   % Vector to matrix
% e1 = lpcfilt(X_s1,X_lpc);                 % error signal
% X1 = lpcifilt2(e1,Y_lpc);                 % reconstructed matrix
% x1 = concat_overlap(X1,len,anal_x-len);
% X_s2 = split_pm(x,len_x(1),anal_x,skip_x);  % Vector to matrix
% X_s2 
% e2 = lpcfilt(X_s2,X_lpc);                   % error signal
% X2 = lpcifilt2(e2,Y_lpc);                   % reconstructed matrix

% x2 = concat_pm(X2,pm_y,length(y));
%%
% temp = X2';
% x2 = temp(:);                               % matrix to vector
% e3 = 0.003*randn(length(Y_lpc),len);      % noise
% X3 = lpcifilt2(e3,Y_lpc);                 % reconstructed matrix
% temp = X3';
% x3 = temp(:);    

% temp = e';
% e2 = temp(:);
% temp = e3';
% e3 = temp(:);
 
% figure(1)
% subplot(4,1,1)
% plot(x1)
% title('Converted overlap')
% subplot(4,1,2)
% plot(x2)
% title('Converted without overlap')
% subplot(4,1,3)
% plot(x3)
% title('Converted with noise')
% subplot(4,1,4)
% plot(y)
% title('Target')

figure(1)
subplot(2,1,1)
plot(x1)
title('Converted')
subplot(2,1,2)
plot(y)
title('Source')
%%
% x1(x1>max(x)) = max(x);
% x1(x1<min(x)) = min(x);
wavwrite(x1,fs,'data/test_mem.wav')

%% draw selections
% close all;
% a = zeros(1000,1);
% a(pm_x(1:7)) = 1e-3;
% b = zeros(1000,1);
% c = zeros(1000,1);
% b(start(1):stop(1)) = 5e-4;
% c(start(4):stop(4)) = 5e-4;
% 
% figure(2)
% plot(x(1:1000))
% hold on;
% stem(a,'r')
% plot(b,'g')
% plot(c,'k')