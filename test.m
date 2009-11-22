%% TEST
close all;
clear all;

%%
[pm_x,~] = textread('data/t03s000228.txt','%f%f');
[x,fs]=wavread('data/t03s000228.wav');  % source
y=wavread('data/t01s000228.wav');       % target
[pm_y,~] = textread('data/t01s000228.txt','%f%f');
pm_x = pm_x*fs;
pm_y = pm_y*fs;
window_size = 10e-3;                    % 10ms
p = 16;

%% Find t = [len anal skip]
N_pm = length(pm_x);
len = window_size*fs;
% len = pm_x(1);
for i = 2:N_pm
   temp = pm_x(i)-pm_x(i-1);
   if temp<len
       len = floor(temp);
   end
end

skip = zeros(N_pm,1);
anal = zeros(N_pm,1);
anal(1) = pm_x(2);
skip(1) = min(0,round(pm_x(1)-anal(1)/2));

for i = 2:N_pm-1
    anal(i) = round(pm_x(i+1)-pm_x(i-1));
    skip(i) = round(pm_x(i)-anal(i)/2-(i-1)*len);
end
anal(N_pm) = length(x)-pm_x(N_pm);
skip(N_pm) = round(pm_x(N_pm)-anal(N_pm)/2-(N_pm-1)*len);
len = len*ones(N_pm,1);


% skip = zeros(N_pm,1);
% anal = 2*len;
% skip(1) = round(pm_x(1)-anal/2);
% for i = 2:N_pm-1
%     skip(i) = round(pm_x(i)-i*anal/2); % pm_x(i)-anal/2-(i-1)*len
% end
% skip(N_pm) = round(pm_x(N_pm)-pm_x(N_pm-1)-anal);
% len = len*ones(N_pm,1);
% anal = anal*ones(N_pm,1);

t_x = [len anal skip];
X_lpc = lpcauto(x,p,t_x);

%% Find t = [len anal skip]
N_pm = length(pm_y);
len = window_size*fs;
% len = pm_y(1);
for i = 2:N_pm
   temp = pm_y(i)-pm_y(i-1);
   if temp<len
       len = floor(temp);
   end
end

skip = zeros(N_pm,1);
anal = zeros(N_pm,1);
anal(1) = pm_y(2);
skip(1) = min(0,round(pm_y(1)-anal(1)/2));

for i = 2:N_pm-1
    anal(i) = round(pm_y(i+1)-pm_y(i-1));
    skip(i) = round(pm_y(i)-anal(i)/2-(i-1)*len);
end

anal(N_pm) = length(y)-pm_y(N_pm);
skip(N_pm) = round(pm_y(N_pm)-anal(N_pm)/2-(N_pm-1)*len);

len = len*ones(N_pm,1);
t = [len anal skip];
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
len = window_size*fs;
% [X_lpc,Y_lpc,index] = lpcdtw(x,y,fs);
% Y_lpc2 = lpcauto(y,p,len);

% X_s1 = split_overlap(x,len,anal-len);                     % Vector to matrix
% e1 = lpcfilt(X_s1,X_lpc);                 % error signal
% X1 = lpcifilt2(e1,Y_lpc);                % reconstructed matrix
% x1 = concat_overlap(X1,len,anal-len);
X_s2 = split(x,len);                     % Vector to matrix
e2 = lpcfilt(X_s2,X_lpc);                 % error signal
X2 = lpcifilt2(e2,Y_lpc);                % reconstructed matrix
temp = X2';
x2 = temp(:);                           % matrix to vector
% e3 = 0.003*randn(length(Y_lpc),len);         % noise
% X3 = lpcifilt2(e3,Y_lpc);                % reconstructed matrix
% temp = X3';
% x3 = temp(:);    

% temp = e';
% e2 = temp(:);
% temp = e3';
% e3 = temp(:);
 
figure(1)
% subplot(4,1,1)
% plot(x1)
% title('Converted overlap')
% subplot(4,1,2)
plot(x2)
% title('Converted without overlap')
% subplot(4,1,3)
% plot(x3)
% title('Converted with noise')
% subplot(4,1,4)
% plot(y)
% title('Target')

%%
% wavwrite(x1,fs,'data/test.wav')
wavwrite(x2,fs,'data/test3.wav')