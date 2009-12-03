%% TEST
close all;
clear all;

%%
[pm_x,~] = textread('data/t03s000228.txt','%f%f');
[x,fs] = wavread('data/t03s000228.wav');  % source
y = wavread('data/t01s000228.wav');       % target
[pm_y,~] = textread('data/t01s000228.txt','%f%f');
pm_x = round(pm_x*fs)+1;
pm_y = round(pm_y*fs)+1;

%% NY METODE
% lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
% analx = [pm_x(2); pm_x(3:nfx-1)-pm_x(1:nfx-3); size(x,1)-pm_x(nfx-1)]-1;
% skipx = zeros(nfx-1,1);
% tfx = [lenx analx skipx];
% 
% leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
% analy = [pm_y(2); pm_y(3:nfy-1)-pm_y(1:nfy-3); size(y,1)-pm_y(nfy-1)]-1;
% skipy = zeros(nfy-1,1);
% tfy = [leny analy skipy];
% 
% X_lpc = lpcauto(x,p,tfx);~
% Y_lpc = lpcauto(y,p,tfy);


%% DTW
% SM = distitar(X_lpc,Y_lpc); % Construct the scores matrix 
% SM = SM./(max(SM(:))+0.1);  % Normalise
% 
% [p1,q1,~] = dp2(1-SM);    % Dynamic programming to find the lowest-cost path
% 
% m = length(X_lpc);
% index = zeros(m,1);     % Matching indecies
% for i = 1:m
%     index(i) = q1(find(p1 >= i,1));
% end
% Y_lpc = Y_lpc(index,:); % Update Y_lpc

%%
% [X_lpc,Y_lpc,~,tfx] = lpcdtw(x,y,pm_x,pm_y);
e1 = lpcifilt2(x,X_lpc,tfx);     % Exitation
% x_y = lpcfilt2(e1,Y_lpc,tfx);    % Synthesis

%% PLOTS

x_y = x_y-mean(x_y);
% x_y(x_y>0.5) = 0.5;
% x_y(x_y<-0.5) = -0.5;
% 
figure(1)
subplot(3,1,1)
plot(x)
title('Source')
subplot(3,1,2)
plot(y)
title('Target')
subplot(3,1,3)
plot(x_y)
title('Converted')
%%


wavwrite(x_y,fs,'data/test_2.wav')
% soundsc(x_y,fs);
