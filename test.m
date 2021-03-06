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
[X_lpc,Y_lpc,index] = lpcdtw(x,y,pm_x,pm_y);
e_x = lpcifilt2(x,X_lpc,pm_x);     % Exitation
e_y = lpcifilt2(y,Y_lpc,pm_y);     % Exitation
% x_y = lpcfilt2(e_x,Y_lpc(index,:),pm_x);    % Synthesis
[y_yx,exct]=psolasynth(length(e_x),e_y,pm_x,pm_y,length(X_lpc),Y_lpc,index);

% e2 = 0.02*randn(length(e_x),1);
% e_x = lpcfilt2(e2,X_lpc,pm_x);    % Synthesis
% e_y = lpcfilt2(e2,Y_lpc,pm_x);    % Synthesis

%% PLOTS

x_y = x_y-mean(x_y);
y_yx = y_yx-mean(y_yx);
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
plot(y_yx)
title('Converted')
%%


wavwrite(y_yx,fs,'data/test_psola.wav')
% soundsc(x_y,fs);
% soundsc(y_yx,fs);
