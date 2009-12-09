%% DTW ALIGN
% Terje Gundersen 01.11.2009
close all;
clear all;

%% Load two speech waveforms of the same utterance
[pm_x,~] = textread('data/t03s000228.txt','%f%f');
[d1,sr] = wavread('data/t03s000228.wav');  % source
d2 = wavread('data/t01s000228.wav');       % target
[pm_y,~] = textread('data/t01s000228.txt','%f%f');
pm_x = pm_x*sr;
pm_y = pm_y*sr;
%% Calculate LPC features for both sounds

p = 16; % LPC order (Fs/1000)
nfx = length(pm_x);
nfy = length(pm_y);

lenx = [pm_x(1); pm_x(2:nfx-1)-pm_x(1:nfx-2)];
analx = max(256*ones(nfx-1,1),[pm_x(2);pm_x(3:nfx-1)-pm_x(1:nfx-3);length(d1)-pm_x(nfx-2)]-1);
skipx = zeros(nfx-1,1);
tfx = [lenx analx skipx];

leny = [pm_y(1); pm_y(2:nfy-1)-pm_y(1:nfy-2)];
analy = max(256*ones(nfy-1,1),[pm_y(2);pm_y(3:nfy-1)-pm_y(1:nfy-3);length(d2)-pm_y(nfy-2)]-1);
skipy = zeros(nfy-1,1);
tfy = [leny analy skipy];

[D1,E1,K1] = lpcauto(d1,p,tfx);
[D2,E2,K2] = lpcauto(d2,p,tfy);

%% Construct the 'local match' scores matrix 
SM = distitar(D1,D2);
SM = SM./(max(max(SM))+0.001); % scale values to [0 0.9999]

figure(1)
subplot(121)
imagesc(SM);
colormap(1-gray);

%% Shortest path 
% Use dynamic programming to find the lowest-cost path
% Note that we use 1-SM because dp will find the *lowest* total cost
[p,q,C] = dp2(1-SM);

% Overlay the path on the local similarity matrix
hold on; 
plot(q,p,'r'); 
hold off

subplot(122)
imagesc(C)
hold on; 
plot(q,p,'r');
hold off
%% Calculate the frames in D2 that are indicated to match each frame
% in D1, so we can resynthesize a warped, aligned version
[m,n] = size(D1);
index = zeros(m,1);
for i = 1:m
    index(i) = q(find(p >= i,1));
end
D2_new = D2(index,:);
