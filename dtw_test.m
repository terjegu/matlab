%% DTW ALIGN
% Terje Gundersen 01.11.2009
close all;
clear all;

%% Load two speech waveforms of the same utterance
[d1,sr] = wavread('data/t01s000228.wav');
[d2,sr2] = wavread('data/t03s000228.wav');

%% Calculate LPC features for both sounds
window_size = 10e-3; % 15ms
len = floor(sr*window_size); % samples per frame
anal = round(1.5*len); % samples per analysis frame (overlapping)
p = 16; % LPC order (Fs/1000)

[D1,E1,K1]=lpcauto(d1,p,[len anal]);
[D2,E2,K2]=lpcauto(d2,p,[len anal]);


%% Construct the 'local match' scores matrix 
SM = distitar(D1,D2);
SM = SM./(max(max(SM))+0.001); % scale values to [0 0.9999]

figure(1)
subplot(121)
imagesc(SM);
colormap(1-gray);

%% Shortest path 
% Use dynamic programming to find the lowest-cost path between the 
% opposite corners of the cost matrix
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

%% reconstruct
overlap=K1(1,2)-K1(2,1)+1; % end frame1- start frame2 + 1
d2_s = split(d2,len,floor(overlap/2));
d2_new = d2_s(index,:);

e2 = lpcfilt(d2_new,D2_new); % error signal
X2 = lpcifilt2(e2,D2_new); % reconstructed matrix
x2 = concat(X2,len,floor(overlap/2)); % matrix to vector

%% Write file
wavwrite(x2,sr,'data/test_dtw.wav')