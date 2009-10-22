%% DTW ALIGN
close all;
clear all;

%% Load two speech waveforms of the same utterance (from TIMIT)
[d1,sr] = wavread('data/t01s000228.wav');
[d2,sr2] = wavread('data/t03s000228.wav');

%% Calculate LPC features for both sounds
[D1,E1,Fx1] = split_lpc(d1,sr,12);
[D2,E2,Fx2] = split_lpc(d2,sr,12);

%% Construct the 'local match' scores matrix 
SM = distitar(D1,D2);
SM = SM./(max(max(SM))+0.1);

figure(1)
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

%% Calculate the frames in D2 that are indicated to match each frame
% in D1, so we can resynthesize a warped, aligned version
[m,n] = size(D1);
index = zeros(m,1);
for i = 1:m
    index(i) = q(find(p >= i,1));
end
D2_new = D2(index,:);

%% reconstruct
% window_size = 10e-3; % 20ms
% inc = sr * window_size; % samples per frame
% 
% % Source vector
% N_x = length(d2);
% N_loop = floor(N_x/inc);
% F_x2 = zeros(N_loop,inc);
% for i=1:N_loop
%    F_x2(i,:) = d2(1+inc*(i-1):inc*i); 
% end
% 
% if N_loop < N_x/inc
%    F_x2(N_loop+1,:) = zeros(inc,1);
%    F_x2(N_loop+1,1:(N_x-inc*N_loop)) = d2(inc*N_loop+1:end); 
% end
% 
% Fx_new = F_x2(index,:);
% 
% [n,m] = size(Fx2);
% error = zeros(n,m);
% for i=1:n
%     error(i,:) = filter(D2(i,:),1,Fx2(i,:));
% end
% 
% % inverse filter with error to get y
% [n,m] = size(Fx_new);
% Y = zeros(n,m);
% error2 = error(index,:);
% for i=1:n
%     Y(i,:) = filter(1,D2_new(i,:),error2(i,:));
% end
% 
% X_final = [];
% for i=1:n
%    X_final = [X_final; Y(i,:)']; 
% end
%% Write file
% figure(3)
% plot(X_final)
% figure(4)
% plot(d2)
% wavwrite(X_final,sr,'data/test_dtw.wav')