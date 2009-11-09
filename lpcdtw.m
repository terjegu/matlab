function [D1,D2_new,K2,index,len] = lpcdtw(d1,d2,sr)

% Calculate LPC features for both sounds
window_size = 10e-3;            % 10ms
len = floor(sr*window_size);	% samples per frame
anal = round(1.5*len);          % samples per analysis frame (overlapping)
p = 16;                         % LPC order (Fs/1000)

[D1]=lpcauto(d1,p,[len anal]);
[D2,~,K2]=lpcauto(d2,p,[len anal]);


% Construct the 'local match' scores matrix 
SM = distitar(D1,D2);
SM = SM./(max(max(SM))+0.1);

% Shortest path 
% Use dynamic programming to find the lowest-cost path between the 
% opposite corners of the cost matrix
% Note that we use 1-SM because dp will find the *lowest* total cost
[p,q,~] = dp2(1-SM);

% Calculate the frames in D2 that are indicated to match each frame
% in D1, so we can resynthesize a warped, aligned version
m = length(D1);
index = zeros(m,1);
for i = 1:m
    index(i) = q(find(p >= i,1));
end
D2_new = D2(index,:);

end