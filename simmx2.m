function M = simmx2(A,B)
% M = simmx2(A,B) calculate a sim matrix 
% between MFCC feature matrices A and B.

[m,na] = size(A);
nb = size(B,2);

A_cep = zeros(m,na);
B_cep = zeros(m,nb);

for i=1:m
    A_cep(i,:) = ifft(log(A(i,:)));
    B_cep(i,:) = ifft(log(B(i,:)));
end

M = zeros(na,nb);
for i = 1:na
	for j = 1:nb
        M(i,j) = (norm(A_cep(:,i)-B_cep(:,j)))/m;
	end
end

end






% ncA = size(A,2);
% ncB = size(B,2);
% M = zeros(ncA, ncB);
%for i = 1:ncA
%  for j = 1:ncB
%    % normalized inner product i.e. cos(angle between vectors)
%    M(i,j) = (A(:,i)'*B(:,j))/(EA(i)*EB(j));
%  end
%end

% M = (A'*B)./(EA'*EB);






%% SIMMX2
% close all;
% clear all;

%% Load two speech waveforms of the same utterance (from TIMIT)
% [x,fs] = wavread('data/t01s000228.wav');
% [y,fs2] = wavread('data/t03s000228.wav');

% % Enframe
% 
% window_size = 20e-3; % 20ms
% inc = fs * window_size; % samples per frame
% 
% 
% Source vector
% N_loop = floor(length(x)/inc);
% F_x = zeros(N_loop,inc);
% for i=1:N_loop
%    F_x(i,:) = x(1+inc*(i-1):inc*i); 
% end
% 
% if N_loop < length(x)/inc
%    F_x(N_loop+1,:) = zeros(inc,1);
%    F_x(N_loop+1,1:(length(x)-inc*N_loop)) = x(inc*N_loop+1:end); 
% end
% F_x = F_x';
% 
% Target vector
% N_loop = floor(length(y)/inc);
% F_y = zeros(N_loop,inc);
% for i=1:N_loop
%    F_y(i,:) = y(1+inc*(i-1):inc*i); 
% end
% 
% if N_loop < length(y)/inc
%    F_y(N_loop+1,:) = zeros(inc,1);
%    F_y(N_loop+1,1:(length(y)-inc*N_loop)) = y(inc*N_loop+1:end); 
% end
% F_y = F_y';
% 
% 
% %
% x_cep = zeros(size(F_x));
% y_cep = zeros(size(F_y));
% for i=1:size(F_x,2)
%    [temp,n,x_cep(:,i)] = cceps(F_x(:,i)); 
% end
% for i=1:size(F_y,2)
%    y_cep(:,i) = cceps(F_y(:,i)); 
% end

% test = cceps(F_x);
% back = icceps(test);
% plot(test(1:20))


%%

% wavwrite(back,sr,'data/test_ccep.wav')
