function [p,q,D] = dp2(M)
% [p,q,D] = dp2(M) 
%    Use dynamic programming to find a min-cost path through matrix M.
%    Return state sequence in p,q
%    This version has limited slopes [2/1] .. [1/2]
% 2003-03-15 dpwe@ee.columbia.edu

% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

[r,c] = size(M);

% costs
D = zeros(r+1, c+1);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
D(2:(r+1), 2:(c+1)) = M;

% traceback
N_i = r+1;
N_j = c+1;
lim_1 = 2/3*(N_j-N_i/2);
lim_2 = 2/3*(2*N_i-N_j);
open_ends = 20;
phi = zeros(N_i,N_j);

for i = 2:N_i;    
    border_a = floor(i/2)-open_ends;
    border_b = 2*i+open_ends;
    border_c = floor((i-N_i)/2)+N_j+open_ends;
    border_d = 2*i+N_j-2*N_i-open_ends;
    
    if i<lim_1 && i<lim_2
        k = max(2,border_a);
        l = min(c+1,border_b);
    elseif (i<lim_1 && i>lim_2)||(i<lim_2 && i>lim_1)
        k = max(2,border_a);
        l = min(N_j,border_c);
    else
        k = border_d;
        l = min(N_j,border_c);
    end
	D(i,2:k-1) = NaN;
    D(i,l+1:N_j) = NaN;
	for j = k:l
        % Scale the steps to discourage skipping ahead
        kk1 = 2; % long
        kk2 = 1; % diagonal
        kk3 = 5; % vertical and horizontal
        dd = D(i,j);
        [dmax, tb] = min([D(i-1, j-1)+dd*kk2, D(max(1,i-2), j-1)+dd*kk1,...
            D(i-1, max(1,j-2))+dd*kk1, D(i-1,j)+kk3*dd, D(i,j-1)+kk3*dd]);
        D(i,j) = dmax;
        phi(i,j) = tb;
	end
end

% Traceback from top left
i = r+1; 
j = c+1;
p = i;
q = j;
while i > 2 && j > 2
	tb = phi(i,j);
	if (tb == 1)
        i = i-1;
        j = j-1;
	elseif (tb == 2)
        i = i-2;
        j = j-1;
	elseif (tb == 3)
        j = j-2;
        i = i-1;
	elseif (tb == 4)
        i = i-1;
        %j = j;
	elseif (tb == 5)
        j = j-1;
        %i = i;
	else    
        error('dp error');
	end
    p = [i,p];
    q = [j,q];
end

% Strip off the edges of the D matrix before returning
D = D(2:(r+1),2:(c+1));

% map down p and q
p = p-1;
q = q-1;
