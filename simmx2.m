function M = simmx2(A,B)
% M = simmx2(A,B) calculate a sim matrix 
% between MFCC feature matrices A and B.

[m,na] = size(A); % m segments of length na
nb = size(B,2); % m segments of length nb

A_cep = zeros(m,na);
B_cep = zeros(m,nb);

for i=1:m
    A_cep(i,:) = ifft(log(A(i,:)));
    B_cep(i,:) = ifft(log(B(i,:)));
end

M = zeros(na,nb);
for i = 1:na
	for j = 1:nb
        temp = A_cep(i,1:20)-B_cep(1:20,j);
        M(i,j) = norm(temp);
	end
end
M = M/m;

end


% EA = sqrt(sum(A.^2));
% EB = sqrt(sum(B.^2));
% denom = EA'*EB

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







