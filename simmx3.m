function M = simmx3(A,B)
% M = simmx2(A,B) calculate a sim matrix 
% between MFCC feature matrices A and B.

na = length(A); % m segments of length na
nb = length(B); % m segments of length nb

M = zeros(na,nb);
for i = 1:na
	for j = 1:nb
        M(i,j) = norm(A(i,:)-B(j,:));
	end
end
M = M./(max(max(M))+1);

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







