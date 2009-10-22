function M = simmx_lpc(A,B)
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

end