function [V,Gamma] = param(k,n,m,P,X_lsf,Y_lsf,gm_obj,sigma_diag)

% Calculate the matrix D = P(C|x) * (x-mu)^T * Sigma^-1
D = zeros(n,m);
for l=1:n
	D(l,:) = (X_lsf(l,k)-gm_obj.mu(:,k)).*sigma_diag(:,k);
end
D = D.*P;

% Conversion variables
param_vg = ([P';D']*[P D])\[P';D']*Y_lsf(:,k);
V = param_vg(1:m);
Gamma = param_vg((m+1):end);

end