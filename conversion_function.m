function y = conversion_function(X_rc,P,sigma,mu,v,Gamma)

[n,m] = size(P);
X_rc = zeros(size(X_rc));

for i=1:n
    for j=1:m
        y() = sum(P(i,:) * (v(i,:) + Gamma(i,:)*sigma(i)(x-mu_i)));
    end
end

end