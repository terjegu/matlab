function y = conversion_function(x)

y = sum(P(:,x) * (v_i + Gamma_i*sigma_i(x-mu_i)));

end