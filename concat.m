function y = concat(X,fl)
% y = concat(X,fl)
%   fl is the new frame length used in the resulting vector y

[fn,len] = size(X);
Y = zeros(fn,fl);

for i=1:fn
    Y(i,:) = X(i:len);
end

Y = Y';
y = Y(:);
end