function X = split(x,fl)
% X = split(x,fl)
%   Convert vector to matrix
N = length(x);
fn = floor(N/fl);
X = zeros(fn,fl);

for i=1:fn
    X(i,:) = x(1+fl*(i-1):fl*i);
end

end
