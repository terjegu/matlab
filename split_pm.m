function X = split_pm(x,fl,anal,skip)
% X = split(x,fl)
%   Convert vector to matrix
N = length(x);
fn = length(anal);
X = zeros(fn,max(anal)-1);
for i=1:fn
    start = max(1,ceil(fl*(i-1)+skip(i)));     % start (len + skip)
    stop = ceil(fl*(i-1)+skip(i)+anal(i)-1);   % stop (len + anal + skip)
    X(i,1:stop-start+1) = x(start:stop);
end

end
