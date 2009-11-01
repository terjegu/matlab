function X = split(x,fl,overlap)

fn = floor(length(x)/fl);
X = zeros(fn,fl+2*overlap);

X(1,overlap+1:end) = x(1:fl+overlap);
for i=2:fn-1
    X(i,:) = x(1+fl*(i-1)-overlap:fl*i+overlap);
end

X(fn,1:fl+overlap) = x(1+fl*(fn-1)-overlap:fl*fn);

hw = (hamming(fl+2*overlap))';
for i = 1:fn
    X(i,:) = X(i,:).*hw;
end

end