function y = concat_pm(X,pm,len)

[fn,fl] = size(X); 
y = zeros(1,len+1);

y(1:pm(1)) = X(1,1:pm(1));
for i=2:fn
    lastframe = X(i-1,(pm(i)-pm(i-1)+1):(pm(i+1)-pm(i-1)));
    thisframe = X(i,1:(pm(i)-pm(i-1)));
    test = lastframe+thisframe;
    sum = size(test)
    target = size(pm(i-1)+1:pm(i))
    y(1,pm(i-1)+1:pm(i)) = lastframe+thisframe;
end
% y(pm(length(pm)+1:end) = X(1,1:pm(1));
 

% [fn,len] = size(X);
% Y = zeros(fn,fl);
% 
% Y(1,:) = X(1:fl);
% Y(1,fl-overlap+1:fl) = Y(1,fl-overlap+1:fl) + X(2,1:overlap);
% for i=2:fn-1
%     Y(i,:) = X(i,overlap+1:len-overlap);
%     Y(i,1:overlap) = Y(i,1:overlap) + X(i-1,len-overlap+1:len);
%     Y(i,fl-overlap+1:fl) = Y(i,fl-overlap+1:fl) + X(i+1,1:overlap);
% end
% Y(fn,:) = X(fn,overlap+1:len-overlap);
% Y(fn,1:overlap) = Y(fn,1:overlap) + X(fn-1,len-overlap+1:len);
% 
% Y = Y';
% y = Y(:);
end