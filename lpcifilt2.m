function X = lpcifilt2(E,ar)
% X = lpcifilt2(E,ar)
% Inverse LP filtering with memory
% X, E and ar are matrices

% Terje Gundersen 01.11.2009

[fn,fl] = size(E);
X = zeros(fn,fl);

[X(1,:),mem] = filter(1,ar(1,:),E(1,:));
for i=2:fn
	[X(i,:),mem] = filter(1,ar(i,:),E(i,:),mem);
end

end

