function E = lpcfilt(X,ar)
% E = lpcfilt(X,ar)
% LP filtering with memory
% E, X and ar are matrices

% Terje Gundersen 01.11.2009

[fn,fl] = size(X);
E = zeros(fn,fl);

[E(1,:),mem] = filter(ar(1,:),1,X(1,:));
for i=2:fn
	[E(i,:),mem] = filter(ar(i,:),1,X(i,:),mem);
end

end

