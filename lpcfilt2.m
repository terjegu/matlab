function u=lpcfilt2(s,ar,t)
% u=lpcfilt2(s,ar,t)
%   LP filtering with memory
%   E, X and ar are matrices
%   Used in conversion_function.m to resynthesize converted signal

% Terje Gundersen 01.11.2009

[nf,p1]=size(ar);
dc=zeros(nf,1);

p=p1-1;
ns=length(s);

u=zeros(ns,1);
x0=(max(1,ceil(t(nf)-p)):ns)';

[u(x0),mem] = filter(1,ar(nf,:),s(x0)-dc(nf));
for i=nf-1:-1:2
	x0=(max(1,ceil(t(i)-p)):ceil(t(i+1)-1))';
	[u(x0),mem]=filter(1,ar(i,:),s(x0)-dc(i),mem);
end
x0=(1:ceil(t(2)-1))';
u(x0)=filter(1,ar(1,:),s(x0)-dc(1),mem);

end

