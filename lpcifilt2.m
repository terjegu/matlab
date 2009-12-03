function e=lpcifilt2(s,ar,t)
% e=lpcifilt2(s,ar,t)
%   LP filtering with memory
%   s is a vector
%   t and ar are matrices
%   Used in conversion_function.m to filter signal

% Terje Gundersen 01.11.2009

nf = length(ar);
ns = length(s);
e = zeros(ns,1);
mem = zeros(16,1);

start = 1;
endp = round(t(1,1)+0.5*t(2,1));
for i=1:nf-2
    [e(start:endp),mem] = filter(ar(i,:),1,s(start:endp),mem);
    start = endp+1;
    endp = min(round(start+0.5*(t(i+1,1)+t(i+2,1))-1),ns);
%     pm_t = round(start+(endp-start)/2)
%     delta = endp-start
end
e(start:endp) = filter(ar(nf-1,:),1,s(start:endp),mem);

end