function x=lpcfilt2(e,ar,t)
% x=lpcfilt2(e,ar,t)
%   LP filtering with memory
%   e is a vector
%   t and ar are matrices
%   Used in conversion_function.m to resynthesize converted signal

% Terje Gundersen 01.11.2009

nf = length(ar);
ne = length(e);
x = zeros(ne,1);
mem = zeros(16,1);

start = 1;
endp = round(t(1,1)+0.5*t(2,1));
for i=1:nf-2
    [x(start:endp),mem] = filter(1,ar(i,:),e(start:endp),mem); 
    start = endp+1;
    endp = min(round(start+0.5*(t(i+1,1)+t(i+2,1))-1),ne);
end
x(start:endp) = filter(1,ar(nf-1,:),e(start:endp),mem);

end