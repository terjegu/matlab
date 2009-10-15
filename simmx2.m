% function M = simmx2(A,B)
% M = simmx2(A,B)
%    calculate a sim matrix between MFCC feature matrices A and B.



%% SIMMX2
close all;
clear all;

%% Load two speech waveforms of the same utterance (from TIMIT)
[d1,sr] = wavread('data/t01s000228.wav');
% [d2,sr2] = wavread('data/t03s000228.wav');

test = rceps(d1);
% back = icceps(test);
plot(test(1:20))




% wavwrite(back,sr,'data/test_ccep.wav')

% EA = sqrt(sum(A.^2));
% EB = sqrt(sum(B.^2));
% 
% M = (A'*B)./(EA'*EB);