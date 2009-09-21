function makelpc(filename)
% [x,f_s,WMODE,FIDX]=readwav('data/000358_JF.wav');

window_size = 20e-3; % 20ms
inc = f_s * window_size; % number of frames

F = enframe(x,inc); 
% x_lp = lpcauto(X,12);

end
