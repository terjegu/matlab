%% Read files to LSF matrix
% Terje Gundersen 14.11.2009
close all;
clear all;

%% Read filenames
list = dir('data/source');

%% Calculate LPC features for both sounds
fs = 16e3;                      % Sampling frequency
window_size = 10e-3;            % 15ms frame length
len = floor(fs*window_size);	% samples per frame
anal = round(1*len);            % samples per analysis frame (overlapping)
p = 16;                         % LPC order (Fs/1000)
X_lsf = [];                     % Feature matrix used in training
N_iter = ceil(200/window_size/585); % 10ms
% N_iter = ceil(300/window_size/585) % 15ms
% N_iter = ceil(400/window_size/585) % 20ms
for i=4:(5+N_iter)
    filename = {list(i,1).name};
    x=wavread(['data/source/',filename{1}]);	% Read wav file
    X=lpcauto(x,p,[len anal]);                  % Make LPC matrix

    fn = length(X);
    X_lsf_temp = zeros(fn,p);
    for j=1:fn
        X_lsf_temp(j,:) = poly2lsf(X(j,:));     % Convert LPC to LSF
    end
    
    % Add to matrix
    X_lsf = [X_lsf;X_lsf_temp];
end

%%
save('wavfiles','X_lsf');