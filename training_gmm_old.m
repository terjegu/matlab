%% TRAINING of GMM
% Terje Gundersen 30.10.2009
close all;
clear all;

%% Read filenames
list = dir('data/source');

%% Calculate LPC features for both sounds
fs = 16e3;                      % Sampling frequency
window_size = 10e-3;            % 10ms frame length
len = floor(fs*window_size);	% samples per frame
anal = round(1.5*len);          % samples per analysis frame (overlapping)
p = 16;                         % LPC order (Fs/1000)
m = 8;                         % Number of mixtures

for i=3:5
    filename = {list(i,1).name};
    filename{1}                                 % print filename
    x=wavread(['data/source/',filename{1}]);	% Read wav file
    X=lpcauto(x,p,[len anal]);                  % Make LPC matrix

    [fn,fl] = size(X);
    p = fl-1;
    X_lsf = zeros(fn,p);
    for j=1:fn
        X_lsf(j,:) = poly2lsf(X(j,:));          % Convert LPC to LSF
    end

    % Train GMM with EM-algorithm
    if exist('gm_obj')
        % If not first iteration
        S.mu = gm_obj.mu;
        S.Sigma = gm_obj.Sigma;
        S.PComponents = gm_obj.PComponents;
        gm_obj = gmdistribution.fit(X_lsf,m,'CovType','diagonal','Start',S,'Regularize',1e-4);
    else
        % First iteration
        gm_obj = gmdistribution.fit(X_lsf,m,'CovType','diagonal','Regularize',1e-4);
    end
end


%% Save variables
save('gmm','gm_obj','m');