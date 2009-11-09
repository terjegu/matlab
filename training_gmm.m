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
m = 32;                          % Number of mixtures
X_lsf = [];                     % Feature matrix used in training

for i=3:50
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

%% Train GMM with EM-algorithm and kmeans for initialisation 
% VQ for initialisation
[S.mu,~,J]=kmeans(X_lsf,m);

% Variance of each cluster
S.Sigma = zeros(1,p,m);
S.PComponents = zeros(1,m);
for i=1:m
    S.Sigma(1,:,i) = max(1e-5,var(X_lsf(J==i)));
    S.PComponents(1,i) = sum(J==i)/length(J);
end

% GMM with EM
gm_obj = gmdistribution.fit(X_lsf,m,'CovType','diagonal','Start',S,'Regularize',1e-5);


%% Save variables
save('gmm','gm_obj','m');