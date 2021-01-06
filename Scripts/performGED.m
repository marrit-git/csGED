function [ evecs, evals, covS, covR ] = performGED( Xs, Xr, varargin )
%PERFORMGED takes two matrices Xs and Xr (channels x time x trials). It 
% returns the results of a GED for covS > covR as well as the constructed 
% covariance matrices covS, covR. Returned matrices have no shrinkage 
% regularization ( because they're intended for topoplot construction), 
% but 1% shrinkage regularization is automatically applied to covR prior to 
% computing the GED.
% Function takes optional boolean "fast". If fast is true, covariance
% matrices are constructed from concatenated trials, instead of computing
% each individual covariance matrix per trial and averaging them. May
% sacrifice accuracy (i.e., assumes covariance stationarity across trials).
% Default is 0.

if isempty(varargin)
    fast = 0;
else
    fast = varargin{1};
end


if size(Xs,1) ~= size(Xr,1)
    error('GED signal inputs have different number of channels.');
end

[nbchan, pntsS, trialsS] = size(Xs);
[~, pntsR, trialsR] = size(Xr);

% Compute first covariance matrix
covS = zeros(nbchan, nbchan);
Xs = Xs - mean(Xs,2); % mean-center all trials
if fast == 1
    Xs = reshape(Xs, nbchan, pntsS * trialsS); % concatenate trials
    covS = Xs * Xs' / pntsS;
elseif fast == 0
    for i = 1:trialsS
        temp = Xs(:,:,i);
        covS = covS + temp*temp' / pntsS;
    end
end
covS = covS ./ trialsS;

% Compute second covariance matrix
covR = zeros(nbchan, nbchan);
Xr = Xr - mean(Xr,2); % mean-center all trials
if fast == 1
    Xr = reshape(Xr, nbchan, pntsR * trialsR); % concatenate trials
    covR = Xr * Xr' / pntsR;
elseif fast == 0
    for i = 1:trialsR
            temp = Xr(:,:,i);
            covR = covR + temp*temp' / pntsR;
    end
end
covR = covR ./ trialsR;

% Apply 1% shrinkage regularization to improve separability
g = 0.01;
covRr = (1-g)*covR + g*mean(eig(covR))*eye(nbchan);

% Apply GED
[evecs, evals] = eig(covS, covRr);
[evals, sidx] = sort(diag(evals), 'descend');
evecs = evecs(:, sidx);

% Normalize eigenvectors to unit length
for v = 1:size(evecs,2)
    evecs(:,v) = evecs(:,v)/norm(evecs(:,v));
end

end

