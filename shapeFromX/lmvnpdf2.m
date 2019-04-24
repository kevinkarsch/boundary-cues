function [LL, dLL] = lmvnpdf2(X, mu, Sigma, LOSS_MODE)
    
if nargin < 4
  LOSS_MODE = 0;
end

if all(mu == 0)
  X_norm = X;
else
  X_norm = X - repmat(mu, size(X,1), 1);
end

X_norm_iSigma = X_norm * inv(Sigma);

if LOSS_MODE
  LL = 0.5 * sum((X_norm_iSigma) .* X_norm,2);
  dLL = X_norm_iSigma;
else
  logZ = 0.5*(log((2*pi).^size(X,2)) + log(det(Sigma)));
  
  LL = -logZ - 0.5 * sum((X_norm_iSigma) .* X_norm,2);
  dLL = -X_norm_iSigma;
end