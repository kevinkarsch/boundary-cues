function Sigma = vec2PSD(X)

d = sqrt(1/4 + 2 * size(X,1)) - 0.5;
L_diag_log = X(1:d);
L_triu = X(d+1:end);

L_recon = diag(exp(L_diag_log));
L_recon(triu(true(size(L_recon)),1)) = L_triu;
Sigma = L_recon' * L_recon;
