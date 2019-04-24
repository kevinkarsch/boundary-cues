function [X, dX] = PSD2vec(Sigma)

L = chol(Sigma);
% L = cholcov(Sigma);

L_diag_log = log(diag(L));
% L(tril(true(size(L)))) = nan;
% L_triu = L(~isnan(L));
L_triu = L(triu(true(size(L)),1));

X = [L_diag_log; L_triu];

if nargout >= 2
  
  dX = {};
  step = 10^-5;
  for ii = find(triu(ones(size(Sigma))))'
    
    j = ceil(ii / size(Sigma,1));
    i = ii - (j-1) * size(Sigma,1);
%     [i2,j2] = ind2sub(size(Sigma), ii);
%     assert( (i2 == i) && (j2 == j) )
    
    Sigma2 = Sigma;
    Sigma2(i,j) = Sigma(i,j) + step;
    Sigma2(j,i) = Sigma2(i,j);
    X2 = PSD2vec(Sigma2);
    dX{ii} = (X2 - X) / step;
    
    jj = (i-1)*size(Sigma,1) + j;
%     jj2 = sub2ind(size(Sigma), j, i);
%     assert(jj == jj2)
    
    dX{jj} = dX{ii};
  end
  dX = cat(2,dX{:});
  
  assert(size(dX,1) == size(X,1))
  assert(size(dX,2) == numel(Sigma))
  
end