function [loglike, d_loglike, ps] = MOG_mvn_pdf(data, model, USE_COST)

if nargin < 3
  USE_COST = 0;
end

K = length(model.pis);
D = size(data,2);
N = size(data,1);


P = 0;
if nargout >= 2
  d_loglike = zeros(N,D);
end

if ~isempty(data)
  for k = 1:K

    for d = 1:D
      data_diff(:,d) = data(:,d) - model.mus(k,d);
    end
    
%     tmp = data_diff * model.iSigmas(:,:,k);
%     p = exp(-0.5 * sum(tmp .* data_diff,2) + model.logZ(k));
    
    p = model.pis(k) * exp(lmvnpdf(data, model.mus(k,:), model.Sigmas(:,:,k)));
    tmp = data_diff * inv(model.Sigmas(:,:,k));
    
    P = P + p;
    
    if nargout >= 2
      for d = 1:D
        d_loglike(:,d) = d_loglike(:,d) + -p .* tmp(:,d);
      end
    end
    
    if nargout >= 3
      ps{k} = p;
    end

  end
end

% P = max(eps, P);
loglike = log(P);

if nargout >= 2
  for d = 1:D
    d_loglike(:,d) = d_loglike(:,d) ./ P;
  end
end

if USE_COST
  loglike = model.LL_max-loglike;
  if nargout >= 2
    d_loglike = -d_loglike;
  end
end
