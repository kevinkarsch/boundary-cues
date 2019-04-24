function [loglike, d_loglike, LLs, mahal_dist, d_loglike_mult] = GSM_pdf(model, data, USE_LUT)

RETURN_GRADIENT = (nargout >= 2);

if nargin < 3
  USE_LUT = true;
end

K = length(model.pis);
M = size(data,2);
N = size(data,1);

if USE_LUT && isfield(model, 'lut')
  
  if M == 1
    
    if RETURN_GRADIENT
      try
        [loglike, d_loglike] = interp1_fixed_sum_fast(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      catch
        fprintf('executing interp1_fixed_sum_fast failed\n');
        [loglike, d_loglike] = interp1_fixed_sum(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      end
    else
      try
        [loglike] = interp1_fixed_sum_fast(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      catch
        fprintf('executing interp1_fixed_sum_fast failed\n');
        [loglike] = interp1_fixed_sum(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      end
    end
    
  else
      
    if isfield(model, 'mu')
      for j = 1:M
        data(:,j) = data(:,j) - model.mu(j);
      end
    end
    
    if RETURN_GRADIENT
      
      data_icov = data * model.Sigma_inv;
      mahal_dist = 0.5 * sum(data_icov .* data,2);
      
%       epsilon = 10^-2;
%       invalid = mahal_dist < epsilon;
%       mahal_dist = max(mahal_dist, epsilon);
      
      try
        [loglike, d_loglike_mult] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      catch
        fprintf('executing interp1_fixed_sum_fast failed\n');
        [loglike, d_loglike_mult] = interp1_fixed_sum(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      end
      
%       d_loglike_mult(invalid) = 0;
      
      d_loglike = bsxfun(@times, d_loglike_mult, data_icov);
      
    else

      mahal_dist = 0.5 * sum((data * model.Sigma_inv) .* data,2);
      
      try
        [loglike] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      catch
        fprintf('executing interp1_fixed_sum_fast failed\n');
        [loglike] = interp1_fixed_sum(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, model.lut.F);
      end
      
    end
  end
  
  if USE_LUT == 2
    loglike = N*model.LL_zero - loglike;
    if RETURN_GRADIENT
      d_loglike = -d_loglike;
    end
  end
  
  
  LLs = [];

%   if M == 1
%     val = (data-model.lut.bin_range(1)) / model.lut.bin_width;
%   else
%     data_icov = data * model.Sigma_inv;
%     mahal_dist = 0.5 * sum(data_icov .* data,2);
%     val = (mahal_dist-model.lut.bin_range(1)) / model.lut.bin_width;
%   end
% 
%   val = max(0, min(length(model.lut.F)-2, val));
%   below_bin = floor(val);
% 
%   val = val - below_bin;
%   val2 = 1 - val;
%   
%   below_bin = below_bin+1;
%   above_bin = below_bin+1;
% 
%   loglike = model.lut.F(below_bin) * val2 + model.lut.F(above_bin) * val;
% 
%   if M == 1
% 
%     if nargout >= 2
%       d_loglike = model.lut.dF(below_bin) .* val2 + model.lut.dF(above_bin) .* val;
%     end
% 
%   else
% 
%     if nargout >= 2
%       d_loglike_mult = model.lut.dF(below_bin) .* val2 + model.lut.dF(above_bin) .* val;
%       d_loglike = zeros(N,M);
%       for m = 1:M
%         d_loglike(:,m) = d_loglike_mult .* data_icov(:,m);
%       end
%     end
% 
%   end

  
else

  if ~isfield(model, 'Sigma_inv')
    model.Sigma_inv = inv(model.Sigma);
  end

  if ~isfield(model, 'Sigma_R')
    model.Sigma_R = cholcov(model.Sigma,0);
  end

  if isfield(model, 'mu')
    for j = 1:M
      data(:,j) = data(:,j) - model.mu(j);
    end
  end

  data_icov = data * model.Sigma_inv;
  mahal_dist = 0.5 * sum(data_icov .* data,2);

  P = 0;
  p_div = 0;
  if nargout >= 3
    LLs = zeros(N,K);
  end

  for k = 1:K

    if isfield(model, 'logmults')
      logmult = model.logmults(k);
    else
      log_len = log((2*pi).^(M/2));

      %     [R] = cholcov(model.Sigma*model.vars(k),0);
      R = model.Sigma_R * sqrt(model.vars(k));

      logmult = log(model.pis(k)) - log_len - sum(log(diag(R)));
    end

    %     log(sqrt(det(model.Sigma*model.vars(k))))
    %     mult = (model.pis(k) ./ ((2*pi).^(size(data,2)/2) * sqrt(det(model.Sigma*model.vars(k)))));

    if nargout >= 3
      l = logmult + mahal_dist./-model.vars(k);
      LLs(:,k) = l;
      p = exp(l);
      P = P + p;
    else
      p = exp(logmult + mahal_dist./-model.vars(k));
      P = P + p;
    end

    if nargout >= 2
      p_div = p_div + p./-model.vars(k);
    end
  end
  
%   P = max(eps, P);
  loglike = log(P);
  loglike(loglike < -1000000) = -1000000;

  if RETURN_GRADIENT
    p_div = p_div ./ P;
    d_loglike = zeros(N,M);
    for m = 1:M
      d_loglike(:,m) = p_div .* data_icov(:,m);
    end

    %   d_loglike = repmat(p_div ./ P, 1, M) .* data_icov;
  end
end
