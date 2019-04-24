function [loglike, d_loglike] = GSM_pdf(model, data, USE_LUT)

RETURN_GRADIENT = (nargout >= 2);

if nargin < 3
  USE_LUT = true;
end

K = length(model.pis);
M = size(data,2);
N = size(data,1);

assert(M == 1)

if USE_LUT && isfield(model, 'lut')
  if RETURN_GRADIENT
    try
      [loglike, d_loglike] = interp1_fixed_sum_fast(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F_LL);
    catch
      fprintf('executing interp1_fixed_sum_fast failed\n');
      [loglike, d_loglike] = interp1_fixed_sum(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F_LL);
    end
  else
    try
      [loglike] = interp1_fixed_sum_fast(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F_LL);
    catch
      fprintf('executing interp1_fixed_sum_fast failed\n');
      [loglike] = interp1_fixed_sum(data, model.lut.bin_range(1), model.lut.bin_width, model.lut.F_LL);
    end
  end
  
else

  P = 0;
  for k = 1:K
    P = P + exp(log(model.pis(k)) + lnormpdf(data, model.mu, model.sigs(k)));%model.pis(k) * normpdf(data, model.mu, model.sigs(k));
  end

  loglike = log(max(exp(-300),P));
  
end
