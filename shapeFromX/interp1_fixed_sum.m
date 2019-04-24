function [F, dF] = interp1_fixed(X, low_val, bin_width, LUT_F, LUT_dF)

QUADRATIC_EXTRAPOLATION = false;

X = (X-low_val) / bin_width;

above_idx = (X >= (length(LUT_F)-1));
below_idx = (X <= 0);
within_idx = ~(below_idx | above_idx);

below_bin = floor(X(within_idx));

V1 = X(within_idx) - below_bin;
V2 = 1 - V1;

below_bin = below_bin+1;
above_bin = below_bin+1;

if size(LUT_F,1) > size(LUT_F,2)
  LUT_F = LUT_F';
end
F = LUT_F(below_bin) * V2 + LUT_F(above_bin) * V1;


if nargout >= 2
  dF = zeros(size(X));
  dF(within_idx) = (LUT_F(above_bin) - LUT_F(below_bin))/(bin_width);
%   dF(within_idx) = LUT_dF(below_bin) .* V2 + LUT_dF(above_bin) .* V1;  
end


above_F = LUT_F(end)';
above_dF = LUT_dF(end)';
above_ddF = (above_dF - LUT_dF(end-1));


delta = X(above_idx) - (length(LUT_F) - 1);
if QUADRATIC_EXTRAPOLATION
  F = F + length(delta)*above_F + (above_dF*bin_width)*sum(delta) + (0.5*above_ddF*bin_width)*sum(delta.^2);
else
  F = F + length(delta)*above_F + (above_dF*bin_width)*sum(delta);
end

if nargout >= 2
  if QUADRATIC_EXTRAPOLATION
    dF(above_idx) = above_dF + (above_ddF)*delta;
  else
    dF(above_idx) = above_dF;
  end
end

below_F = LUT_F(1);
below_dF = LUT_dF(1);
below_ddF = (LUT_dF(2) - below_dF);


delta = X(below_idx);
if QUADRATIC_EXTRAPOLATION
  F = F + length(delta)*below_F + (below_dF*bin_width)*sum(delta) + (0.5*below_ddF*bin_width)*sum(delta.^2);
else
  F = F + length(delta)*below_F + (below_dF*bin_width)*sum(delta);
end

if nargout >= 2
  if QUADRATIC_EXTRAPOLATION
    dF(below_idx) = below_dF + (below_ddF)*delta;
  else
    dF(below_idx) = below_dF;
  end
end
