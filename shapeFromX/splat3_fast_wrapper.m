function X_splat = splat3_fast_wrapper(X, bin_range_low, bin_range_high)

% try
  
  X_splat = splat3_fast(X, bin_range_low, bin_range_high);
  X_splat.bin_range_low = bin_range_low;
  X_splat.bin_range_high = bin_range_high;
  
% catch
%   
%   fprintf('Fast splatting failed\n');
%   X_splat = splat3(X, bin_range_low, bin_range_high);
%   
% end