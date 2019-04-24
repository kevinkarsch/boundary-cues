function d_loss_X = splat3_backprop_fast_wrapper(d_loss_N, X_splat)

try
  
  d_loss_X = splat3_backprop_fast(d_loss_N/X_splat.bin_width, X_splat.idx111, X_splat.idx112, X_splat.idx121, X_splat.idx122, X_splat.idx211, X_splat.idx212, X_splat.idx221, X_splat.idx222, X_splat.f1, X_splat.f2, X_splat.valid);
  
catch
  
  fprintf('Fast splatting failed\n');
  d_loss_X = splat3_backprop(d_loss_N, X_splat)
  
end