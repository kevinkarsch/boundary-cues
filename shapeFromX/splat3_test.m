
clear all;

mex splat3_fast.c
mex splat3_backprop_fast.c

% for ii = 1:100
  
  
  X = randn(100,3);
  bin_low = [-2, -3, -1];
  bin_high = [2, 4, 1.5];
  sigma = 1;
  
  for ii = 1
    X_splat = splat3(X, sigma, bin_low, bin_high);
    % imagesc(visStack_mean(X_splat.N)); imtight; drawnow;
    
    
    
    X_splat2 = splat3_fast(X, sigma, bin_low, bin_high);
  end
  
% end
  
%   X_splat
%   X_splat2
%   
%   names = fieldnames(X_splat2);
%   for name_i = 1:length(names)
%     name = names{name_i};
%     max(abs(X_splat2.(name)(:) - X_splat.(name)(:)))
%     %   sum(X_splat2.(name)(:) ~= X_splat.(name)(:))
%   end
% 
  
M = randn(size(X_splat.N));
D1 = splat3_backprop(M, X_splat);
  
D2 = splat3_backprop_fast_wrapper(M, X_splat);
% D2 = splat3_backprop_fast(M/X_splat.bin_width, X_splat.idx111, X_splat.idx112, X_splat.idx121, X_splat.idx122, X_splat.idx211, X_splat.idx212, X_splat.idx221, X_splat.idx222, X_splat.f1, X_splat.f2, X_splat.valid);


max(max(abs(D1 - D2)))
