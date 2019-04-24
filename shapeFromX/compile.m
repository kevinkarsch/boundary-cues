
curdir = pwd;

cd minFunc_2012/mex/
try
mex lbfgsAddC.c
mex lbfgsC.c
mex lbfgsProdC.c
mex mcholC.c
catch
  fprintf('minFunc compile failed\n');
end
cd(curdir)

try
  convnfft_install
catch
  fprintf('convnfft installation failed\n');
end

try
  mex getK_fast.c
catch
  fprintf('getK_fast compile failed\n');
end

try
  mex getK_fast_down.c
catch
  fprintf('getK_fast_down compile failed\n');
end


try
  mex splat3_fast.c
catch
  fprintf('splat3_fast compile failed\n');
end

try
  mex splat3_backprop_fast.c
catch
  fprintf('splat3_backprop_fast compile failed\n');
end



try
  mex hullDist.c
catch
  fprintf('hullDist compile failed\n');
end


try
  mex splineHist_fast.c;
catch
  fprintf('splineHist_fast compile failed\n');
end

try
  mex accumarray_fast.c;
catch
  fprintf('accumarray_fast compile failed\n');
end


try
  mex interp2_fixed_fast.c;
catch
  fprintf('interp2_fixed_fast compile failed\n');
end

try
  mex interp2_fixed_sum_fast.c;
catch
  fprintf('interp2_fixed_sum_fast compile failed\n');
end

try
  mex interp2_fixed_sum_weighted.c;
catch
  fprintf('interp2_fixed_sum_fast compile failed\n');
end

try
  mex interp1_fixed_sum_weighted.c;
catch
  fprintf('interp1_fixed_sum_weighted compile failed\n');
end

try
  mex interp1_fixed_sum_fast.c;
catch
  fprintf('interp1_fixed_sum_fast compile failed\n');
end

try
  mex interp1_fixed_fast.c;
catch
  fprintf('interp1_fixed_fast compile failed\n');
end


try
  mex renderSH_helper_mex.c
catch
  fprintf('renderSH_helper_mex compile failed\n');
end

% try
%   mex getNormals_fast.c;
% catch
%   fprintf('getNormals_fast compile failed\n');
% end

try
  mex mex_imresize.c;
catch
  fprintf('mex_imresize compile failed\n');
end

try
  D = pwd;
  cd binaryops
  mex bdil6.c;
  mex bero4.c
  cd(D)
catch
  fprintf('binaryops compile failed\n');
  cd(D)
end

