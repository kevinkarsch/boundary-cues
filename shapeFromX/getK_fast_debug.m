
clear all
close all;

mex getK_fast.c

Z = rand(1000,500);

dL_K = rand(size(Z));

tic;
for ii = 1:10
  [K, dK] = getK(Z, 3);
  dL_Z = getK_backprop( dL_K, dK, 3);
end
toc

tic;
for ii = 1:10
  [K_fast, dK_fast] = getK_fast(Z);
  dL_Z_fast = getK_backprop_fast_mat(dL_K, dK_fast);
end
toc

max(max(abs(K - K_fast)))
max(max(abs(dK.denom(:) - dK_fast(1,:)')))
max(max(abs(dK.f1(:) - dK_fast(2,:)')))
max(max(abs(dK.f2(:) - dK_fast(3,:)')))
max(max(abs(dK.f11(:) - dK_fast(4,:)')))
max(max(abs(dK.f22(:) - dK_fast(5,:)')))
max(max(abs(dK.f12(:) - dK_fast(6,:)')))
max(max(abs(dL_Z_fast - dL_Z)))

