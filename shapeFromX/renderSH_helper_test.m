clear all;

mex renderSH_helper_mex.c

NO_SHADOW = false;
L = sphereLight(1,1);
sz = [16,16];

min_sz = min(sz);
j = ([1:sz(1)]-min_sz/2 - (sz(1)-min_sz)/2)*(2/min_sz);
i = ([1:sz(2)]-min_sz/2 - (sz(2)-min_sz)/2)*(2/min_sz);
[Y,X] = meshgrid(i,j);
Z = sqrt(max(0, 1-(X.^2 + Y.^2)));
valid = Z ~= 0;
N = [X(valid), Y(valid), Z(valid)];

tic;
for i = 1:10000
  [E1, dL1, dN1]= renderSH_helper_mat(N, L);
end
toc


tic;
for i = 1:10000
  [E2, dL2, dN2] = renderSH_helper_mex(N, L);
end
toc

err = max(abs(E1 - E2))
assert(err < 10^-14)

err = max(abs(dL2(:) - dL1(:)))
assert(err < 10^-14)

err = max(abs(dN1(:) - dN2(:)))
assert(err < 10^-14)
