function [V, M] = renyi3_fixed_sub(X, X_splat, sigma, USE_LOG)

K = 3;

X_count = size(X,1);

if USE_LOG == 2
  Z = X_count^(2/K);
else
  Z = X_count^(2/K) * sqrt(4*pi*sigma^2);
end

N_SIGMAS = 6;
half_width = min(X_splat.n_bins, ceil(N_SIGMAS * sigma / X_splat.bin_width));
blur_kernel = exp(-([-half_width:half_width]'*X_splat.bin_width).^2 / (4*sigma^2))/Z;

N1 = squeeze(any(any(X_splat.N, 2),3));
N2 = squeeze(any(any(X_splat.N, 1),3))';
N3 = squeeze(any(any(X_splat.N, 1),2));
h = half_width;
N1idx = max(1, find(N1, 1, 'first') - h) : min(size(X_splat.N,1), find(N1, 1, 'last') + h);
N2idx = max(1, find(N2, 1, 'first') - h) : min(size(X_splat.N,2), find(N2, 1, 'last') + h);
N3idx = max(1, find(N3, 1, 'first') - h) : min(size(X_splat.N,3), find(N3, 1, 'last') + h);

N_sub = X_splat.N(N1idx, N2idx, N3idx);
M = zeros(size(X_splat.N));
M(N1idx, N2idx, N3idx) = permute(convn(permute(convn(permute(convn(N_sub, blur_kernel, 'same'), [3,1,2]), blur_kernel, 'same'), [3,1,2]), blur_kernel, 'same'),[3,1,2]);

% M = permute(convn(permute(convn(permute(convn(X_splat.N, blur_kernel, 'same'), [3,1,2]), blur_kernel, 'same'), [3,1,2]), blur_kernel, 'same'),[3,1,2]);

if nargout == 0
  imagesc([mean(M,3), squeeze(mean(M,2)); squeeze(mean(M,1))', nan(size(M,3), size(M,3))].^(1/4)); imtight; drawnow;
end

M = max(eps, M);

V = sum(sum(sum(M .* X_splat.N)));

if USE_LOG
  M = -M ./ V;
  V = -log(V);
end

M = 2*M;
