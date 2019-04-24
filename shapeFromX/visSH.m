function [out, N, E] = visSH(L, sz, DO_FLIP)

if nargin < 2
  sz = [256,256];
end

min_sz = min(sz);
j = ([1:sz(1)]-min_sz/2 - (sz(1)-min_sz)/2)*(2/min_sz);
i = ([1:sz(2)]-min_sz/2 - (sz(2)-min_sz)/2)*(2/min_sz);
[Y,X] = meshgrid(i,j);
Z = sqrt(max(0, 1-(X.^2 + Y.^2)));
valid = Z ~= 0;

if nargin >= 3 && DO_FLIP
  Z = -Z;
end

N = [X(valid), Y(valid), Z(valid)];
E = renderSH_helper(N, L);

% E = E - min(E(:));
% E = E ./ max(E(:));

% keyboard
% V = [0,0,1];
% R = repmat(2*(N*V'), [1,3]) .* N - repmat(V, size(N,1), 1);
% Es = renderSH_helper(R, L, NO_SHADOW)/2;

% m = max(E(:));
% E = E ./ m;
% E = E.^20;
% E = E * m;
% E = E + E.^50;

vis = nan(size(Z));
vis(valid) = E;% + Es.^30;

if nargout == 0
  imagesc(vis);
  imtight; colormap('gray');
else
  out = vis;
end


