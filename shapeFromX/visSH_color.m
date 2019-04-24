function out = visSH(L, sz, color_enhance)

if nargin < 2
  sz = [256,256];
end

if nargin < 3
  color_enhance = 1;
end

if nargin < 3
  NO_SHADOW = 1;
end

min_sz = min(sz);
j = ([1:sz(1)]-min_sz/2 - (sz(1)-min_sz)/2)*(2/min_sz);
i = ([1:sz(2)]-min_sz/2 - (sz(2)-min_sz)/2)*(2/min_sz);
[Y,X] = meshgrid(i,j);
Z = sqrt(max(0, 1-(X.^2 + Y.^2)));
valid = Z ~= 0;

V = {};
for c = 1:3
  E = renderSH_helper([X(valid), Y(valid), Z(valid)], L(:,c));
  v = nan(size(Z));
  v(valid) = E;
  V{c} = v;
end
V = cat(3, V{:});

V_avg = mean(V,3);
V = repmat(V_avg, [1,1,3]) + color_enhance * (V - repmat(V_avg, [1,1,3]));

V = exp(V);


% V = V ./ max(V(:));

if nargout == 0
  imagesc(min(1,V), [0,1]);
  imtight; colormap('gray');
else
  out = V;
end

