function out = visualizeDEM(Z, pad_percent, contrast)

if nargin < 3
  contrast = 0.75;
end

N = getNormals_conv(Z);

if nargin < 2
  pad_percent = 0;
end
range = prctile(Z(~isnan(Z)), [pad_percent, 100-pad_percent]);

Z = (Z - range(1))./max(eps,range(2) - range(1));
Z2 = min(1, max(0, Z));
Z2(isnan(Z)) = nan;
Z = Z2;
Z = mod(.75 - Z * .75, 1);
% Z = Z - min(Z(:));

S = N(:,:,3);
S = (S - min(S(:)))./max(eps,max(S(:)) - min(S(:)));
S = S*contrast + (1-contrast);



% keyboard

vis = max(0, min(1, hsv2rgb(cat(3, Z, ones(size(S))*.75, S))));



if nargout == 0
  imagesc(vis);
  imtight;
else
  out = vis;
end
