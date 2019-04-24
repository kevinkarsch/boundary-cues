function pyr = buildGpyr_color(im, depth, filt, edges)

if nargin < 3
  filt = [1;4;6;4;1]/16;
end

if nargin < 4
  edges = 'repeat';
end

pyr = {im};
for d = 2:depth
  
  im_sz = size(pyr{d-1});
  
  pyr{d} = {};
  for c = 1:3
    
    if (im_sz(2) == 1)
      pyr{d}{c} = corrDn(pyr{d-1}(:,:,c), filt, edges, [2 1], [1 1]);
    elseif (im_sz(1) == 1)
      pyr{d}{c} = corrDn(pyr{d-1}(:,:,c), filt', edges, [1 2], [1 1]);
    else
      lo     = corrDn(pyr{d-1}(:,:,c), filt', edges, [1 2], [1 1]);
      pyr{d}{c} = corrDn(lo, filt, edges, [2 1], [1 1]);
    end
    
  end
  
  pyr{d} = cat(3, pyr{d}{:});
  
end
  
