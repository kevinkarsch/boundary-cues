function pyr = buildGpyr(im, depth, filt, edges)

pyr = {im};
for d = 2:depth
  
  im_sz = size(pyr{d-1});
  
  if (im_sz(2) == 1)
    pyr{d} = corrDn(pyr{d-1}, filt, edges, [2 1], [1 1]);
  elseif (im_sz(1) == 1)
    pyr{d} = corrDn(pyr{d-1}, filt', edges, [1 2], [1 1]);
  else
    lo     = corrDn(pyr{d-1}, filt', edges, [1 2], [1 1]);
    pyr{d} = corrDn(lo, filt, edges, [2 1], [1 1]);
  end
  
end
  
