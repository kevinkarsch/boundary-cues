function pyr = buildGpyr(im, depth, filt, edges)

[pyr2, pind] = buildLpyr(im, depth, filt, filt, edges);

pyr = {};
for s = 1:size(pind,1)
  pyr{s} = pyrBand(pyr2, pind, s);
end
