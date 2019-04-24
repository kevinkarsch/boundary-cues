function pyr2 = buildSpyr_simple(im, depth, scale)

[pyr, pind] = buildSFpyr(im, depth);

pyr2 = {};
for s = 1:size(pind,1)
  band = pyrBand(pyr, pind, s);
  mult = scale.^(log2(pind(1,1) / pind(s,1)));
  pyr2{s} = mult * band;
end

