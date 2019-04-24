function im = reconSpyr_simple(pyr2, scale)

pind = cellfun(@(x) size(x), pyr2, 'UniformOutput', false);
pind = cat(1, pind{:});

pyr = zeros(sum(prod(pind,2)),1);
for s = 1:size(pind,1)
  band = pyr2{s};
  mult = scale.^(log2(pind(1,1) / pind(s,1)));
  pyr(pyrBandIndices(pind, s)) = band(:) / mult;
end

im = reconSFpyr(pyr, pind);
