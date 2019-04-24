function pyr = buildGpyr(vol, depth, filt)

pyr = {vol};
for d = 2:depth
  
  pyr{d} = corrDn3(pyr{d-1}, filt);
  
end
  
