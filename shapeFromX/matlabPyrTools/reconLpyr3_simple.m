function res = reconLpyr3_simple(pyr, filt)

res = pyr{end};
for d = (length(pyr)-1):-1:1
  
  res = upConv3(res, filt, size(pyr{d}));
  
  res = res + pyr{d};  
end

