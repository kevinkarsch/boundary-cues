% RES = corrDn3(IM, FILT)
function res = corrDn3(vol, filt)

tmp = permute(convn(permute(convn(permute(convn(vol, filt, 'same'), [3,1,2]), filt, 'same'), [3,1,2]), filt, 'same'),[3,1,2]);
res = tmp(1:2:end,1:2:end,1:2:end);

