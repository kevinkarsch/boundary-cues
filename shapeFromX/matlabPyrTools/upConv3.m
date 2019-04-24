% RES = upConv(IM, FILT, SZ)

function res = upConv3(vol,filt,sz)

tmp = zeros(sz);
tmp(1:2:end, 1:2:end, 1:2:end) = vol;

res = permute(convn(permute(convn(permute(convn(tmp, filt, 'same'), [3,1,2]), filt, 'same'), [3,1,2]), filt, 'same'),[3,1,2]);
