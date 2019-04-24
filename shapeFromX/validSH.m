function valid = validSH(L)

L = mean(L,2);

Vfront = visSH(L, [64,64], 0);
Vback = visSH(L, [64,64], 1);

valid = max(Vfront(:)) >= max(Vback(:));

% C = conv2(Vfront, laplacian3, 'valid');
% valid = valid & mean(C(~isnan(C)) < 0) > 0.98;

valid = valid & (Vfront(32,32) > -2);
