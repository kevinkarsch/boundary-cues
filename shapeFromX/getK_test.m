
addpath(genpath('./matlabPyrTools/'));

clear all;

Z = 20*pyrBlur(rand(300,400),2, 'reflect1');

tic;
for i = 1:200
  [K, dK_Z] = getK(Z);
  d_loss = getK_backprop(sign(K), dK_Z);
end
toc


Z = 20*pyrBlur(rand(40,32),2, 'reflect1');
[K, dK_Z] = getK(Z);

d_loss = getK_backprop(sign(K), dK_Z);

step = 10^-5;

% loss = 0.5 * sum(K(:).^2);
loss = sum(abs(K(:)));
n_loss = zeros(size(Z));
for ii = 1:numel(Z)
  Z2 = Z;
  Z2(ii) = Z2(ii) + step;
  
  Kb = getK(Z2);
  lossb = sum(abs(Kb(:)));
%   lossb = 0.5 * sum(Kb(:).^2);
  n_loss(ii) = (lossb - loss)/step;
  
end

imagesc([d_loss, n_loss; d_loss - n_loss, 0*d_loss]); imtight;
err = max(max(abs([d_loss - n_loss])))

assert(err < (10*step))
