function d_loss_Z = getK_backprop(d_loss_K, dKZ, tap)

getK_filters;

% d_loss_K = d_loss_K(2:end-1, 2:end-1);

d_loss_K = d_loss_K ./ dKZ.denom;

% d_loss_Z = ...
%   conv3d(d_loss_K .* dKZ.f1,  f1m) + ...
%   conv3d(d_loss_K .* dKZ.f2,  f2m) + ...
%   conv3d(d_loss_K .* dKZ.f11, f11m) + ...
%   conv3d(d_loss_K .* dKZ.f22, f22m) + ...
%   conv3d(d_loss_K .* dKZ.f12, f12m);


d_loss_Z = unpad1( ...
  conv2(d_loss_K .* dKZ.f1,  f1m,  'full') + ...
  conv2(d_loss_K .* dKZ.f2,  f2m,  'full') + ...
  conv2(d_loss_K .* dKZ.f11, f11m, 'full') + ...
  conv2(d_loss_K .* dKZ.f22, f22m, 'full') + ...
  conv2(d_loss_K .* dKZ.f12, f12m, 'full') );


% 
% M = conv2(d_loss_K .* dKZ.f1,  f1m,  'full');
% M = M + conv2(d_loss_K .* dKZ.f2,  f2m,  'full');
% M = M + conv2(d_loss_K .* dKZ.f11, f11m, 'full');
% M = M + conv2(d_loss_K .* dKZ.f22, f22m, 'full');
% M = M + conv2(d_loss_K .* dKZ.f12, f12m, 'full');
% d_loss_Z = unpad1(M);

if nargin < 3
  tap = GETK_TAP;
end

assert( (tap == 3) || (tap == 4) || (tap == 5) )
if tap == 4
  assert(1==0);
  d_loss_Z = conv2(d_loss_Z, ones(2,2)/4, 'full');
  d_loss_Z = d_loss_Z(1:end-1,1:end-1);
elseif tap == 5
%   d_loss_Z = unpad1(conv2(d_loss_Z, ([1;2;1] * [1,2,1])/12, 'full'));
%   d_loss_Z = conv3d(d_loss_Z, ([1;2;1] * [1,2,1])/16);
  d_loss_Z = conv3d(d_loss_Z, ([1;2;1] * [1,2,1])/12);
end


% function d_loss_Z = getK_backprop(d_loss_K, dKZ, tap)
% 
% getK_filters;
% 
% % d_loss_K = d_loss_K(2:end-1, 2:end-1);
% 
% d_loss_K = d_loss_K ./ dKZ.denom;
% 
% d_loss_Z = ...
%   conv2(d_loss_K .* dKZ.f1,  f1m, 'same') + ...
%   conv2(d_loss_K .* dKZ.f2,  f2m, 'same') + ...
%   conv2(d_loss_K .* dKZ.f11, f11m, 'same') + ...
%   conv2(d_loss_K .* dKZ.f22, f22m, 'same') + ...
%   conv2(d_loss_K .* dKZ.f12, f12m, 'same');
% 
% if nargin < 3
%   tap = GETK_TAP;
% end
% 
% assert( (tap == 3) || (tap == 4) || (tap == 5) )
% if tap == 4
%   d_loss_Z = conv2(d_loss_Z, ones(2,2)/4, 'full');
%   d_loss_Z = d_loss_Z(1:end-1,1:end-1);
% elseif tap == 5
%   d_loss_Z = conv2([1,2,1]/4, [1,2,1]'/4, d_loss_Z, 'same');
% end




% function d_loss_Z = getK_backprop(d_loss_K, dKZ, tap)
% 
% getK_filters;
% 
% d_loss_K = d_loss_K(2:end-1, 2:end-1);
% 
% d_loss_K = d_loss_K ./ dKZ.denom;
% 
% d_loss_Z = ...
%   conv2(d_loss_K .* dKZ.f1,  f1m, 'full') + ...
%   conv2(d_loss_K .* dKZ.f2,  f2m, 'full') + ...
%   conv2(d_loss_K .* dKZ.f11, f11m, 'full') + ...
%   conv2(d_loss_K .* dKZ.f22, f22m, 'full') + ...
%   conv2(d_loss_K .* dKZ.f12, f12m, 'full');
% 
% if nargin < 3
%   tap = GETK_TAP;
% end
% 
% assert( (tap == 3) || (tap == 4) || (tap == 5) )
% if tap == 4
%   d_loss_Z = conv2(d_loss_Z, ones(2,2)/4, 'full');
%   d_loss_Z = d_loss_Z(1:end-1,1:end-1);
% elseif tap == 5
%   d_loss_Z = conv2([1,2,1]/4, [1,2,1]'/4, d_loss_Z, 'same');
% end


% d_loss_Z = d_loss_numer;
% d_loss_Z = d_loss_denom;


% d_loss_K_sub = 0.5 * d_loss_K;
% d_loss_N1 = conv2_rep( d_loss_K_sub, dKZ.f1);
% d_loss_N2 = conv2_rep( d_loss_K_sub, dKZ.f2);
% 
% d_loss_Z = conv2( dKZ.F1_1 .* d_loss_N1 - dKZ.F1_2 .* d_loss_N2, dKZ.f1, 'valid') ...
%          + conv2( dKZ.F2_2 .* d_loss_N2 - dKZ.F1_2 .* d_loss_N1, dKZ.f2, 'valid');
