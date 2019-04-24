function err = LMSE(X_us, X_true, X_mask, k)

if nargin < 3
  mask = true(size(X_us));
end

if (nargin < 4) || isempty(k)
 k = 20;
end


X_us(~X_mask) = 0;
X_true(~X_mask) = 0;

C_true = [im2col(X_true, [k,k], 'distinct'), im2col(X_true(:,(k/2):end), [k,k], 'distinct'), im2col(X_true((k/2):end,:), [k,k], 'distinct'), im2col(X_true((k/2):end,(k/2):end), [k,k], 'distinct')];
C_us = [im2col(X_us, [k,k], 'distinct'), im2col(X_us(:,(k/2):end), [k,k], 'distinct'), im2col(X_us((k/2):end,:), [k,k], 'distinct'), im2col(X_us((k/2):end,(k/2):end), [k,k], 'distinct')];
C_mask = [im2col(X_mask, [k,k], 'distinct'), im2col(X_mask(:,(k/2):end), [k,k], 'distinct'), im2col(X_mask((k/2):end,:), [k,k], 'distinct'), im2col(X_mask((k/2):end,(k/2):end), [k,k], 'distinct')];

% C_true = im2col(X_true, [k,k], 'distinct');
% C_us = im2col(X_us, [k,k], 'distinct');
% C_mask = im2col(X_mask, [k,k], 'distinct');

C_alpha = sum(C_mask .* C_true .* C_us, 1) ./ max(eps, sum(C_mask .* C_us .* C_us));
numer = sum(C_mask .* (C_true - repmat(C_alpha, size(C_true,1), 1) .* C_us).^2,1);
denom = sum(C_mask .* (C_true.^2),1);
err = sum(numer) ./ sum(denom);



% X_true_bak = X_true;
% X_us_bak = X_us;
% 
% err = 0;
% for offset_i = [0, 1];
%   for offset_j = [0, 1];
%     
%     X_true = X_true_bak;
%     X_us = X_us_bak;
% 
%     X_true = X_true((1+offset_i/2 * k):end, (1+offset_j/2 * k):end);
%     X_us = X_us((1+offset_i/2 * k):end, (1+offset_j/2 * k):end);
%     
%     new_size = floor(size(X_true)/k) * k;
% 
%     X_true = X_true(1:new_size(1), 1:new_size(2));
%     X_us = X_us(1:new_size(1), 1:new_size(2));
% 
%     C_true = im2col(X_true, [k,k], 'distinct');
%     C_us = im2col(X_us, [k,k], 'distinct');
%     C_us = im2col(X_us, [k,k], 'distinct');
% 
%     C_alpha = sum(C_true .* C_us, 1) ./ max(eps, sum(C_us.^2));
%     C_alpha = repmat(C_alpha, size(C_true,1), 1);
%     e = sum(mean((C_true - C_alpha .* C_us).^2,1),2);
%     
%     err = err + e;
% 
%   end
% end
% err = err / 4;
% 
% 
