function d_loss_Z = getK_backprop_fast_mat(d_loss_K, dKZ)

getK_filters;

d_loss_K = d_loss_K ./ dKZ(:,:,1);

fs = {f1m, f2m, f11m, f22m, f12m};

d_loss_Z = 0;
for fi = 1:length(fs)
  d_loss_Z = d_loss_Z + conv2(d_loss_K .* dKZ(:,:,fi+1),  fs{fi},  'full');
end
d_loss_Z = unpad1(d_loss_Z);



% function d_loss_Z = getK_backprop(d_loss_K, dKZ)
% 
% getK_filters;
% 
% % d_loss_Kv = d_loss_K(:)' ./ dKZ(1,:);
% % d_loss_Kv = reshape(permute(bsxfun(@times, d_loss_Kv, dKZ(2:end,:)), [2,1]), [size(d_loss_K), 5]);
% % 
% % d_loss_Z = unpad1( ...
% %   conv2(d_loss_Kv(:,:,1),  f1m,  'full') + ...
% %   conv2(d_loss_Kv(:,:,2),  f2m,  'full') + ...
% %   conv2(d_loss_Kv(:,:,3), f11m, 'full') + ...
% %   conv2(d_loss_Kv(:,:,4), f22m, 'full') + ...
% %   conv2(d_loss_Kv(:,:,5), f12m, 'full') );
% 
% 
% 
% dKZ = reshape(permute(dKZ, [2,1]), [size(d_loss_K), 6]);
% 
% d_loss_K = d_loss_K ./ dKZ(:,:,1);
% 
% d_loss_Z = unpad1( ...
%   conv2(d_loss_K .* dKZ(:,:,2),  f1m,  'full') + ...
%   conv2(d_loss_K .* dKZ(:,:,3),  f2m,  'full') + ...
%   conv2(d_loss_K .* dKZ(:,:,4), f11m, 'full') + ...
%   conv2(d_loss_K .* dKZ(:,:,5), f22m, 'full') + ...
%   conv2(d_loss_K .* dKZ(:,:,6), f12m, 'full') );
% 
