function [err, us] = getError(us, gt)

V = ~isnan(us.height) & ~isnan(gt.height) & ~any(isnan(us.shading),3) & ~any(isnan(gt.shading),3) & ~any(isnan(us.albedo),3) & ~any(isnan(gt.albedo),3) & ~any(isnan(us.normal),3) & ~any(isnan(gt.normal),3);
V(:,[1,end]) = false;
V([1,end],:) = false;

V3 = repmat(V,[1,1,3]);  

errs_grosse = nan(1, size(us.shading,3));
assert(size(us.shading,3) == size(us.albedo,3))
for c = 1:size(us.shading,3)
  errs_grosse(c) = 0.5 * LMSE(us.shading(:,:,c), gt.shading(:,:,c), V) + 0.5 * LMSE(us.albedo(:,:,c), gt.albedo(:,:,c), V);
end
err.grosse = mean(errs_grosse);

alpha_shading = sum(gt.shading(V3) .* us.shading(V3), 1) ./ max(eps, sum(us.shading(V3) .* us.shading(V3)));
S = us.shading * alpha_shading;

% if isfield(us, 'light')
%   us.light(1,:) = us.light(1,:) + log(alpha_shading)/0.886227;
% end

alpha_albedo = sum(gt.albedo(V3) .* us.albedo(V3), 1) ./ max(eps, sum(us.albedo(V3) .* us.albedo(V3)));
A = us.albedo * alpha_albedo;

% log_shading = log(us.shading);
% w_shading = [log_shading(V3), ones(nnz(V3), 1)] \ log(gt.shading(V3));
% S = exp(log_shading * w_shading(1) + w_shading(2));
% 
% log_albedo = log(us.albedo);
% w_albedo = [log_albedo(V3), ones(nnz(V3), 1)] \ log(gt.albedo(V3));
% A = exp(log_albedo * w_albedo(1) + w_albedo(2));

err.shading = mean((S(V3) - gt.shading(V3)).^2);
err.albedo = mean((A(V3) - gt.albedo(V3)).^2);

us.shading = S;
us.albedo = A;

% n1 = us.normal(:,:,1) ./ us.normal(:,:,3);
% n2 = us.normal(:,:,2) ./ us.normal(:,:,3);
% n1 = n1(V);
% n2 = n2(V);
% n = [n1, n2];
% 
% Nt = reshape(gt.normal(V3), [], 3);
% 
% shift = [0,0];
% for ii = 1:3
%   [shift, f_trace] = minimize_conjugate(shift, 'lossfun_normal_error', 100, n, Nt, struct('USE_NUMERICAL', true));
% end
% 
% [loss, d_loss, nb] = lossfun_normal_error(shift, n, Nt, []);
% n1 = nan(size(us.height));
% n2 = nan(size(us.height));
% n1(V) = n(:,1) + shift(1);
% n2(V) = n(:,2) + shift(2);
% 
% N3 = 1./sqrt(n1.^2 + n2.^2 + 1);
% N1 = n1 .* N3;
% N2 = n2 .* N3;
% 
% us.normal = cat(3, N1, N2, N3);
% err.normal = loss;

% nb = [n(:,1) + shift(1), n(:,2) + shift(2)];
% N3 = 1./sqrt(nb(:,1).^2 + nb(:,2).^2 + 1);
% N1 = nb(:,1) .* N3;
% N2 = nb(:,2) .* N3;
% loss = mean(acos(sum([N1, N2, N3] .* Nt,2)).^2);


% shifts = [-2:.1:2];
% err = nan(length(shifts), length(shifts));
% for shift_i = 1:length(shifts)
%   for shift_j = 1:length(shifts)
%     shift1 = shifts(shift_i);
%     shift2 = shifts(shift_j);
%     n1b = n1 + shift1;
%     n2b = n2 + shift2;
%     N3 = 1./sqrt(n1b.^2 + n2b.^2 + 1);
%     N1 = n1b .* N3;
%     N2 = n2b .* N3;
%     err(shift_i, shift_j) = mean(acos(sum([N1, N2, N3] .* Nt,2)).^2);
%   end
% end

% err.normal = min(err(:));

% [err.height] = getZMSE2(us.height, gt.height);
% [err.height, us.height] = getZMSE1(us.height, gt.height);


% K = getK(us.height);
% 
% if isfield(gt, 'curvature')
%   K_true = gt.curvature;
% else
%   K_true = getK(gt.height);
% end
% v = ~isnan(K) & ~isnan(K_true);
% K_delta = K(v) - K_true(v);
% % err.curvature = mean(abs(K_delta));
% % keyboard
% err.curvature = mean(K_delta.^2);


us.normal = getNormals_conv(us.height);

d = us.normal .* gt.normal;
d = reshape(d, [], 3);
d = d(all(~isnan(d),2),:);
% err.normal = mean((360/pi*acos(sum(d,2))).^2);
err.normal = mean(acos(sum(d,2)).^2);


Lv_gt = visSH_color(gt.light, [256,256]);
Lv_us = visSH_color(us.light, [256,256]);

VL = ~isnan(Lv_gt);
alpha_light = sum(Lv_gt(VL) .* Lv_us(VL), 1) ./ max(eps, sum(Lv_us(VL) .* Lv_us(VL),1));
Lv_us = Lv_us * alpha_light;

err.light = mean((Lv_us(VL) - Lv_gt(VL)).^2);

errs = [];
criteria = {'light', 'normal', 'albedo', 'shading', 'grosse'};
for s = criteria
  if err.(s{1}) > 10^-10 % Assume perfect, and therefore given
    errs(end+1) = err.(s{1});
  end
end
err.avg = exp(mean(log(errs)));








% function [err, us] = getError(us, gt)
% 
% % BE_COLOR_SENSITIVE = true;
% % % BE_COLOR_SENSITIVE = false;
% 
% V = ~isnan(us.height) & ~isnan(gt.height) & ~any(isnan(us.shading),3) & ~any(isnan(gt.shading),3) & ~any(isnan(us.albedo),3) & ~any(isnan(gt.albedo),3) & ~any(isnan(us.normal),3) & ~any(isnan(gt.normal),3);
% V(:,[1,end]) = false;
% V([1,end],:) = false;
% 
% V3 = repmat(V,[1,1,3]);
% 
% % Color sensitive error metrics
% alpha_shading = sum(gt.shading(V3) .* us.shading(V3), 1) ./ max(eps, sum(us.shading(V3) .* us.shading(V3)));
% S = us.shading * alpha_shading;
% 
% alpha_albedo = sum(gt.albedo(V3) .* us.albedo(V3), 1) ./ max(eps, sum(us.albedo(V3) .* us.albedo(V3)));
% A = us.albedo * alpha_albedo;
% 
% % S = us.shading;
% % A = us.albedo;
% 
% err.shading = mean((S(V3) - gt.shading(V3)).^2);
% err.albedo = mean((A(V3) - gt.albedo(V3)).^2);
%   
% 
% % Color insensitive error metrics
% errs_grosse = nan(1, size(us.shading,3));
% % errs_shading = nan(1, size(us.shading,3));
% % errs_albedo = nan(1, size(us.shading,3));
% assert(size(us.shading,3) == size(us.albedo,3))
% for c = 1:size(us.shading,3)
% %   S = us.shading(:,:,c);
% %   St = gt.shading(:,:,c);
% %   a = sum(St(V) .* S(V), 1) ./ max(eps, sum(S(V) .* S(V)));
% %   S = S * a;
% %   errs_shading(c) = mean((S(V) - St(V)).^2);
% % 
% % %     us.shading(:,:,c) = S;
% % 
% %   A = us.albedo(:,:,c);
% %   At = gt.albedo(:,:,c);
% %   a = sum(At(V) .* A(V), 1) ./ max(eps, sum(A(V) .* A(V)));
% %   A = A * a;
% % 
% % %     us.albedo(:,:,c) = A;
% % 
% %   errs_albedo(c) = mean((A(V) - At(V)).^2);
% 
%   errs_grosse(c) = 0.5 * LMSE(us.shading(:,:,c), gt.shading(:,:,c), V) + 0.5 * LMSE(us.albedo(:,:,c), gt.albedo(:,:,c), V);
% end
% 
% 
% err.grosse = mean(errs_grosse);
% 
% % err.inv_shading = mean(errs_shading);
% % err.inv_albedo = mean(errs_albedo);
% 
% % err.survey = getIMSE(us, gt);
% 
% 
% % d = us.normal - gt.normal;
% % d = reshape(d, [], 3);
% % d = d(all(~isnan(d),2),:);
% % err.normal = mean(sum(d.^2,2),1);
% 
% d = us.normal .* gt.normal;
% d = reshape(d, [], 3);
% d = d(all(~isnan(d),2),:);
% err.normal = mean(acos(sum(d,2)).^2);
% % err.normal = mean(acos(sum(d,2)));
% 
% 
% [err.height, us.height] = getZMSE1(us.height, gt.height);
% 
% 
% errs = [];
% % criteria = {'height', 'survey', 'col_albedo', 'col_shading', 'col_grosse'};
% criteria = {'height', 'normal', 'albedo', 'shading', 'grosse'};
% for s = criteria
%   errs(end+1) = err.(s{1});
% end
% err.avg = exp(mean(log(errs)));
% 
% 
% 
% 
% 
% 
% 
% % function [err, us] = getError(us, gt)
% % 
% % % BE_COLOR_SENSITIVE = true;
% % % % BE_COLOR_SENSITIVE = false;
% % 
% % V = ~isnan(us.height) & ~isnan(gt.height) & ~any(isnan(us.shading),3) & ~any(isnan(gt.shading),3) & ~any(isnan(us.albedo),3) & ~any(isnan(gt.albedo),3) & ~any(isnan(us.normal),3) & ~any(isnan(gt.normal),3);
% % V(:,[1,end]) = false;
% % V([1,end],:) = false;
% % 
% % V3 = repmat(V,[1,1,3]);
% % 
% % if BE_COLOR_SENSITIVE
% %   
% %   % Color sensitive error metrics
% % 
% %   alpha_shading = sum(gt.shading(V3) .* us.shading(V3), 1) ./ max(eps, sum(us.shading(V3) .* us.shading(V3)));
% %   S = us.shading * alpha_shading;
% %   
% %   alpha_albedo = sum(gt.albedo(V3) .* us.albedo(V3), 1) ./ max(eps, sum(us.albedo(V3) .* us.albedo(V3)));
% %   A = us.albedo * alpha_albedo;
% %   
% %   err.shading = mean((us.shading(V3) - gt.shading(V3)).^2);
% %   err.albedo = mean((us.albedo(V3) - gt.albedo(V3)).^2);
% %   err.grosse = CLMSE(us.shading, gt.shading, V3);
% %   
% % %   us.shading = S;
% % %   us.albedo = A;
% %   
% % else
% %   
% %   % Color insensitive error metrics
% %   errs_grosse = nan(1, size(us.shading,3));
% %   errs_shading = nan(1, size(us.shading,3));
% %   errs_albedo = nan(1, size(us.shading,3));
% %   assert(size(us.shading,3) == size(us.albedo,3))
% %   for c = 1:size(us.shading,3)
% %     S = us.shading(:,:,c);
% %     St = gt.shading(:,:,c);
% %     a = sum(St(V) .* S(V), 1) ./ max(eps, sum(S(V) .* S(V)));
% %     S = S * a;
% %     errs_shading(c) = mean((S(V) - St(V)).^2);
% %     
% % %     us.shading(:,:,c) = S;
% %     
% %     A = us.albedo(:,:,c);
% %     At = gt.albedo(:,:,c);
% %     a = sum(At(V) .* A(V), 1) ./ max(eps, sum(A(V) .* A(V)));
% %     A = A * a;
% %     
% % %     us.albedo(:,:,c) = A;
% %     
% %     errs_albedo(c) = mean((A(V) - At(V)).^2);
% %     
% %     errs_grosse(c) = 0.5 * LMSE(us.shading(:,:,c), gt.shading(:,:,c), V) + 0.5 * LMSE(us.albedo(:,:,c), gt.albedo(:,:,c), V);
% %   end
% %   err.grosse = mean(errs_grosse);
% %   err.shading = mean(errs_shading);
% %   err.albedo = mean(errs_albedo);
% %   
% % end
% % 
% % err.survey = getIMSE(us, gt);
% % 
% % [err.height, us.height] = getZMSE1(us.height, gt.height);
% % 
% % 
% % % d = reshape(us.normal - gt.normal, [], 3);
% % % d = d(all(~isnan(d),2),:);
% % % err.normal = mean(sum(d.^2,2),1);
% % % 
% % % St = gt.shading;
% % % At = gt.albedo;
% % % St(repmat(isnan(us.height), [1,1,3])) = nan;
% % % At(repmat(isnan(us.height), [1,1,3])) = nan;
% % % gS = {};
% % % gSt = {};
% % % gA = {};
% % % gAt = {};
% % % for c = 1:3
% % %   gS{c} = getGradNorm(us.shading(:,:,c));
% % %   gA{c} = getGradNorm(us.albedo(:,:,c));
% % %   gSt{c} = getGradNorm(St(:,:,c));
% % %   gAt{c} = getGradNorm(At(:,:,c));
% % % end
% % % 
% % % gSt = cat(3,gSt{:});
% % % gAt = cat(3,gAt{:});
% % % gS = cat(3,gS{:});
% % % gA = cat(3,gA{:});
% % % 
% % % vS = ~isnan(gSt) | ~isnan(gS);
% % % err.g_shading = mean((gSt(vS) - gS(vS)).^2);
% % % 
% % % vA = ~isnan(gAt) | ~isnan(gA);
% % % err.g_albedo = mean((gAt(vA) - gA(vA)).^2);
% % 
% % 
% % errs = [];
% % for s = {'albedo', 'shading', 'survey', 'grosse', 'height'}
% %   errs(end+1) = err.(s{1});
% % end
% % err.avg = exp(mean(log(errs)));
% % 
% % % [err.affine_height, us.affine_height] = getZMSE2(us.height, gt.height);
