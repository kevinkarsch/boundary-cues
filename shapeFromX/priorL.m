function [loss_L, d_loss_L, losses_L] = priorL(L, data, params)

losses_L.natural_color.gaussian = 0;
losses_L.lab_color.gaussian = 0;
% losses_L.both_color.gaussian = 0;

d_loss_L = zeros(size(L));

mult_gaussian = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).gaussian{1};
mult_gaussian_mod = 1/numel(L);
if mult_gaussian ~= 0
  [l, dl] = lmvnpdf2(L(:)', data.L_gaussian.mu, data.L_gaussian.Sigma, 1);
  losses_L.(params.L_WHITEN_PARAMS).gaussian = mult_gaussian_mod*l;
  d_loss_L = d_loss_L + (mult_gaussian * mult_gaussian_mod) * reshape(dl, size(L));
end

loss_L = mult_gaussian * (losses_L.natural_color.gaussian + losses_L.lab_color.gaussian);

losses_L.natural_color.GSM= 0;
losses_L.lab_color.GSM = 0;
% losses_L.both_color.GSM = 0;

mult_GSM = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).GSM{1};
mult_GSM_mod = 1/numel(L);
if mult_GSM ~= 0
  [l, dl] = GSM_mvn_pdf(data.L_GSM, L(:)', 2);
  losses_L.(params.L_WHITEN_PARAMS).GSM = mult_GSM_mod*l;
  d_loss_L = d_loss_L + (mult_GSM * mult_GSM_mod) * reshape(dl, size(L));
end

loss_L = loss_L + mult_GSM * (losses_L.natural_color.GSM + losses_L.lab_color.GSM);% + losses_L.both_color.GSM);

losses_L.natural_color.MOG= 0;
losses_L.lab_color.MOG = 0;
% losses_L.both_color.GSM = 0;

mult_MOG = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).MOG{1};
mult_MOG_mod = 1/numel(L);
if mult_MOG ~= 0
  [l, dl] = MOG_mvn_pdf(L(:)', data.L_MOG, 1);
%   checkgrad(L(:)', 'MOG_mvn_pdf', 10^-5, data.L_MOG, 1)
  losses_L.(params.L_WHITEN_PARAMS).MOG = mult_MOG_mod*l;
  d_loss_L = d_loss_L + (mult_MOG * mult_MOG_mod) * reshape(dl, size(L));
end

loss_L = loss_L + mult_MOG * (losses_L.natural_color.MOG + losses_L.lab_color.MOG);% + losses_L.both_color.GSM);


% losses_L.natural_color.smooth = 0;
% losses_L.lab_color.smooth = 0;
% 
% mult_smooth = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).smooth{1};
% 
% mult_mod_smooth = 1/(numel(L) * size(data.light_M,1));
% 
% if mult_smooth ~= 0
%   model = data.prior.ML.(params.L_WHITEN_PARAMS).GSM;
%   sz = [128, 128];
%   
%   min_sz = min(sz);
%   j = ([1:sz(1)]-min_sz/2 - (sz(1)-min_sz)/2)*(2/min_sz);
%   i = ([1:sz(2)]-min_sz/2 - (sz(2)-min_sz)/2)*(2/min_sz);
%   [Y,X] = meshgrid(i,j);
%   Z = sqrt(max(0, 1-(X.^2 + Y.^2)));
%   valid = Z ~= 0;
%   
%   Es = {};
%   dE_Ls = {};
%   for c = 1:3
%     [E, junk, dE_L] = renderSH_helper([X(valid), Y(valid), Z(valid)], L(:,c));
%     Es{c} = E;
%     dE_Ls{c} = dE_L;
%   end
%   Ev = cat(2,Es{:})';
%   
% %   M = data.light_M;
%   M_T = data.light_M_T;
%   M_T = M_T(valid,:);
%   
%   [U, S, V] = svd(model.Sigma_inv);
%   W = U * sqrt(S) * V';
%   
%   WEv = W*Ev;
%   MAW = WEv * M_T;
%   mahal_dist = 0.5 * sum(MAW.^2,1)';
%   
%   F = model.LL_zero - model.lut.F;
%   [LL, dLL_dist] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, F);
%   
%   dLL_E = (M_T * bsxfun(@times, dLL_dist, MAW')) * W;
%   
%   losses_L.(params.L_WHITEN_PARAMS).smooth = mult_mod_smooth * sum(LL(:));
%   
%   dLL_L = zeros(size(L));
%   for c = 1:3
%     dLL_L(:,c) = (dLL_E(:,c)' * dE_Ls{c})';
%   end
%   
%   d_loss_L = d_loss_L + (mult_mod_smooth * mult_smooth) * dLL_L;
% 
% end
  

% loss_L = mult_gaussian * (losses_L.natural_color.gaussian + losses_L.lab_color.gaussian) + mult_smooth * (losses_L.natural_color.smooth + losses_L.lab_color.smooth);


% if mult_smooth ~= 0
%   model = data.prior.ML.(params.L_SMOOTH_MODEL).GSM;
%   sz = [128, 128];
%   
%   min_sz = min(sz);
%   j = ([1:sz(1)]-min_sz/2 - (sz(1)-min_sz)/2)*(2/min_sz);
%   i = ([1:sz(2)]-min_sz/2 - (sz(2)-min_sz)/2)*(2/min_sz);
%   [Y,X] = meshgrid(i,j);
%   Z = sqrt(max(0, 1-(X.^2 + Y.^2)));
%   valid = Z ~= 0;
%   
%   Es = {};
%   dE_Ls = {};
%   for c = 1:3
%     [E, junk, dE_L] = renderSH_helper([X(:), Y(:), Z(:)], L(:,c));
%     Es{c} = E;
%     dE_Ls{c} = dE_L;
%   end
%   Ev = cat(2,Es{:})';
%   
%   M = data.light_M;
%   M_T = data.light_M_T;
%   
%   [U, S, V] = svd(model.Sigma_inv);
%   W = U * sqrt(S) * V';
%   
%   WEv = W*Ev;
%   MAW = WEv * M_T;
%   mahal_dist = 0.5 * sum(MAW.^2,1)';
%   
%   F = model.LL_zero - model.lut.F;
%   [LL, dLL_dist] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, F);
%   
%   dLL_E = (M_T * bsxfun(@times, dLL_dist, MAW')) * W;
%   
%   losses_L.smooth = mult_mod_smooth * sum(LL(:));
%   
%   dLL_L = zeros(size(L));
%   for c = 1:3
%     dLL_L(:,c) = (dLL_E(:,c)' * dE_Ls{c})';
%   end
%   
%   d_loss_L = d_loss_L + (mult_mod_smooth * mult_smooth) * dLL_L;
% 
% end
%   
% %   L = log(visSH_color(L, SZ));
% %   ML = 
%   
% 
% loss_L = mult_gaussian * (losses_L.natural_color + losses_L.lab_color) + mult_smooth * (losses_L.smooth);
