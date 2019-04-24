function [loss_A, d_loss_A, losses_A] = priorA(A, data, params)

d_loss_A = zeros(size(A));

mult_smooth = params.MULT_OPTS.saifs.albedo.smooth{1};
mult_mod_smooth = (3/((params.A_MEDIAN_HALFWIDTH*2+1)^2))/numel(A);


if mult_smooth ~= 0

  if params.ALBEDO_SMOOTH_MODEL > 0
    
    M = data.A_median_filter_mat;
    M_T = data.A_median_filter_mat_T;

    Av = reshape(A, [], 3);
    MA = M * Av;

    low = data.prior.albedo.MA_hist_low;
    high = data.prior.albedo.MA_hist_high;
    MA_splat = splat3_fast_wrapper(MA, low, high);
    
    F = data.prior.albedo.MA_hist;
    v = sum(sum(sum(F .* MA_splat.N)));
    dv = F;
    
    d_loss_MA = splat3_backprop_fast_wrapper(dv, MA_splat);
    
    losses_A.smooth = mult_mod_smooth * v;
    d_loss_A = d_loss_A + reshape(M_T * ((mult_mod_smooth * mult_smooth) * d_loss_MA), size(A));
    
    
    %   F = model.LL_zero - model.lut.F;
    %   [LL, dLL_dist] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, F);
    %
    %   dLL = M_T * bsxfun(@times, dLL_dist, MA_icov);
    %
    %   losses_A.MA3 = mult*sum(LL(:));
    %   loss_A = loss_A + losses_A.MA3;
    %
    %   d_loss_A = d_loss_A + mult * reshape(dLL, size(A));
    
  else
    
    %   M = data.A_median_filter_mat;
    M_T = data.A_median_filter_mat_T;
    
    model = data.prior.albedo.MA.GSM_mvn;
    %   model = data.prior.albedo.MA.GSM_mvn{params.A_MEDIAN_HALFWIDTH};
    
    %   Av = reshape(A, [], 3);
    %   MA = M * Av;
    %   [LL, dLL] = GSM_mvn_pdf(model, MA, 2);
    %
    %   losses_A.MA3 = mult*sum(LL(:));
    %   loss_A = loss_A + losses_A.MA3;
    %
    %   d_loss_A = d_loss_A + mult * reshape(M_T * dLL, size(A));
    
    %   Av = reshape(A, [], 3);
    %   MA = M * Av;
    %
    %   MA_icov = MA * model.Sigma_inv;
    %   mahal_dist = 0.5 * sum(MA_icov .* MA,2);
    %
    %   F = model.LL_zero - model.lut.F;
    %   [LL, dLL_dist] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, F);
    %
    %   dLL = M_T * bsxfun(@times, dLL_dist, MA_icov);
    %
    %   losses_A.MA3 = mult*sum(LL(:));
    %   loss_A = loss_A + losses_A.MA3;
    %
    %   d_loss_A = d_loss_A + mult * reshape(dLL, size(A));
    
    
    %   Av = reshape(A, [], 3);
    %   [U, S, V] = svd(model.Sigma_inv);
    %   W = U * sqrt(S) * V';
    %
    %   MAW = M * (Av * W);
    %   mahal_dist = 0.5 * sum(MAW.^2,2);
    %
    %   F = model.LL_zero - model.lut.F;
    %   [LL, dLL_dist] = interp1_fixed_sum_fast(mahal_dist, model.lut.bin_range(1), model.lut.bin_width, F);
    %
    %   dLL = (M_T * bsxfun(@times, dLL_dist, MAW)) * W;
    %
    %   losses_A.MA3 = mult*sum(LL(:));
    %   loss_A = loss_A + losses_A.MA3;
    %
    %   d_loss_A = d_loss_A + mult * reshape(dLL, size(A));
    
    
    
    %   [U, S, V] = svd(model.Sigma_inv);
    %   W = U * sqrt(S) * V';
    W = model.Sigma_whiten;
    
    Av = reshape(A, [], 3)';
    WAv = W*Av;
    MAW = WAv * M_T;
    mahal_dist = 0.5 * sum(MAW.^2,1)';
    
    epsilon = params.A_SMOOTH_EPSILON{1};
    if epsilon > 0
      mahal_dist_soft = sqrt(mahal_dist.^2 + epsilon.^2);
    else
      mahal_dist_soft = mahal_dist;
    end
    
    F = model.LL_zero - model.lut.F;
%mahal_dist_soft(isnan(mahal_dist_soft)|isinf(mahal_dist_soft)) = 0;
    [LL, dLL_dist] = interp1_fixed_sum_fast(mahal_dist_soft, model.lut.bin_range(1), model.lut.bin_width, F);
    
    if epsilon > 0
      dLL_dist = dLL_dist .* (mahal_dist ./ max(eps,mahal_dist_soft));
    end
    
    dLL = (M_T * bsxfun(@times, dLL_dist, MAW')) * W;
    
    losses_A.smooth = mult_mod_smooth * sum(LL(:));
    
    d_loss_A = d_loss_A + (mult_mod_smooth * mult_smooth) * reshape(dLL, size(A));
    
  end
  
else
  losses_A.smooth = 0;
end


  
mask = data.valid;
mask_rep = repmat(mask, [1,1,3]);
  
Av = reshape(A(mask_rep), [], 3);
d_loss_Av = zeros(size(Av));

mult_entropy = params.MULT_OPTS.saifs.albedo.entropy{1};
mult_mod_entropy = 1;

mult_hist = params.MULT_OPTS.saifs.albedo.hist{1};
mult_mod_hist = 1 / (8*size(Av,1));

if (mult_entropy ~= 0) || (mult_hist ~= 0)
  
  sigma = params.MULT_OPTS.saifs.albedo.entropy_sigma{1};
  C = data.prior.albedo.A_whiten;
  
  Av_white = Av * C;
  
  low =  params.ALBEDO_BIN_LOW;
  high = params.ALBEDO_BIN_HIGH;

%   [prctile(Av_white, 1); prctile(Av_white, 99)]
  
  dV = 0;
  A_splat = splat3_fast_wrapper(Av_white, low, high);
  
  if mult_entropy ~= 0
    [v, dv] = renyi3_fixed_sub(Av_white, A_splat, sigma, 2);
    
    losses_A.entropy = mult_mod_entropy * v;
    dV = dV + (mult_mod_entropy * mult_entropy) * dv;
    
  else
    losses_A.entropy = 0;
  end
  
  if mult_hist ~= 0
    
    F = data.prior.albedo.Aw_hist;
    v = sum(sum(sum(F .* A_splat.N)));
    dv = F;
    
    losses_A.hist = mult_mod_hist * v;
    dV = dV + (mult_mod_hist * mult_hist) * dv;
    
  else
    losses_A.hist = 0;
  end
    
  d_loss_Av = d_loss_Av + splat3_backprop_fast_wrapper(dV, A_splat) * C;
    
else
  losses_A.entropy = 0;
  losses_A.hist = 0;
end


loss_A = mult_smooth * losses_A.smooth + mult_entropy * losses_A.entropy + mult_hist * losses_A.hist;

d_loss_A(mask_rep) = d_loss_A(mask_rep) + d_loss_Av(:);


% mult = params.MULT_OPTS.srfs.albedo.hist{1} / (8*size(Av,1));
% if mult ~= 0
%   
%   X_splat = data.prior.albedo.A_hist{params.ALBEDO_HIST_MODEL};
%   sigma = 0.15;
%   
% %   A_splat = splat3(Av, sigma, X_splat.bin_range_low, X_splat.bin_range_high);
%   A_splat = splat3_fast_wrapper(Av, sigma, X_splat.bin_range_low, X_splat.bin_range_high);
%   
%   V = sum(sum(sum(X_splat.F .* A_splat.N)));
%     
% %   dV = splat3_backprop(X_splat.F, A_splat);
%   dV = splat3_backprop_fast_wrapper(X_splat.F, A_splat);
%   
%   loss_A = loss_A + mult * V;
%   d_loss_Av = d_loss_Av + mult * dV;
%   
% end



% mult = params.MULT_OPTS.srfs.albedo.gamut{1};
% 
% if mult ~= 0
%   
%   hull_A = data.prior.albedo.RGB_hull.A;
%   hull_b = data.prior.albedo.RGB_hull.b;
% 
%   [D_max, D_max_idx] = hullDist(Av', hull_A', hull_b);
%   
%   n = size(D_max,1);
%   
%   if params.GAMUT_POWER == 2
%     
%     loss_A = loss_A + mult/n * ((D_max.^2) * W);
%     d_loss_Av = d_loss_Av + sparse(1:size(D_max,2), D_max_idx, ( (2*mult/n) * D_max .* W'), size(D,2), size(D,1)) * hull_A;
%     
%   elseif params.GAMUT_POWER == 1
%     
%     losses_A.gamut = mult/n * sum(abs(D_max));
%     above = D_max > 0;
%     d_loss_Av(above,:) = d_loss_Av(above,:) + hull_A(D_max_idx(above),:) .* repmat(( (mult/n) * sign(D_max(above))), [1,3]);
%     
%   end
%   
%   loss_A = loss_A + losses_A.gamut;
% 
% end
% 
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.white1{1}/32;
% 
% if mult ~= 0
%  
%   n = size(Av,1);
%   
%   mag_sq = sum(Av.^2 + 10^-4,2);
%   mag = sqrt(mag_sq);
%     
%   losses_A.white1 = mult/n * sum(mag);
%   d_loss_Av = d_loss_Av + (mult/n) * Av .* repmat(1 ./ mag, [1,3]);
% 
%   loss_A = loss_A + losses_A.white1;
% 
% end
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gray{1}/32;
% 
% if mult ~= 0
%   u = ones(3,1)/sqrt(3);
%   R = (eye(3,3) - u*u');
%   Av_resid = Av * R;
%   
%   n = size(Av,1);
% 
%   mag_sq = sum(Av_resid.^2 + 10^-4,2);
%   mag = sqrt(mag_sq);
%     
%   losses_A.gray = mult/n * sum(mag);
%   d_loss_Av_resid = (mult/n) * Av .* repmat(1 ./ mag, [1,3]);
%   d_loss_Av = d_loss_Av + d_loss_Av_resid * R;
% 
%   loss_A = loss_A + losses_A.gray;
% 
% end
% 
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gray2{1}/32;
% 
% if mult ~= 0
%   
%   y = -ones(3,1)/sqrt(3);
%   uv = [[1; -2; 1]/sqrt(6), [-1; 0; 1]/sqrt(2)];
%   
%   Y = max(0, Av* y);
%   UV = Av * uv;
%   
%   Y_soft_sq = Y.^2 + 0.01;
%   Y_soft = sqrt(Y_soft_sq);
%   UV_norm = bsxfun(@rdivide, UV, Y_soft);
%   
%   GSM = data.prior.albedo.gray.UV_GSM;
%   
%   [L, dL] = GSM_mvn_pdf(GSM, UV_norm, 2);
%     
%   losses_A.gray2 = mult/n * sum(L);
%   
%   dL = (mult/n) * dL;
%   d_loss_Av = d_loss_Av + bsxfun(@rdivide, dL, Y_soft) * uv';
%   d_loss_Av = d_loss_Av + sum(dL .* bsxfun(@times, -Y./(Y_soft .* Y_soft_sq), UV),2) * y';
% 
%   loss_A = loss_A + losses_A.gray2;
% 
% end
% 
% d_loss_A(mask_rep) = d_loss_A(mask_rep) + d_loss_Av(:);





% function [loss_A, d_loss_A, losses_A] = priorA(A, data, params)
% 
% Apyr = buildGpyr_color(A, params.APYR_DEPTH{1});
% 
% d_loss_Apyr = {};
% for s = 1:params.APYR_DEPTH{1}
%   d_loss_Apyr{s} = zeros(size(Apyr{s}));
% end
% 
% loss_A = 0;
% 
% losses_A.MA3 = 0;
% 
% for s = 1:params.APYR_DEPTH_SMOOTH{1}
%   
%   Aband = Apyr{s};
% 
%   mult = (3/25)*params.MULT_OPTS.srfs.albedo.MA3{1}/numel(A) * (params.APYR_ROOT{1}.^(s-1)) / params.APYR_DEPTH_SMOOTH{1};
%   
%   if mult ~= 0
%     
%     M = data.median_filter_mats{s};
%     MA = M * reshape(Aband, [], 3);
% 
%     [LL, dLL] = GSM_mvn_pdf(data.prior.albedo.MA.GSM_mvn{s}, MA, 2);
% %     loss_A = loss_A + mult * sum(LL(:));
%     losses_A.MA3 = losses_A.MA3 + mult*sum(LL(:));
%     dL_A = mult * reshape(M' * dLL, size(Aband));
% 
%     d_loss_Apyr{s} = d_loss_Apyr{s} + dL_A;
%     
%   end
% 
% end
% 
% loss_A = loss_A + losses_A.MA3;
% 
% 
% 
% 
% 
% Av = {};
% W = {};
% for s = 1:params.APYR_DEPTH_ENTROPY{1}
%   
%   Aband = Apyr{s};
%   mask = data.valid_pyr{s};
%   mask_rep = repmat(mask, [1,1,3]);
%   
%   Av{s} = reshape(Aband(mask_rep), [], 3);
%   W{s} = (4.^(s-1)) * ones(size(Av{s},1), 1);
% 
% end
% 
% Av = cat(1, Av{:});
% W  = cat(1, W{:});
% 
% d_loss_Av = zeros(size(Av));
% 
% mult = params.MULT_OPTS.srfs.albedo.multi_entropy{1};
% 
% if mult ~= 0
%   
%   sigma = params.RENYI_SIGMA{1};
%   C = data.prior.albedo.A_whiten;
%   
%   if params.WHITEN_ALBEDO_ENTROPY
%     Av_white = Av * C;
%   else
%     Av_white = Av;
%   end
%     
%   [V, dV1, dV2, dV3] = renyi3W(Av_white(:,1), Av_white(:,2), Av_white(:,3), W, sigma);
%   
%   if params.WHITEN_ALBEDO_ENTROPY
%     dV = [dV1, dV2, dV3] * C;
%   else
%     dV = [dV1, dV2, dV3];
%   end
%   
%   losses_A.mult_entropy = mult * V;
%   d_loss_Av = d_loss_Av + mult * dV;
% 
%   loss_A = loss_A + losses_A.mult_entropy;
% 
% end
% 
% count = 0;
% for s = 1:params.APYR_DEPTH_ENTROPY{1}
%   
%   mask = data.valid_pyr{s};
%   mask_rep = repmat(mask, [1,1,3]);
%   
%   n = nnz(mask);
%   dA_band = d_loss_Av(count + [1:n],:);
%   count = count + n;
%   
%   d_loss_Apyr{s}(mask_rep) = d_loss_Apyr{s}(mask_rep) + dA_band(:);
%   
% end
% 
% 
% 
% 
% 
% Av = {};
% W = {};
% for s = 1:params.APYR_DEPTH_WHITE{1}
%   
%   Aband = Apyr{s};
%   mask = data.valid_pyr{s};
%   mask_rep = repmat(mask, [1,1,3]);
%   
%   Av{s} = reshape(Aband(mask_rep), [], 3);
%   W{s} = (4.^(s-1)) * ones(size(Av{s},1), 1);
% 
% end
% 
% Av = cat(1, Av{:});
% W  = cat(1, W{:});
% 
% d_loss_Av = zeros(size(Av));
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gamut{1};
% 
% if mult ~= 0
%   
%   hull_A = data.prior.albedo.RGB_hull.A;
%   hull_b = data.prior.albedo.RGB_hull.b;
%   
% %   rand('twister',5489)
% %   randn('state',0)
% % 
% %   hull_A = randn(size(data.prior.albedo.RGB_hull.A));
% %   hull_b = randn(size(data.prior.albedo.RGB_hull.b));
% 
% %   keyboard
%   
% %   hull_Ab = [hull_A, -hull_b];
% %   Av_aux = [Av'; ones(1, size(Av,1))];
% %   D = hull_Ab * Av_aux;
%   
% %   D = hull_A * Av';
% %   D = bsxfun(@minus, D, hull_b);
% %   [D_max, D_max_idx] = max(D, [], 1);
% %   D_max = max(0, D_max);
% %   D_max = D_max';
% %   D_max_idx = D_max_idx';
%   
%   [D_max, D_max_idx] = hullDist(Av', hull_A', hull_b);
%   
%   n = size(D_max,1);
%   
%   if params.GAMUT_POWER == 2
%     loss_A = loss_A + mult/n * ((D_max.^2) * W);
%     d_loss_Av = d_loss_Av + sparse(1:size(D_max,2), D_max_idx, ( (2*mult/n) * D_max .* W'), size(D,2), size(D,1)) * hull_A;
%   elseif params.GAMUT_POWER == 1
%     
%     losses_A.gamut = mult/n * sum(abs(D_max) .* W);
%     
%     above = D_max > 0;
%     d_loss_Av(above,:) = d_loss_Av(above,:) + hull_A(D_max_idx(above),:) .* repmat(( (mult/n) * sign(D_max(above)) .* W(above)), [1,3]);
% %     d_loss_Av = d_loss_Av + sparse(find(above), D_max_idx(above), ( (mult/n) * sign(D_max(above)) .* W(above)'), size(D,2), size(D,1)) * hull_A;
% %     d_loss_Av = d_loss_Av + sparse(1:size(D_max,2), D_max_idx, ( (mult/n) * sign(D_max) .* W'), size(D,2), size(D,1)) * hull_A;
%   end
%   
%   loss_A = loss_A + losses_A.gamut;
% 
% end
% 
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.white1{1}/32;
% 
% if mult ~= 0
%  
%   n = size(Av,1);
%   
%   mag_sq = sum(Av.^2 + 10^-4,2);
%   mag = sqrt(mag_sq);
%     
%   losses_A.white1 = mult/n * sum(mag .* W);
%   d_loss_Av = d_loss_Av + (mult/n) * Av .* repmat(W ./ mag, [1,3]);
% 
%   loss_A = loss_A + losses_A.white1;
% 
% end
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gray{1}/32;
% 
% if mult ~= 0
%   u = ones(3,1)/sqrt(3);
%   R = (eye(3,3) - u*u');
%   Av_resid = Av * R;
%   
%   n = size(Av,1);
% 
%   mag_sq = sum(Av_resid.^2 + 10^-4,2);
%   mag = sqrt(mag_sq);
%     
%   losses_A.gray = mult/n * sum(mag .* W);
%   d_loss_Av_resid = (mult/n) * Av .* repmat(W ./ mag, [1,3]);
%   d_loss_Av = d_loss_Av + d_loss_Av_resid * R;
% 
%   loss_A = loss_A + losses_A.gray;
% 
% end
% 
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gray2{1}/128;
% 
% if mult ~= 0
%   
%   y = -ones(3,1)/sqrt(3);
%   uv = [[1; -2; 1]/sqrt(6), [-1; 0; 1]/sqrt(2)];
%   
%   Y = max(0, Av* y);
%   UV = Av * uv;
%   
%   Y_soft = sqrt(Y.^2 + 0.01);
%   Y_soft_rep = repmat(Y_soft, [1,2]);
%   UV_norm = UV ./ Y_soft_rep;
%   
%   GSM = data.prior.albedo.gray.UV_GSM;
%   
%   W_unique = unique(W);
%   L = 0;
%   dL = nan(size(UV_norm));
%   for wi = 1:length(W_unique)
%     w = W_unique(wi);
%     idx = W == w;
%     [l, dl] = GSM_mvn_pdf(GSM, UV_norm(idx,:), 2);
%     L = L + w*l;
%     dL(idx,:) = w*dl;
%   end
%   
% %   % This is wrong, the weights are omitted
% %   [L, dL] = GSM_mvn_pdf(GSM, UV_norm, 2);
%   
%   losses_A.gray2 = mult/n * sum(L);
%   
%   dL = (mult/n) * dL;
%   d_loss_Av = d_loss_Av + (dL ./ Y_soft_rep) * uv';
%   d_loss_Av = d_loss_Av + sum((dL .* (repmat(-Y./(Y_soft.^3), [1,2]) .* UV)),2) * y';
% 
%   loss_A = loss_A + losses_A.gray2;
% 
% end
% 
% 
% 
% count = 0;
% for s = 1:params.APYR_DEPTH_WHITE{1}
%   
%   mask = data.valid_pyr{s};
%   mask_rep = repmat(mask, [1,1,3]);
%   
%   n = nnz(mask);
%   dA_band = d_loss_Av(count + [1:n],:);
%   count = count + n;
%   
%   d_loss_Apyr{s}(mask_rep) = d_loss_Apyr{s}(mask_rep) + dA_band(:);
%   
% end
% 
% 
% d_loss_A = {};
% for c = 1:3
%   d_loss_A{c} = reconLpyr_simple(cellfun(@(x) x(:,:,c), d_loss_Apyr, 'UniformOutput', false), [1;4;6;4;1]/16, 'zero');
% end
% d_loss_A = cat(3, d_loss_A{:});

















% figure(10); imagesc_norm(d_loss_A); imtight; drawnow;





% function [loss_A, d_loss_A, Apyr] = priorA(A, data, params)
% 
% MULT = params.MULT_OPTS.srfs;
% 
% GSM_RGB = data.prior.albedo.RGB.GSM;
% 
% Apyr = buildGpyr_color(A, params.APYR_DEPTH);
% 
% pind_mass = (params.APYR_ROOT{1}.^([1:length(Apyr)]'-1))/numel(A);
% 
% loss_A = 0;
% d_loss_Apyr = {};
% for s = 1:length(Apyr)
%   
%   d_loss_Apyr{s} = zeros(size(Apyr{s}));
%   
%   Aband = Apyr{s};
%   mask = data.valid_pyr{s};
% 
%   multMA = 3/25*MULT.albedo.MA{1} * pind_mass(s) / length(Apyr);
%   multMA3 = 3/25*MULT.albedo.MA3{1} * pind_mass(s) / length(Apyr);
%   multMYC = 3/25*MULT.albedo.MYC{1} * pind_mass(s) / length(Apyr);
%   multYC = 3*MULT.albedo.YC{1} * pind_mass(s) / length(Apyr);
%   
%   
%   if multMA3 ~= 0
%     
%     M = data.median_filter_mats{s};
%     MA = M * reshape(Aband, [], 3);
% 
%     [LL, dLL] = GSM_mvn_pdf(data.prior.albedo.MA.GSM_mvn{s}, MA, 2);
%     loss_A = loss_A + multMA3 * sum(LL(:));
%     dL_A = multMA3 * reshape(M' * dLL, size(Aband));
% 
%     d_loss_Apyr{s} = d_loss_Apyr{s} + dL_A;
%     
%   end
%   
%   
%   if multMA ~= 0
%     
%     M = data.median_filter_mats{s};
%     MA = M * reshape(Aband, [], 3);
% 
%     
%     if strcmp(params.PRIOR_PDF, 'GSM')
%       
%       W = data.prior.albedo.(params.MA_MODEL).whiten_map{s};
%       C = data.prior.albedo.(params.MA_MODEL).GSM{s};
%       
%       WMA = MA * W;
%       
%       epsilon = 10^-2;
%       MA_mag_sq = max(epsilon, sum(WMA.^2,2));
%       [L, dL] = interp1_fixed_sum_fast(MA_mag_sq, C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%       dL( MA_mag_sq <= epsilon ) = 0;
%       
%       loss_A = loss_A + multMA * sum(L(:));
%       dL_A = reshape(M' * ((repmat(2*multMA * dL, [1,3]) .* WMA) * W), size(Aband));
% 
%     elseif strcmp(params.PRIOR_PDF, 'myGaussian')
%       
%       W = data.prior.albedo.(params.MA_MODEL).whiten_map{s};
%       C = data.prior.albedo.(params.MA_MODEL).GSM{s};
%       
%       WMA = MA * W;
% 
%       MA_mag = sum(WMA.^2,2);
%       C = data.prior.albedo.MA.myGaussian{s};
%       [L, dL] = myGaussian_loss(MA_mag, C.alpha, C.beta, C.epsilon);
%       
%       loss_A = loss_A + multMA * sum(L(:));
%       dL_A = reshape(M' * ((repmat(2*multMA * dL, [1,3]) .* WMA) * W), size(Aband));
%       
%     end
%     
% %     MA_mag = sqrt(sum(WMA.^2,2) + 10^-4);
% %     [L, dL] = interp1_fixed_sum_fast(MA_mag, C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%     
% %     loss_A = loss_A + multMA * sum(L(:));
% %     dL_A = reshape(M' * ((repmat((multMA * dL) ./ MA_mag, [1,3]) .* WMA) * W), size(Aband));
%     
%     d_loss_Apyr{s} = d_loss_Apyr{s} + dL_A;
%     
%   end
%   
%   
%   
%   if multMYC ~= 0
%     
%     M = data.median_filter_mats{s};
%     
% %     [Mi, Mj, Mv] = find(M');
% %     Mi = reshape(Mi, 2, [])';
% %     MY2 = Y(Mi(:,1)) - Y(Mi(:,2));
%     
%     Y = mean(Aband,3);
%     MY = M * Y(:);
%     abs_MY = 3/4*abs(MY);
%     
%     C1 = (Aband(:,:,1) - Aband(:,:,2));
%     C2 = (Aband(:,:,2) - Aband(:,:,3));
%     MC1 = M * C1(:);
%     MC2 = M * C2(:);
%     MC = sqrt(MC1.^2 + MC2.^2 + 10^-4);
%     
%     C = data.prior.albedo.MYC.model1{s};
%     
%     [L, dL_MY, dL_MC] = interp2_fixed_sum_fast(abs_MY, MC, C.F_cost, C.edges{:});
%     
%     loss_A = loss_A + multMYC * sum(L(:));
%     
%     dL_Y = multMYC * reshape((dL_MY .* (3/4*sign(MY)))' * M, size(Y));
%     dL_C1 = multMYC * reshape(((MC1 ./ MC) .* dL_MC)' * M, size(Y));
%     dL_C2 = multMYC * reshape(((MC2 ./ MC) .* dL_MC)' * M, size(Y));
%    
%     dL_A = repmat(dL_Y/3, [1,1,3]);
%     dL_A(:,:,1) = dL_A(:,:,1) + dL_C1;
%     dL_A(:,:,2) = dL_A(:,:,2) - dL_C1 + dL_C2;
%     dL_A(:,:,3) = dL_A(:,:,3) - dL_C2;
%     
%     
%     d_loss_Apyr{s} = d_loss_Apyr{s} + dL_A;
% 
%     
%   end
%   
%   if multYC ~= 0
%         
%     [Yg, Cg, dYg, dCg] = getColorGrad(Aband);
%     
%     if params.ONE_ALBEDO_MODEL
%       C = spline_color{1};
%     else
%       C = spline_color{s};
%     end
%     [L, dL_Yg_mask, dL_Cg_mask] = interp2_fixed_sum_fast(Yg(mask), Cg(mask), C.F_cost, C.edges{:});    
%     
%     loss_A = loss_A + multYC * sum(L(:));
%     
%     dL_Yg = zeros(size(Yg));
%     dL_Yg(mask) = multYC * dL_Yg_mask;
%     
%     dL_Cg = zeros(size(Cg));
%     dL_Cg(mask) = multYC * dL_Cg_mask;
%     
%     dL_Aband = getColorGrad_backprop(dL_Yg, dL_Cg, dYg, dCg);
%     
%     d_loss_Apyr{s} = d_loss_Apyr{s} + dL_Aband;
%   
%   end
%     
%   
%   multRGB = MULT.albedo.RGB{1} * pind_mass(s) / length(Apyr);
%   
%   if multRGB ~= 0
%     
%     C = GSM_RGB{s};
%     for c = 1:3
%       
%       [G, dG] = getGradNorm(Aband(:,:,c));
%       
%       epsilon = 5*10^-3;
%       G = max(epsilon, G);
% 
%       [L, dL_mask] = interp1_fixed_sum_fast(G(mask), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%       
%       dL_mask(G(mask) <= epsilon) = 0;
%       
%       loss_A = loss_A + multRGB * sum(L(:));
%       
%       dL_G = zeros(size(G));
%       dL_G(mask) = multRGB * dL_mask;
%       
%       d_loss_Apyr{s}(:,:,c) = d_loss_Apyr{s}(:,:,c) + getGradNorm_backprop(dL_G, dG);
%       
%     end
%     
%       
%   end
%     
%   
%   mult = params.MULT_OPTS.srfs.albedo.entropy{1}/length(Apyr);
%   
%   if (mult ~= 0) && (s >= params.ALBEDO_ENTROPY_MINSCALE);
%     
%     mask_rep = repmat(mask, [1,1,3]);
%     Av = reshape(Aband(mask_rep), [], 3);
%     sigma = params.RENYI_SIGMA{1};
%     
%     if params.WHITEN_ALBEDO_ENTROPY
%       C = data.prior.albedo.A_whiten{1};
%       Av_white = Av * C;
%       
%       if params.USE_APPROX_RENYI
%         [V, dV_white] = renyiK_approx(Av_white, sigma, params.ENTROPY_PROJECTIONS_ORDER);
%       else
%         [V, dV1, dV2, dV3] = renyi3(Av_white(:,1), Av_white(:,2), Av_white(:,3), sigma);
%         dV_white = [dV1, dV2, dV3];
%       end
%       
%       dV = dV_white * C;
%     else
%       
%       if params.USE_APPROX_RENYI
%         [V, dV] = renyiK_approx(Av, sigma, params.ENTROPY_PROJECTIONS_ORDER);
%       else
%         [V, dV1, dV2, dV3] = renyi3(Av(:,1), Av(:,2), Av(:,3), sigma);
%         dV = [dV1, dV2, dV3];
%       end
%       
%     end
%         
%     loss_A = loss_A + mult*V;
%     d_loss_Apyr{s}(mask_rep) = d_loss_Apyr{s}(mask_rep) + mult*dV(:);
%     
%   end
%   
% end
% 
% 
% 
% 
% Av = {};
% W = {};
% for s = 1:length(Apyr)
%   
%   Aband = Apyr{s};
%   mask = data.valid_pyr{s};
%   mask_rep = repmat(mask, [1,1,3]);
%   
%   Av{s} = reshape(Aband(mask_rep), [], 3);
%   W{s} = (4.^(s-1)) * ones(size(Av{s},1), 1);
% 
% end
% 
% Av = cat(1, Av{:});
% W  = cat(1, W{:});
% 
% d_loss_Av = zeros(size(Av));
% 
% mult = params.MULT_OPTS.srfs.albedo.multi_entropy{1};
% 
% if mult ~= 0
%   
%   sigma = params.RENYI_SIGMA{1};
%   C = data.prior.albedo.A_whiten;
%   
%   if params.WHITEN_ALBEDO_ENTROPY
%     Av_white = Av * C;
%   else
%     Av_white = Av;
%   end
%   
% %   Av_white = single(Av_white);
% %   W = single(W);
%   
%   [V, dV1, dV2, dV3] = renyi3W(Av_white(:,1), Av_white(:,2), Av_white(:,3), W, sigma);
%   
%   if params.WHITEN_ALBEDO_ENTROPY
%     dV = [dV1, dV2, dV3] * C;
%   else
%     dV = [dV1, dV2, dV3];
%   end
%   
% %   V = double(V);
% %   dV = double(dV);
%   
%   loss_A = loss_A + mult * V;
%   d_loss_Av = d_loss_Av + mult * dV;
%   
% end
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gamut{1};
% 
% if mult ~= 0
%   
%   hull_A = data.prior.albedo.RGB_hull.A;
%   hull_b = data.prior.albedo.RGB_hull.b;
%   
% %   rand('twister',5489)
% %   randn('state',0)
% % 
% %   hull_A = randn(size(data.prior.albedo.RGB_hull.A));
% %   hull_b = randn(size(data.prior.albedo.RGB_hull.b));
% 
% %   keyboard
%   
% %   hull_Ab = [hull_A, -hull_b];
% %   Av_aux = [Av'; ones(1, size(Av,1))];
% %   D = hull_Ab * Av_aux;
%   
% %   D = hull_A * Av';
% %   D = bsxfun(@minus, D, hull_b);
% %   [D_max, D_max_idx] = max(D, [], 1);
% %   D_max = max(0, D_max);
% %   D_max = D_max';
% %   D_max_idx = D_max_idx';
%   
%   [D_max, D_max_idx] = hullDist(Av', hull_A', hull_b);
%   
%   n = size(D_max,1);
%   
%   if params.GAMUT_POWER == 2
%     loss_A = loss_A + mult/n * ((D_max.^2) * W);
%     d_loss_Av = d_loss_Av + sparse(1:size(D_max,2), D_max_idx, ( (2*mult/n) * D_max .* W'), size(D,2), size(D,1)) * hull_A;
%   elseif params.GAMUT_POWER == 1
%     
%     loss_A = loss_A + mult/n * sum(abs(D_max) .* W);
%     
%     above = D_max > 0;
%     d_loss_Av(above,:) = d_loss_Av(above,:) + hull_A(D_max_idx(above),:) .* repmat(( (mult/n) * sign(D_max(above)) .* W(above)), [1,3]);
% %     d_loss_Av = d_loss_Av + sparse(find(above), D_max_idx(above), ( (mult/n) * sign(D_max(above)) .* W(above)'), size(D,2), size(D,1)) * hull_A;
% %     d_loss_Av = d_loss_Av + sparse(1:size(D_max,2), D_max_idx, ( (mult/n) * sign(D_max) .* W'), size(D,2), size(D,1)) * hull_A;
%   end
%   
% end
% 
% mult = params.MULT_OPTS.srfs.albedo.white1{1}/32;
% 
% if mult ~= 0
%  
%   n = size(Av,1);
%   
%   mag_sq = sum(Av.^2 + 10^-4,2);
%   mag = sqrt(mag_sq);
%     
%   loss_A = loss_A + mult/n * sum(mag .* W);
%   d_loss_Av = d_loss_Av + (mult/n) * Av .* repmat(W ./ mag, [1,3]);
%   
% end
% 
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gray{1}/32;
% 
% if mult ~= 0
%   u = ones(3,1)/sqrt(3);
%   R = (eye(3,3) - u*u');
%   Av_resid = Av * R;
%   
%   n = size(Av,1);
% 
%   mag_sq = sum(Av_resid.^2 + 10^-4,2);
%   mag = sqrt(mag_sq);
%     
%   loss_A = loss_A + mult/n * sum(mag .* W);
%   d_loss_Av_resid = (mult/n) * Av .* repmat(W ./ mag, [1,3]);
%   d_loss_Av = d_loss_Av + d_loss_Av_resid * R;
% 
% end
% 
% 
% mult = params.MULT_OPTS.srfs.albedo.gray2{1}/128;
% 
% if mult ~= 0
%   
%   y = -ones(3,1)/sqrt(3);
%   uv = [[1; -2; 1]/sqrt(6), [-1; 0; 1]/sqrt(2)];
%   
%   Y = max(0, Av* y);
%   UV = Av * uv;
%   
%   Y_soft = sqrt(Y.^2 + 0.01);
%   Y_soft_rep = repmat(Y_soft, [1,2]);
%   UV_norm = UV ./ Y_soft_rep;
%   
%   GSM = data.prior.albedo.gray.UV_GSM;
%   
%   W_unique = unique(W);
%   L = 0;
%   dL = nan(size(UV_norm));
%   for wi = 1:length(W_unique)
%     w = W_unique(wi);
%     idx = W == w;
%     [l, dl] = GSM_mvn_pdf(GSM, UV_norm(idx,:), 2);
%     L = L + w*l;
%     dL(idx,:) = w*dl;
%   end
%   
% %   % This is wrong, the weights are omitted
% %   [L, dL] = GSM_mvn_pdf(GSM, UV_norm, 2);
%   
%   loss_A = loss_A + mult/n * sum(L);
%   
%   dL = (mult/n) * dL;
%   d_loss_Av = d_loss_Av + (dL ./ Y_soft_rep) * uv';
%   d_loss_Av = d_loss_Av + sum((dL .* (repmat(-Y./(Y_soft.^3), [1,2]) .* UV)),2) * y';
% 
% end
% 
% 
% % mult = params.MULT_OPTS.srfs.albedo.white2{1};
% % 
% % if mult ~= 0
% %  
% %   n = size(Av,1);
% % 
% %   loss_A = loss_A + mult/n * 0.5*sum(W' * (Av.^2));
% %   d_loss_Av = d_loss_Av + (mult/n) * repmat(W, [1,3]) .* Av;
% %   
% % end
% % 
% % 
% % mult = params.MULT_OPTS.srfs.albedo.gaussian{1};
% % 
% % if mult ~= 0
% %  
% %   n = size(Av,1);
% % 
% %   [LL, dLL] = lmvnpdf2(Av, data.prior.albedo.RGB_gaussian.mu, data.prior.albedo.RGB_gaussian.Sigma);
% %   LL = LL .* W;
% %   dLL = dLL .* repmat(W, [1,3]);
% %   
% %   loss_A = loss_A - mult/n * sum(LL);
% %   d_loss_Av = d_loss_Av - (mult/n) * dLL;
% %   
% % end
% 
% 
% 
% count = 0;
% for s = 1:length(Apyr)
%   
%   mask = data.valid_pyr{s};
%   mask_rep = repmat(mask, [1,1,3]);
%   
%   n = nnz(mask);
%   dA_band = d_loss_Av(count + [1:n],:);
%   count = count + n;
%   
%   d_loss_Apyr{s}(mask_rep) = d_loss_Apyr{s}(mask_rep) + dA_band(:);
%   
% end
% 
% 
% 
% d_loss_A = {};
% for c = 1:3
%   d_loss_A{c} = reconLpyr_simple(cellfun(@(x) x(:,:,c), d_loss_Apyr, 'UniformOutput', false), [1;4;6;4;1]/16, 'zero');
% end
% d_loss_A = cat(3, d_loss_A{:});
% 
% 
% 
