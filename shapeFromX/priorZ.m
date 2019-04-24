function [loss_Z, d_loss_Z, losses_Z] = priorZ(Z, data, params, N, dN_Z)

d_loss_Z = 0;

mult_smooth = params.MULT_OPTS.saifs.height.smooth{1};
mult_mod_smooth = (1/((params.Z_MEDIAN_HALFWIDTH*2+1)^2))/numel(Z);

mult_mod_contour = 60 / numel(Z);
mult_contour = params.MULT_OPTS.saifs.height.contour{1};

mult_mod_slant = 1 / numel(Z);
mult_slant = params.MULT_OPTS.saifs.height.slant{1};


if mult_smooth ~= 0

  [KZ, dKZ] = getK_fast(Z);
  C = data.prior.height.MKZ.GSM;
  
  M = data.Z_median_filter_mat;
  M_T = data.Z_median_filter_mat_T;
  
  MKZ = M * KZ(:);
  
  abs_MKZ = abs(MKZ);
  sign_MKZ = sign(MKZ);
  
  epsilon = params.Z_SMOOTH_EPSILON{1};
  if epsilon > 0
    soft_MKZ = sqrt(abs_MKZ.^2 + epsilon.^2);
  else
    soft_MKZ = abs_MKZ;
  end
  
%   epsilon = 10^-4;
%   abs_MKZ = max(epsilon, abs_MKZ);
  [L, dL_mask] = interp1_fixed_sum_fast(soft_MKZ(:), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%   dL_mask(abs_MKZ <= epsilon) = 0;
  
  if epsilon > 0
    dL_mask = dL_mask .* (abs_MKZ ./ max(eps,soft_MKZ));
  end

  dL_mask = dL_mask .* sign_MKZ;
  
  dL_KZ = reshape(M_T * dL_mask, size(KZ));%reshape(dL_mask' * M, size(KZ));
  
  losses_Z.smooth = mult_mod_smooth * sum(L(:));
  d_loss_Z = d_loss_Z + mult_smooth * mult_mod_smooth * getK_backprop_fast_mat(dL_KZ, dKZ);
  
%   if params.USE_KZUP
%     d_loss_Z = d_loss_Z + mult_smooth * mult_mod_smooth * getKup_backprop(dL_KZ, dKZ);
%   else
%     d_loss_Z = d_loss_Z + mult_smooth * mult_mod_smooth * getK_backprop(dL_KZ, dKZ);
%   end
  
else
  losses_Z.smooth = 0;
end



if nargin < 4
  [N, dN_Z] = getNormals_conv(Z);
end
N1 = N(:,:,1);
N2 = N(:,:,2);



mask = data.border.mask;
normal = data.border.normal;

d1 = N1(mask) - normal(:,1);
d2 = N2(mask) - normal(:,2);

d_loss_N1 = zeros(size(N1));
d_loss_N2 = zeros(size(N2));

dmag = sqrt(d1.^2 + d2.^2);
losses_Z.contour = mult_mod_contour * sum(dmag);

d_loss_N1(mask) = mult_contour * mult_mod_contour * (d1 ./ dmag);
d_loss_N2(mask) = mult_contour * mult_mod_contour * (d2 ./ dmag);


N3 = N(:,:,3);
mask = data.valid;

losses_Z.slant = mult_mod_slant * sum(-log(N3(mask))); % Actual PDF is sum(-log(2*N3(mask)))
d_loss_N3 = mask .* (-mult_mod_slant * mult_slant)./max(eps,N3);

d_loss_D1 = d_loss_N1 .* dN_Z.F1_1 + d_loss_N2 .* dN_Z.F2_1 + d_loss_N3 .* dN_Z.F1_3;
d_loss_D2 = d_loss_N1 .* dN_Z.F1_2 + d_loss_N2 .* dN_Z.F2_2 + d_loss_N3 .* dN_Z.F2_3;

d_loss_Z = d_loss_Z + conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);



loss_Z = mult_smooth * losses_Z.smooth + mult_contour * losses_Z.contour + mult_slant * losses_Z.slant;







% function [loss_Z, d_loss_Z, losses_Z] = priorZ(Z, data, params, N, dN_Z)
% 
% [KZpyr, dKZ] = buildKpyr(Z, params.ZPYR_DEPTH{1});
% 
% loss_Z = 0;
% d_loss_Zpyr = {};
% 
% % mult = (1/25) * params.MULT_OPTS.srfs.height.L1L2_MKZ{1} / numel(Z);
% % 
% % if mult ~= 0
% % 
% %   epsilon = 0.1;
% % 
% %   KZ = KZpyr{1};
% %   Mi = data.median_filter_pairs{1}(:,1);
% %   Mj = data.median_filter_pairs{1}(:,2);
% %   
% %   x = KZ(Mi);
% %   y = KZ(Mj);
% %   
% %   del = x - y;
% %   abs_del = abs(del);
% %   mag_sq = x.^2 + y.^2 + epsilon;
% %   mag = sqrt(mag_sq);
% %   mag3 = mag .* mag_sq;
% %   dx = -(x .* abs_del) ./ mag3 + sign(del)./mag;
% %   dy = -(y .* abs_del) ./ mag3 - sign(del)./mag;
% %   
% %   dL_KZ = reshape(accumarray(Mi, dx, [numel(KZ),1]) + accumarray(Mj, dy, [numel(KZ),1]), size(KZ));
% %   
% %   L = sum(abs_del ./ mag);
% %   
% % %   del_KZ = M * KZ(:);
% % %   abs_del_KZ = abs(del_KZ);
% % %   sum_sq_KZ = abs(M) * (KZ(:).^2) + epsilon.^2;
% % %   mag_KZ = sqrt(sum_sq_KZ);
% % %   
% % %   L = sum(abs_del_KZ ./ mag_KZ);
% % %   dL_numer = sign(del_KZ) ./ mag_KZ;
% % %   dL_denom = -(del_KZ) ./ mag_KZ;
% % %   
% % %   dL_KZ = reshape(dL_mask * M, size(KZ));
% %   
% % %   mult * sum(L(:))
% %   losses_Z.L1L2_MKZ = mult * sum(L(:));
% %   d_loss_Zpyr{1} = mult * getK_backprop(dL_KZ, dKZ{1});
% % 
% %   loss_Z = loss_Z + losses_Z.L1L2_MKZ;
% % end
% 
% % mult = (1/25)*params.MULT_OPTS.srfs.height.MlogKZ{1}/numel(Z);
% % 
% % if mult ~= 0
% % 
% %   KZ = KZpyr{1};
% %   [logKZ, d_logKZ] = getLogK(KZ);
% %   
% %   C = data.prior.height.MlogKZ.GSM{1};
% %   M = data.median_filter_mats{1};
% %   
% %   MKZ = M * logKZ(:);
% %   
% %   abs_MKZ = abs(MKZ);
% %   sign_MKZ = sign(MKZ);
% %   
% %   epsilon = 10^-4;
% %   abs_MKZ = max(epsilon, abs_MKZ);
% %   [L, dL_mask] = interp1_fixed_sum_fast(abs_MKZ(:), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
% %   dL_mask(abs_MKZ <= epsilon) = 0;
% %   
% %   dL_mask = dL_mask .* sign_MKZ;
% %   
% %   dL_logKZ = reshape(dL_mask' * M, size(KZ));
% %   dL_KZ = dL_logKZ .* d_logKZ;
% %   
% %   losses_Z.MlogKZ = mult * sum(L(:));
% %   d_loss_Zpyr{1} = mult * getK_backprop(dL_KZ, dKZ{1});
% % 
% %   loss_Z = loss_Z + losses_Z.MlogKZ;
% % end
% 
% 
% 
% 
% losses_Z.MKZ = 0;
% 
% for s = 1:length(KZpyr)
% 
%   d_loss_Zpyr{s} = zeros(size(KZpyr{s}));
% 
%   KZ = KZpyr{s};
%   
%   mult = (1/25)*params.MULT_OPTS.srfs.height.MKZ{1}/numel(Z) * (params.ZPYR_ROOT{1}.^(s-1)) / length(KZpyr);
%   
%   if mult ~= 0
% 
%     C = data.prior.height.MKZ.GSM{s};
%     M = data.median_filter_mats{s};
%     
%     MKZ = M * KZ(:);
%     
%     abs_MKZ = abs(MKZ);
%     sign_MKZ = sign(MKZ);
% 
%     epsilon = 10^-4;
%     abs_MKZ = max(epsilon, abs_MKZ);
%     [L, dL_mask] = interp1_fixed_sum_fast(abs_MKZ(:), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%     dL_mask(abs_MKZ <= epsilon) = 0;
%     
%     dL_mask = dL_mask .* sign_MKZ;
%     
%     dL_KZ = reshape(dL_mask' * M, size(KZ));
%     
%     losses_Z.MKZ = losses_Z.MKZ + mult * sum(L(:));
%     d_loss_Zpyr{s} = mult * getK_backprop(dL_KZ, dKZ{s});
%     
%   end
%   
%   d_loss_Z = reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% 
% end
% 
% loss_Z = loss_Z + losses_Z.MKZ;
% 
% 
% if nargin < 4
%   [N, dN_Z] = getNormals_conv(Z);
% end
% N1 = N(:,:,1);
% N2 = N(:,:,2);
% 
% 
% mult = params.MULT_OPTS.srfs.height.contour{1}*60 / numel(Z);
% 
% mask = data.border.mask;
% normal = data.border.normal;
% 
% d1 = N1(mask) - normal(:,1);
% d2 = N2(mask) - normal(:,2);
% 
% d_loss_N1 = zeros(size(N1));
% d_loss_N2 = zeros(size(N2));
% 
% dmag = sqrt(d1.^2 + d2.^2);
% losses_Z.contour = mult * sum(dmag);
% loss_Z = loss_Z + losses_Z.contour;
% 
% d_loss_N1(mask) = mult * (d1 ./ dmag);
% d_loss_N2(mask) = mult * (d2 ./ dmag);
% 
% 
% 
% mult = params.MULT_OPTS.srfs.height.slant{1} / numel(Z);
% 
% N3 = N(:,:,3);
% mask = data.valid;
% 
% losses_Z.slant = mult * sum(-log(N3(mask))); % Actual PDF is sum(-log(2*N3(mask)))
% loss_Z = loss_Z + losses_Z.slant;
% d_loss_N3 = mask .* (-mult)./max(eps,N3);
% 
% 
% d_loss_D1 = d_loss_N1 .* dN_Z.F1_1 + d_loss_N2 .* dN_Z.F2_1 + d_loss_N3 .* dN_Z.F1_3;
% d_loss_D2 = d_loss_N1 .* dN_Z.F1_2 + d_loss_N2 .* dN_Z.F2_2 + d_loss_N3 .* dN_Z.F2_3;
% 
% d_loss_Z = d_loss_Z + conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);
% 








% function [loss_Z, d_loss_Z] = priorZ(Z, data, params, N, dN_Z)
% 
% MULT = params.MULT_OPTS.srfs;
% 
% if strcmp(params.CURVATURE_MODE, 'K')
% 
%   GSM_GKZ = data.prior.height.GKZ.(params.PRIOR_PDF);
%   [KZpyr, dKZ, Zpyr] = buildKpyr(Z, params.ZPYR_DEPTH);
%   
% elseif strcmp(params.CURVATURE_MODE, 'IK')
% 
%   GSM_GKZ = data.prior.height.GIKZ.(params.PRIOR_PDF);
%   [KZpyr, dKZ, Zpyr] = buildIKpyr(Z, params.ZPYR_DEPTH);
%   
% end
% 
% pind_mass = params.ZPYR_ROOT{1}.^([1:length(KZpyr)]'-1);
% 
% mult_height_contour = MULT.height.contour{1}*100;
% mult_GKZ = MULT.height.GKZ{1};
% 
% loss_Z = 0;
% d_loss_Z = zeros(size(Z));
% d_loss_Zpyr = {};
% 
%   
% for s = 1:length(KZpyr)
%   
%   mask = data.valid_pyr{s};
%   
%   mult = mult_GKZ * pind_mass(s);
%   
%   if mult ~= 0
% 
%     KZ = KZpyr{s};
%     
% %     [GKZ, dGKZ] = getGrad2(KZ);
% %     mask = mask(1:end-1,1:end-1);
% %     
% %     epsilon = 5*10^-3;
% %     GKZ = max(epsilon, GKZ);
% %     
% %     if strcmp(params.PRIOR_PDF, 'myGaussian')
% %       [LL, dLL_mask] = myGaussian_pdf(GKZ(mask), GSM_GKZ{s}.alpha, GSM_GKZ{s}.beta, GSM_GKZ{s}.epsilon);
% %       L = -LL;
% %       dL_mask = -dLL_mask;
% %     else
% %       [L, dL_mask] = interp1_fixed_sum_fast(GKZ(mask), GSM_GKZ{s}.lut.bin_range(1), GSM_GKZ{s}.lut.bin_width, GSM_GKZ{s}.lut.F_cost);
% %     end
% %     
% %     dL_mask(GKZ(mask) <= epsilon) = 0;
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_GKZ = zeros(size(GKZ));
% %     dL_GKZ(mask) = mult * dL_mask;
% %     
% %     if strcmp(params.CURVATURE_MODE, 'K')
% %       d_loss_Zpyr{s} = getK_backprop(getGrad2_backprop(dL_GKZ, dGKZ), dKZ{s});
% %     elseif strcmp(params.CURVATURE_MODE, 'IK')
% %       d_loss_Zpyr{s} = getIK_backprop(getGrad2_backprop(dL_GKZ, dGKZ), dKZ{s});
% %     end
%     
%     
%     [GKZ, dGKZ] = getGradNorm(KZ);
%     
%     epsilon = 5*10^-3;
%     GKZ = max(epsilon, GKZ);
% 
%     if strcmp(params.PRIOR_PDF, 'myGaussian')
%       [LL, dLL_mask] = myGaussian_pdf(GKZ(mask), GSM_GKZ{s}.alpha, GSM_GKZ{s}.beta, GSM_GKZ{s}.epsilon);
%       L = -LL;
%       dL_mask = -dLL_mask;
%     else
%       [L, dL_mask] = interp1_fixed_sum_fast(GKZ(mask), GSM_GKZ{s}.lut.bin_range(1), GSM_GKZ{s}.lut.bin_width, GSM_GKZ{s}.lut.F_cost);
%     end
%     
%     dL_mask(GKZ(mask) <= epsilon) = 0;
%     
%     loss_Z = loss_Z + mult * sum(L(:));
%     dL_GKZ = zeros(size(KZ));
%     dL_GKZ(mask) = mult * dL_mask;
%     
%     if strcmp(params.CURVATURE_MODE, 'K')
%       d_loss_Zpyr{s} = getK_backprop(getGradNorm_backprop(dL_GKZ, dGKZ), dKZ{s});
%     elseif strcmp(params.CURVATURE_MODE, 'IK')
%       d_loss_Zpyr{s} = getIK_backprop(getGradNorm_backprop(dL_GKZ, dGKZ), dKZ{s});
%     end
%     
% 
%   else
%     d_loss_Zpyr{s} = zeros(size(mask));
%   end
% end
% 
% 
% mult = MULT.height.lowflat{1};
% 
% if mult ~= 0
%   s = params.ZPYR_DEPTH;
%   
%   delta = Zpyr{s};
%   mask = data.valid_pyr{s};
%   
%   loss_Z = loss_Z + mult*0.5*sum(delta(mask).^2);
%   d_loss_Zpyr{s} = d_loss_Zpyr{s} + mult*(mask.*delta);
%   
% end
% 
% 
% d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% 
% 
% 
% mult = mult_height_contour;
% 
% if nargin < 4
%   [N, dN_Z] = getNormals_conv(Z);
% end
% N1 = N(:,:,1);
% N2 = N(:,:,2);
% 
% 
% mask = data.border.mask;
% normal = data.border.normal;
% 
% d1 = N1(mask) - normal(:,1);
% d2 = N2(mask) - normal(:,2);
% 
% d_loss_N1 = zeros(size(N1));
% d_loss_N2 = zeros(size(N2));
% 
% dmag = sqrt(d1.^2 + d2.^2);
% loss_Z = loss_Z + mult * sum(dmag);
% d_loss_N1(mask) = mult * (d1 ./ dmag);
% d_loss_N2(mask) = mult * (d2 ./ dmag);
% 
% 
% 
% mult = params.MULT_OPTS.srfs.height.slant{1};
% 
% if mult ~= 0
%     
%   N3 = N(:,:,3);
%   mask = data.valid;
%   
% %   mask = true(size(mask));
% %   mask([1,end],:) = false;
% %   mask(:,[1,end]) = false;
%   
% %   loss_Z = loss_Z + mult * sum(-log(2*N3(mask)));  % Actual PDF
%   loss_Z = loss_Z + mult * sum(-log(N3(mask))); % Don't need the factor of 2.
%   d_loss_N3 = mask .* (-mult)./max(eps,N3);
%   
% else
%   d_loss_N3 = 0;
% end
% 
% 
% d_loss_D1 = d_loss_N1 .* dN_Z.F1_1 + d_loss_N2 .* dN_Z.F2_1 + d_loss_N3 .* dN_Z.F1_3;
% d_loss_D2 = d_loss_N1 .* dN_Z.F1_2 + d_loss_N2 .* dN_Z.F2_2 + d_loss_N3 .* dN_Z.F2_3;
% 
% d_loss_Z = d_loss_Z + conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);
% 
% 
% 
% 
% mult = params.MULT_OPTS.srfs.height.multislant{1};
% 
% if mult ~= 0
%   
%   d_loss_Zpyr = {};
%   for s = 1:length(Zpyr)
%   
%     mask = data.valid_pyr{s};
%     
% %     mask = true(size(mask));
% %     mask([1,end],:) = false;
% %     mask(:,[1,end]) = false;
%     
%     mult_band = mult * pind_mass(s);
%     Zband = Zpyr{s};
% 
%     if (s >= 2)
%       [N, dN_Z] = getNormals_conv(Zband);
%     end
%     N3 = N(:,:,3);
% 
%     loss_Z = loss_Z + mult_band * sum(-log(2*N3(mask)));
%     dL = mask .* (-mult_band)./max(eps,N3);
%     
%     d_loss_Zpyr{s} = conv3d( dL .* dN_Z.F1_3, dN_Z.f1) + conv3d( dL .* dN_Z.F2_3, dN_Z.f2);
%     
%   end
%   
%   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
%   
% end
% 
% 
% 
% 
% mult = params.MULT_OPTS.srfs.height.frontness{1};
% 
% if mult ~= 0
% 
%   [GZ, dGZ] = getGradNorm(Z);
%     
%   mask = data.valid;
%     
%   GSM = data.prior.height.GZ.(params.PRIOR_PDF){1};
%   
%   [L, dL_mask] = interp1_fixed_sum_fast(GZ(mask), GSM.lut.bin_range(1), GSM.lut.bin_width, GSM.lut.F_cost);
%   dL = zeros(size(Z));
%   dL(mask) = dL_mask;
%   
%   loss_Z = loss_Z + mult*L;
%   d_loss_Z = d_loss_Z + getGradNorm_backprop(mult*dL, dGZ);
%   
% end
% 
% 
% mult_GZ = params.MULT_OPTS.srfs.height.multifront{1};
% 
% if mult_GZ ~= 0
%   GSM_GZ = data.prior.height.GZ.(params.PRIOR_PDF);
%   
%   d_loss_Zpyr = {};
% 
%   for s = 1:length(Zpyr)
%     
%     mask = data.valid_pyr{s};
%     
%     mult = mult_GZ * pind_mass(s);
%     Zband = Zpyr{s};
%     
%     [GZ, dGZ] = getGradNorm(Zband);
%     
%     epsilon = 5*10^-3;
%     GZ = max(epsilon, GZ);
%     
%     [L, dL_mask] = interp1_fixed_sum_fast(GZ(mask), GSM_GZ{s}.lut.bin_range(1), GSM_GZ{s}.lut.bin_width, GSM_GZ{s}.lut.F_cost);
%     
%     dL_mask(GZ(mask) <= epsilon) = 0;
%     
%     loss_Z = loss_Z + mult * sum(L(:));
%     dL_GZ = zeros(size(GZ));
%     dL_GZ(mask) = mult * dL_mask;
%     
%     d_loss_Zpyr{s} = getGradNorm_backprop(dL_GZ, dGZ);
%   
%   end
%   
%   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% 
% end
% 
% 
% mult = params.MULT_OPTS.srfs.height.flat{1}/20000;
% if mult ~= 0
%   
%   delta = Z;
%   
%   mask = data.valid;
%   loss_Z = loss_Z + mult*0.5*sum(delta(mask).^2);
%   d_loss_Z = d_loss_Z + mult*(mask.*delta);
%   
% end
% 
% loss_Z = loss_Z / numel(Z);
% d_loss_Z = d_loss_Z / numel(Z);
% 






% function [loss_Z, d_loss_Z] = priorZ(Z, data, params)
% 
% MULT = params.MULT_OPTS.srfs;
% 
% loss_Z = 0;
% d_loss_Z = zeros(size(Z));
% 
% mult_GKZ = MULT.height.GKZ{1};
% mult_GIKZ = MULT.height.GIKZ{1};
% mult_MKZ = MULT.height.MKZ{1};
% mult_MIKZ = MULT.height.MIKZ{1};
% mult_NKZ = MULT.height.NKZ{1};
% mult_NIKZ = MULT.height.NIKZ{1};
% mult_KZ2 = MULT.height.KZ2{1};
% mult_dN = MULT.height.dN{1};
% 
% epsilon = 10^-4;
% 
% 
% % if mult_dN ~= 0
% %   GSM_dN = data.prior.height.dN.(params.PRIOR_PDF);
% %   
% %   Zpyr = buildGpyr_simple(Z, params.ZPYR_DEPTH, [1;4;6;4;1]/(16*sqrt(2)), 'repeat');
% %   Npyr = {};
% %   dNpyr = {};
% %   for s = 1:length(Zpyr)
% %     [Npyr{s}, dNpyr{s}] = getNormals_conv(Zpyr{s});
% %   end
% %   
% %   pind_mass = params.ZPYR_ROOT{1}.^([1:length(Zpyr)]'-1);
% %   pind_mass = pind_mass / ((params.Z_MEDIAN_HALFWIDTH*2+1)^2);
% %   
% %   d_loss_Zpyr = {};
% %   
% %   for s = 1:length(Zpyr)
% %     
% %     M = data.median_filter_mats_weighted{s};
% %     
% %     C = GSM_dN{s};
% %     
% %     mult = mult_dN * pind_mass(s) / length(Zpyr);
% %     
% %     N = Npyr{s};
% %     Nv = reshape(N, [], 3);
% %     
% %     MN = M * Nv;
% %     dN = sqrt(sum(MN.^2,2) + 10^-5);
% %     
% %     [L, dL] = interp1_fixed_sum_fast(dN, C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
% %     
% %     dL_N = reshape(((repmat((dL./dN), [1,3]) .* MN)' * M)', size(N));
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_N = mult * dL_N;
% %     
% %     dN_Z = dNpyr{s};
% %     dL_F1 = dN_Z.F1_1 .* dL_N(:,:,1) + dN_Z.F1_2 .* dL_N(:,:,2) + dN_Z.F1_3 .* dL_N(:,:,3);
% %     dL_F2 = dN_Z.F2_1 .* dL_N(:,:,1) + dN_Z.F2_2 .* dL_N(:,:,2) + dN_Z.F2_3 .* dL_N(:,:,3);
% %     dL_Zband = conv3d( dL_F1, dN_Z.f1) + conv3d( dL_F2, dN_Z.f2);
% % 
% %     d_loss_Zpyr{s} = dL_Zband;
% %     
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% % 
% % 
% % if mult_KZ2 ~= 0
% %   
% %   Cs = data.prior.height.KZ2.model1;
% %   [KZpyr, dKZ] = buildKpyr(Z, params.ZPYR_DEPTH);
% %   
% %   pind_mass = params.ZPYR_ROOT{1}.^([1:length(KZpyr)]'-1);
% %   
% %   d_loss_Zpyr = {};
% %   
% %   for s = 1:length(KZpyr)
% %     
% %     mask = data.valid_pyr{s};
% %   
% %     mult = mult_KZ2 * pind_mass(s) / length(KZpyr);
% %     
% %     C = Cs{s};
% %     
% %     KZ = KZpyr{s};
% %     [GKZ, dGKZ] = getGradNorm(KZ);
% %     KZb = blur3(KZ);
% %     
% %     abs_KZb = 3*abs(KZb);
% %     sqrt_GKZ = 2*sqrt(GKZ + 10^-4);
% % 
% %     
% %     [L, dL1_mask, dL2_mask] = interp2_fixed_sum_fast(abs_KZb(mask), sqrt_GKZ(mask), C.F_cost, C.edges{:});    
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     
% %     dL1 = zeros(size(mask));
% %     dL2 = zeros(size(mask));
% %     dL1(mask) = mult * dL1_mask;
% %     dL2(mask) = mult * dL2_mask;
% %     
% %     dL_KZb = 3 * dL1 .* sign(KZb);
% %     dL_GKZ = 2*dL2./sqrt_GKZ;
% %     
% %     dL_KZ = blur3(dL_KZb) + getGradNorm_backprop(dL_GKZ, dGKZ);
% %     
% %     d_loss_Zpyr{s} = getK_backprop(dL_KZ, dKZ{s});
% %     
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% % 
% % 
% % 
% % if mult_NIKZ ~= 0
% %   
% %   GSM_NIKZ = data.prior.height.NIKZ.(params.PRIOR_PDF);
% %   [IKZpyr, dIKZ] = buildIKpyr(Z, params.ZPYR_DEPTH);
% %   
% %   pind_mass = params.ZPYR_ROOT{1}.^([1:length(IKZpyr)]'-1);
% %   
% %   d_loss_Zpyr = {};
% %   
% %   for s = 1:length(IKZpyr)
% %     
% %     mask = data.valid_pyr{s};
% %   
% %     mult = mult_NIKZ * pind_mass(s) / length(IKZpyr);
% %     
% %     C = GSM_NIKZ{s};
% %     
% %     IKZ = IKZpyr{s};
% %     
% %     [GIKZ, dGIKZ] = getGradNorm(IKZ);
% %     
% %     IKZb = blur3(IKZ);
% %     
% %     IKZ_pad = sqrt(IKZb.^2 + 10^-4);
% %     NIKZ = GIKZ ./ IKZ_pad;
% %     
% %     NIKZ = max(epsilon, NIKZ);
% %     
% %     [L, dL_mask] = interp1_fixed_sum_fast(NIKZ(mask), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
% %     
% %     dL_mask(NIKZ(mask) <= epsilon) = 0;
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_NIKZ = zeros(size(IKZ));
% %     dL_NIKZ(mask) = mult * dL_mask;
% %     
% %     dL_GIKZ = dL_NIKZ ./ IKZ_pad;
% %     dL_IKZb = -(GIKZ .* IKZb .* dL_NIKZ) ./ (IKZ_pad.^3);
% %     dL_IKZ = blur3(dL_IKZb) + getGradNorm_backprop(dL_GIKZ, dGIKZ);
% %     
% %     d_loss_Zpyr{s} = getIK_backprop(dL_IKZ, dIKZ{s});
% %     
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% % 
% % 
% % 
% % if mult_NKZ ~= 0
% %   
% %   GSM_NKZ = data.prior.height.NKZ.(params.PRIOR_PDF);
% %   [KZpyr, dKZ] = buildKpyr(Z, params.ZPYR_DEPTH);
% %   
% %   pind_mass = params.ZPYR_ROOT{1}.^([1:length(KZpyr)]'-1);
% %   
% %   d_loss_Zpyr = {};
% %   
% %   for s = 1:length(KZpyr)
% %     
% %     mask = data.valid_pyr{s};
% %   
% %     mult = mult_NKZ * pind_mass(s) / length(KZpyr);
% %     
% %     C = GSM_NKZ{s};
% %     
% %     KZ = KZpyr{s};
% %     
% %     [GKZ, dGKZ] = getGradNorm(KZ);
% %     
% %     KZb = blur3(KZ);
% %     
% %     KZ_pad = sqrt(KZb.^2 + 10^-4);
% %     NKZ = GKZ ./ KZ_pad;
% %     
% %     NKZ = max(epsilon, NKZ);
% %     
% %     [L, dL_mask] = interp1_fixed_sum_fast(NKZ(mask), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
% %     
% %     dL_mask(NKZ(mask) <= epsilon) = 0;
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_NKZ = zeros(size(KZ));
% %     dL_NKZ(mask) = mult * dL_mask;
% %     
% %     dL_GKZ = dL_NKZ ./ KZ_pad;
% %     dL_KZb = -(GKZ .* KZb .* dL_NKZ) ./ (KZ_pad.^3);
% %     dL_KZ = blur3(dL_KZb) + getGradNorm_backprop(dL_GKZ, dGKZ);
% %     
% %     d_loss_Zpyr{s} = getK_backprop(dL_KZ, dKZ{s});
% %     
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% 
% 
% 
% if mult_MKZ ~= 0
%   GSM_MKZ = data.prior.height.MKZ.(params.PRIOR_PDF);
%   [KZpyr, dKZ] = buildKpyr(Z, params.ZPYR_DEPTH);
%   pind_mass = params.ZPYR_ROOT{1}.^([1:length(KZpyr)]'-1);
%   pind_mass = pind_mass / ((params.Z_MEDIAN_HALFWIDTH*2+1)^2);
%   
%   d_loss_Zpyr = {};
%   
%   for s = 1:length(KZpyr)
%     
%     M = data.median_filter_mats{s};
%     
%     if params.ONE_SHAPE_MODEL
%       C = GSM_MKZ{1};
%     else
%       C = GSM_MKZ{s};
%     end
%     
%     mult = mult_MKZ * pind_mass(s) / length(KZpyr);
%     
%     KZ = KZpyr{s};
%     
%     MKZ = M * KZ(:);
%     
%     abs_MKZ = abs(MKZ);
%     sign_MKZ = sign(MKZ);
% 
%     if strcmp(params.PRIOR_PDF, 'GSM')
%       
%       abs_MKZ = max(epsilon, abs_MKZ);
%       [L, dL_mask] = interp1_fixed_sum_fast(abs_MKZ(:), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%       dL_mask(abs_MKZ <= epsilon) = 0;
% 
%     elseif strcmp(params.PRIOR_PDF, 'myGaussian')
%     
%       C = data.prior.height.MKZ.myGaussian{s};
%       [L, dL_mask] = myGaussian_loss(abs_MKZ, C.alpha, C.beta, C.epsilon);
%       
%     end
%     
%     dL_mask = dL_mask .* sign_MKZ;
%     
%     loss_Z = loss_Z + mult * sum(L(:));
%     dL_MKZ = mult * reshape(dL_mask' * M, size(KZ));
%     
%     d_loss_Zpyr{s} = getK_backprop(dL_MKZ, dKZ{s});
%     
%   end
%   
%   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% 
% end
% 
% 
% 
% if mult_MIKZ ~= 0
%   GSM_MIKZ = data.prior.height.MIKZ.(params.PRIOR_PDF);
%   [IKZpyr, dIKZ] = buildIKpyr(Z, params.ZPYR_DEPTH);
%   pind_mass = params.ZPYR_ROOT{1}.^([1:length(IKZpyr)]'-1);
%   pind_mass = pind_mass / ((params.Z_MEDIAN_HALFWIDTH*2+1)^2);
%   
%   d_loss_Zpyr = {};
%   
%   for s = 1:length(IKZpyr)
%     
%     M = data.median_filter_mats{s};
%     
%     if params.ONE_SHAPE_MODEL
%       C = GSM_MIKZ{1};
%     else
%       C = GSM_MIKZ{s};
%     end
%     
%     mult = mult_MIKZ * pind_mass(s) / length(IKZpyr);
%     
%     IKZ = IKZpyr{s};
%     
%     MIKZ = M * IKZ(:);
%     
%     abs_MIKZ = abs(MIKZ);
%     sign_MIKZ = sign(MIKZ);
% 
%     if strcmp(params.PRIOR_PDF, 'GSM')
%       
%       abs_MIKZ = max(epsilon, abs_MIKZ);
%       [L, dL_mask] = interp1_fixed_sum_fast(abs_MIKZ(:), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
%       dL_mask(abs_MIKZ <= epsilon) = 0;
% 
%     elseif strcmp(params.PRIOR_PDF, 'myGaussian')
%     
%       C = data.prior.height.MIKZ.myGaussian{s};
%       [L, dL_mask] = myGaussian_loss(abs_MIKZ, C.alpha, C.beta, C.epsilon);
%       
%     end
%     
%     dL_mask = dL_mask .* sign_MIKZ;
%     
%     loss_Z = loss_Z + mult * sum(L(:));
%     dL_MIKZ = mult * reshape(dL_mask' * M, size(IKZ));
%     
%     d_loss_Zpyr{s} = getIK_backprop(dL_MIKZ, dIKZ{s});
%     
%   end
%   
%   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% 
% end
% 
% 
% % if mult_GKZ ~= 0
% % 
% %   GSM_GKZ = data.prior.height.GKZ.(params.PRIOR_PDF);
% %   [KZpyr, dKZ] = buildKpyr(Z, params.ZPYR_DEPTH);
% %   
% %   pind_mass = params.ZPYR_ROOT{1}.^([1:length(KZpyr)]'-1);
% %   
% %   d_loss_Zpyr = {};
% %   
% %   for s = 1:length(KZpyr)
% %     
% %     mask = data.valid_pyr{s};
% %     
% %     mult = mult_GKZ * pind_mass(s) / length(KZpyr);
% %     
% %     if params.ONE_SHAPE_MODEL
% %       C = GSM_GKZ{1};
% %     else
% %       C = GSM_GKZ{s};
% %     end
% %     
% %     if mult ~= 0
% %       
% %       KZ = KZpyr{s};
% %       
% %       [GKZ, dGKZ] = getGradNorm(KZ);
% %       
% %       GKZ = max(epsilon, GKZ);
% %       
% %       [L, dL_mask] = interp1_fixed_sum_fast(GKZ(mask), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
% %       
% %       dL_mask(GKZ(mask) <= epsilon) = 0;
% %       
% %       loss_Z = loss_Z + mult * sum(L(:));
% %       dL_GKZ = zeros(size(KZ));
% %       dL_GKZ(mask) = mult * dL_mask;
% %       
% %       d_loss_Zpyr{s} = getK_backprop(getGradNorm_backprop(dL_GKZ, dGKZ), dKZ{s});
% %       
% %     else
% %       d_loss_Zpyr{s} = zeros(size(mask));
% %     end
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% %   
% % 
% % 
% % if mult_GIKZ ~= 0
% % 
% %   GSM_GIKZ = data.prior.height.GIKZ.(params.PRIOR_PDF);
% %   [IKZpyr, dIKZ] = buildIKpyr(Z, params.ZPYR_DEPTH);
% %   
% %   pind_mass = params.ZPYR_ROOT{1}.^([1:length(IKZpyr)]'-1);
% %   
% %   d_loss_Zpyr = {};
% %   
% %   for s = 1:length(IKZpyr)
% %     
% %     mask = data.valid_pyr{s};
% %     
% %     mult = mult_GIKZ * pind_mass(s) / length(IKZpyr);
% %     
% %     IKZ = IKZpyr{s};
% %     
% %     [GIKZ, dGIKZ] = getGradNorm(IKZ);
% %     
% %     GIKZ = max(epsilon, GIKZ);
% %     
% %     [L, dL_mask] = interp1_fixed_sum_fast(GIKZ(mask), GSM_GIKZ{s}.lut.bin_range(1), GSM_GIKZ{s}.lut.bin_width, GSM_GIKZ{s}.lut.F_cost);
% %     
% %     dL_mask(GIKZ(mask) <= epsilon) = 0;
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_GIKZ = zeros(size(IKZ));
% %     dL_GIKZ(mask) = mult * dL_mask;
% %     
% %     d_loss_Zpyr{s} = getIK_backprop(getGradNorm_backprop(dL_GIKZ, dGIKZ), dIKZ{s});
% %     
% %   end
% %   
% %   
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% 
% 
% 
% 
% 
% loss_Z = loss_Z / numel(Z);
% d_loss_Z = d_loss_Z / numel(Z);
% 
% 
% % mult = mult_height_contour;
% % 
% % if nargin < 4
% %   [N, dN_Z] = getNormals_conv(Z);
% % end
% % N1 = N(:,:,1);
% % N2 = N(:,:,2);
% % 
% % 
% % mask = data.border.mask;
% % normal = data.border.normal;
% % 
% % d1 = N1(mask) - normal(:,1);
% % d2 = N2(mask) - normal(:,2);
% % 
% % d_loss_N1 = zeros(size(N1));
% % d_loss_N2 = zeros(size(N2));
% % 
% % dmag = sqrt(d1.^2 + d2.^2);
% % loss_Z = loss_Z + mult * sum(dmag);
% % d_loss_N1(mask) = mult * (d1 ./ dmag);
% % d_loss_N2(mask) = mult * (d2 ./ dmag);
% % 
% % 
% % 
% % mult = params.MULT_OPTS.srfs.height.slant{1};
% % 
% % if mult ~= 0
% %     
% %   N3 = N(:,:,3);
% %   mask = data.valid;
% %   
% %   loss_Z = loss_Z + mult * sum(-log(2*N3(mask)));  % Actual PDF
% % %   loss_Z = loss_Z + mult * sum(-log(N3(mask))); % Don't need the factor of 2.
% %   d_loss_N3 = mask .* (-mult)./max(eps,N3);
% %   
% % else
% %   d_loss_N3 = 0;
% % end
% % 
% % 
% % d_loss_D1 = d_loss_N1 .* dN_Z.F1_1 + d_loss_N2 .* dN_Z.F2_1 + d_loss_N3 .* dN_Z.F1_3;
% % d_loss_D2 = d_loss_N1 .* dN_Z.F1_2 + d_loss_N2 .* dN_Z.F2_2 + d_loss_N3 .* dN_Z.F2_3;
% % 
% % d_loss_Z = d_loss_Z + conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);
% % 
% % 
% % 
% % loss_Z = loss_Z / numel(Z);
% % d_loss_Z = d_loss_Z / numel(Z);
% 
% 
% 
% 
% % function [loss_Z, d_loss_Z] = priorZ(Z, data, params, N, dN_Z)
% % 
% % MULT = params.MULT_OPTS.srfs;
% % 
% % if strcmp(params.CURVATURE_MODE, 'K')
% % 
% %   GSM_GKZ = data.prior.height.GKZ.(params.PRIOR_PDF);
% %   [KZpyr, dKZ, Zpyr] = buildKpyr(Z, params.ZPYR_DEPTH);
% %   
% % elseif strcmp(params.CURVATURE_MODE, 'IK')
% % 
% %   GSM_GKZ = data.prior.height.GIKZ.(params.PRIOR_PDF);
% %   [KZpyr, dKZ, Zpyr] = buildIKpyr(Z, params.ZPYR_DEPTH);
% %   
% % end
% % 
% % pind_mass = params.ZPYR_ROOT{1}.^([1:length(KZpyr)]'-1);
% % 
% % mult_height_contour = MULT.height.contour{1}*100;
% % mult_GKZ = MULT.height.GKZ{1};
% % 
% % loss_Z = 0;
% % d_loss_Z = zeros(size(Z));
% % d_loss_Zpyr = {};
% % 
% %   
% % for s = 1:length(KZpyr)
% %   
% %   mask = data.valid_pyr{s};
% %   
% %   mult = mult_GKZ * pind_mass(s);
% %   
% %   if mult ~= 0
% % 
% %     KZ = KZpyr{s};
% %     
% % %     [GKZ, dGKZ] = getGrad2(KZ);
% % %     mask = mask(1:end-1,1:end-1);
% % %     
% % %     epsilon = 5*10^-3;
% % %     GKZ = max(epsilon, GKZ);
% % %     
% % %     if strcmp(params.PRIOR_PDF, 'myGaussian')
% % %       [LL, dLL_mask] = myGaussian_pdf(GKZ(mask), GSM_GKZ{s}.alpha, GSM_GKZ{s}.beta, GSM_GKZ{s}.epsilon);
% % %       L = -LL;
% % %       dL_mask = -dLL_mask;
% % %     else
% % %       [L, dL_mask] = interp1_fixed_sum_fast(GKZ(mask), GSM_GKZ{s}.lut.bin_range(1), GSM_GKZ{s}.lut.bin_width, GSM_GKZ{s}.lut.F_cost);
% % %     end
% % %     
% % %     dL_mask(GKZ(mask) <= epsilon) = 0;
% % %     
% % %     loss_Z = loss_Z + mult * sum(L(:));
% % %     dL_GKZ = zeros(size(GKZ));
% % %     dL_GKZ(mask) = mult * dL_mask;
% % %     
% % %     if strcmp(params.CURVATURE_MODE, 'K')
% % %       d_loss_Zpyr{s} = getK_backprop(getGrad2_backprop(dL_GKZ, dGKZ), dKZ{s});
% % %     elseif strcmp(params.CURVATURE_MODE, 'IK')
% % %       d_loss_Zpyr{s} = getIK_backprop(getGrad2_backprop(dL_GKZ, dGKZ), dKZ{s});
% % %     end
% %     
% %     
% %     [GKZ, dGKZ] = getGradNorm(KZ);
% %     
% %     epsilon = 5*10^-3;
% %     GKZ = max(epsilon, GKZ);
% % 
% %     if strcmp(params.PRIOR_PDF, 'myGaussian')
% %       [LL, dLL_mask] = myGaussian_pdf(GKZ(mask), GSM_GKZ{s}.alpha, GSM_GKZ{s}.beta, GSM_GKZ{s}.epsilon);
% %       L = -LL;
% %       dL_mask = -dLL_mask;
% %     else
% %       [L, dL_mask] = interp1_fixed_sum_fast(GKZ(mask), GSM_GKZ{s}.lut.bin_range(1), GSM_GKZ{s}.lut.bin_width, GSM_GKZ{s}.lut.F_cost);
% %     end
% %     
% %     dL_mask(GKZ(mask) <= epsilon) = 0;
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_GKZ = zeros(size(KZ));
% %     dL_GKZ(mask) = mult * dL_mask;
% %     
% %     if strcmp(params.CURVATURE_MODE, 'K')
% %       d_loss_Zpyr{s} = getK_backprop(getGradNorm_backprop(dL_GKZ, dGKZ), dKZ{s});
% %     elseif strcmp(params.CURVATURE_MODE, 'IK')
% %       d_loss_Zpyr{s} = getIK_backprop(getGradNorm_backprop(dL_GKZ, dGKZ), dKZ{s});
% %     end
% %     
% % 
% %   else
% %     d_loss_Zpyr{s} = zeros(size(mask));
% %   end
% % end
% % 
% % 
% % mult = MULT.height.lowflat{1};
% % 
% % if mult ~= 0
% %   s = params.ZPYR_DEPTH;
% %   
% %   delta = Zpyr{s};
% %   mask = data.valid_pyr{s};
% %   
% %   loss_Z = loss_Z + mult*0.5*sum(delta(mask).^2);
% %   d_loss_Zpyr{s} = d_loss_Zpyr{s} + mult*(mask.*delta);
% %   
% % end
% % 
% % 
% % d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % 
% % 
% % mult = mult_height_contour;
% % 
% % if nargin < 4
% %   [N, dN_Z] = getNormals_conv(Z);
% % end
% % N1 = N(:,:,1);
% % N2 = N(:,:,2);
% % 
% % 
% % mask = data.border.mask;
% % normal = data.border.normal;
% % 
% % d1 = N1(mask) - normal(:,1);
% % d2 = N2(mask) - normal(:,2);
% % 
% % d_loss_N1 = zeros(size(N1));
% % d_loss_N2 = zeros(size(N2));
% % 
% % dmag = sqrt(d1.^2 + d2.^2);
% % loss_Z = loss_Z + mult * sum(dmag);
% % d_loss_N1(mask) = mult * (d1 ./ dmag);
% % d_loss_N2(mask) = mult * (d2 ./ dmag);
% % 
% % 
% % 
% % mult = params.MULT_OPTS.srfs.height.slant{1};
% % 
% % if mult ~= 0
% %     
% %   N3 = N(:,:,3);
% %   mask = data.valid;
% %   
% % %   mask = true(size(mask));
% % %   mask([1,end],:) = false;
% % %   mask(:,[1,end]) = false;
% %   
% % %   loss_Z = loss_Z + mult * sum(-log(2*N3(mask)));  % Actual PDF
% %   loss_Z = loss_Z + mult * sum(-log(N3(mask))); % Don't need the factor of 2.
% %   d_loss_N3 = mask .* (-mult)./max(eps,N3);
% %   
% % else
% %   d_loss_N3 = 0;
% % end
% % 
% % 
% % d_loss_D1 = d_loss_N1 .* dN_Z.F1_1 + d_loss_N2 .* dN_Z.F2_1 + d_loss_N3 .* dN_Z.F1_3;
% % d_loss_D2 = d_loss_N1 .* dN_Z.F1_2 + d_loss_N2 .* dN_Z.F2_2 + d_loss_N3 .* dN_Z.F2_3;
% % 
% % d_loss_Z = d_loss_Z + conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);
% % 
% % 
% % 
% % 
% % mult = params.MULT_OPTS.srfs.height.multislant{1};
% % 
% % if mult ~= 0
% %   
% %   d_loss_Zpyr = {};
% %   for s = 1:length(Zpyr)
% %   
% %     mask = data.valid_pyr{s};
% %     
% % %     mask = true(size(mask));
% % %     mask([1,end],:) = false;
% % %     mask(:,[1,end]) = false;
% %     
% %     mult_band = mult * pind_mass(s);
% %     Zband = Zpyr{s};
% % 
% %     if (s >= 2)
% %       [N, dN_Z] = getNormals_conv(Zband);
% %     end
% %     N3 = N(:,:,3);
% % 
% %     loss_Z = loss_Z + mult_band * sum(-log(2*N3(mask)));
% %     dL = mask .* (-mult_band)./max(eps,N3);
% %     
% %     d_loss_Zpyr{s} = conv3d( dL .* dN_Z.F1_3, dN_Z.f1) + conv3d( dL .* dN_Z.F2_3, dN_Z.f2);
% %     
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% %   
% % end
% % 
% % 
% % 
% % 
% % mult = params.MULT_OPTS.srfs.height.frontness{1};
% % 
% % if mult ~= 0
% % 
% %   [GZ, dGZ] = getGradNorm(Z);
% %     
% %   mask = data.valid;
% %     
% %   GSM = data.prior.height.GZ.(params.PRIOR_PDF){1};
% %   
% %   [L, dL_mask] = interp1_fixed_sum_fast(GZ(mask), GSM.lut.bin_range(1), GSM.lut.bin_width, GSM.lut.F_cost);
% %   dL = zeros(size(Z));
% %   dL(mask) = dL_mask;
% %   
% %   loss_Z = loss_Z + mult*L;
% %   d_loss_Z = d_loss_Z + getGradNorm_backprop(mult*dL, dGZ);
% %   
% % end
% % 
% % 
% % mult_GZ = params.MULT_OPTS.srfs.height.multifront{1};
% % 
% % if mult_GZ ~= 0
% %   GSM_GZ = data.prior.height.GZ.(params.PRIOR_PDF);
% %   
% %   d_loss_Zpyr = {};
% % 
% %   for s = 1:length(Zpyr)
% %     
% %     mask = data.valid_pyr{s};
% %     
% %     mult = mult_GZ * pind_mass(s);
% %     Zband = Zpyr{s};
% %     
% %     [GZ, dGZ] = getGradNorm(Zband);
% %     
% %     epsilon = 5*10^-3;
% %     GZ = max(epsilon, GZ);
% %     
% %     [L, dL_mask] = interp1_fixed_sum_fast(GZ(mask), GSM_GZ{s}.lut.bin_range(1), GSM_GZ{s}.lut.bin_width, GSM_GZ{s}.lut.F_cost);
% %     
% %     dL_mask(GZ(mask) <= epsilon) = 0;
% %     
% %     loss_Z = loss_Z + mult * sum(L(:));
% %     dL_GZ = zeros(size(GZ));
% %     dL_GZ(mask) = mult * dL_mask;
% %     
% %     d_loss_Zpyr{s} = getGradNorm_backprop(dL_GZ, dGZ);
% %   
% %   end
% %   
% %   d_loss_Z = d_loss_Z + reconLpyr_simple(d_loss_Zpyr, [1;4;6;4;1]/(16*sqrt(2)), 'zero');
% % 
% % end
% % 
% % 
% % mult = params.MULT_OPTS.srfs.height.flat{1}/20000;
% % if mult ~= 0
% %   
% %   delta = Z;
% %   
% %   mask = data.valid;
% %   loss_Z = loss_Z + mult*0.5*sum(delta(mask).^2);
% %   d_loss_Z = d_loss_Z + mult*(mask.*delta);
% %   
% % end
% % 
% % loss_Z = loss_Z / numel(Z);
% % d_loss_Z = d_loss_Z / numel(Z);
% % 
