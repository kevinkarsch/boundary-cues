function [loss_Z, d_loss_Z, losses_Z] = myPriorZ(Z, data, params, N, dN_Z)

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
  [L, dL_mask] = interp1_fixed_sum_fast(soft_MKZ(:), C.lut.bin_range(1), C.lut.bin_width, C.lut.F_cost);
  if epsilon > 0
    dL_mask = dL_mask .* (abs_MKZ ./ max(eps,soft_MKZ));
  end

  dL_mask = dL_mask .* sign_MKZ;
  
  dL_KZ = reshape(M_T * dL_mask, size(KZ));%reshape(dL_mask' * M, size(KZ));
  
  losses_Z.smooth = mult_mod_smooth * sum(L(:));
  d_loss_Z = d_loss_Z + mult_smooth * mult_mod_smooth * getK_backprop_fast_mat(dL_KZ, dKZ);
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


%%%%%%%%%%%%%% New terms (KK 10/16/2012)
%Fold term
fold_mult = params.MULT_OPTS.saifs.height.fold{1};
losses_Z.fold = 0;
d_fold_Z = zeros(size(Z));
for i=1:numel(data.fold)
    [fold_loss_i, d_fold_Z_i] = foldLoss(N, dN_Z, data.fold{i});
    losses_Z.fold = losses_Z.fold + fold_loss_i;
    d_fold_Z = d_fold_Z + d_fold_Z_i;
end
loss_Z = loss_Z + fold_mult * losses_Z.fold;
d_loss_Z = d_loss_Z + fold_mult .* d_fold_Z;

%Hinge term (for contact points)
hinge_mult = params.MULT_OPTS.saifs.height.hinge{1};
losses_Z.hinge = 0;
d_hinge_Z = zeros(size(Z));
for hingeDist_i = 1:length(data.hingeDist)
  [sub_loss, d_loss_Z1, d_loss_Z2] = hingeDist(Z(data.hingeDist{hingeDist_i}.mask1), Z(data.hingeDist{hingeDist_i}.mask2), data.hingeDist{hingeDist_i}.delta); 
  losses_Z.hinge = losses_Z.hinge + sub_loss;
  d_hinge_Z(data.hingeDist{hingeDist_i}.mask1) = d_hinge_Z(data.hingeDist{hingeDist_i}.mask1) + d_loss_Z1;
  d_hinge_Z(data.hingeDist{hingeDist_i}.mask2) = d_hinge_Z(data.hingeDist{hingeDist_i}.mask2) + d_loss_Z2;
end
loss_Z = loss_Z + hinge_mult * losses_Z.hinge;
d_loss_Z = d_loss_Z + hinge_mult .* d_hinge_Z;
end


%%%%%%%%%%%%%% Fold loss from Jon (KK 10/16/2012)
function [loss, d_loss_Z] = foldLoss(N, dN_Z, fold)

MULT_FOLD = 1;

FOLD_THRESH = 1/sqrt(2); % Seems to encourage folds to be ~90 degrees
FOLD_OFFSET = 1; % How many pixels from the fold we measure surface normals

N1 = N(:,:,1);
N2 = N(:,:,2);
N3 = N(:,:,3);

loss = 0;

P_left  = round(fold.position + FOLD_OFFSET*fold.normal);
P_right = round(fold.position - FOLD_OFFSET*fold.normal);
T_fold = fold.tangent;

P_left(P_left<1) = 1; P_left(P_left(:,1)>size(N,1),1) = size(N,1); P_left(P_left(:,2)>size(N,2),2) = size(N,2);
P_right(P_right<1) = 1; P_right(P_right(:,1)>size(N,1),1) = size(N,1); P_right(P_right(:,2)>size(N,2),2) = size(N,2);

idx_left = sub2ind(size(N1), P_left(:,1), P_left(:,2));
idx_right = sub2ind(size(N1), P_right(:,1), P_right(:,2));

N_left = [N1(idx_left), N2(idx_left), N3(idx_left)];
N_right = [N1(idx_right), N2(idx_right), N3(idx_right)];

% Cross product of N_left and N_right, dotted with T_fold;
P = T_fold(:,2) .* N_left(:,3) .* N_right(:,1) ...
  - T_fold(:,1) .* N_left(:,3) .* N_right(:,2) ...
  - T_fold(:,2) .* N_left(:,1) .* N_right(:,3) ...
  + T_fold(:,1) .* N_left(:,2) .* N_right(:,3);

% Folds should be *at least* some angle apart
delta = max(0, FOLD_THRESH - P);
loss = loss + MULT_FOLD * sum(delta);
d_loss_P = MULT_FOLD * -(delta > 0);

% Backpropagate through the triple product and onto the normals
d_loss_N_left  = bsxfun(@times, d_loss_P, ...
  [-T_fold(:,2).*N_right(:,3),  ...
    T_fold(:,1).*N_right(:,3),  ...
    T_fold(:,2).*N_right(:,1) - T_fold(:,1).*N_right(:,2)]);
  
d_loss_N_right = bsxfun(@times, d_loss_P, ...
  [ T_fold(:,2).*N_left(:,3), ...
   -T_fold(:,1).*N_left(:,3),  ...
   -T_fold(:,2).*N_left(:,1)  + T_fold(:,1).*N_left(:,2)]);

d_loss_N = {};
for c = 1:3
  d_loss_N{c} = zeros(size(N1));
  d_loss_N{c} = d_loss_N{c} + sparse(P_left(:,1), P_left(:,2), d_loss_N_left(:,c), size(N1,1), size(N1,2));
  d_loss_N{c} = d_loss_N{c} + sparse(P_right(:,1), P_right(:,2), d_loss_N_right(:,c), size(N1,1), size(N1,2));
end
d_loss_N = cat(3, d_loss_N{:});

% Backpropagate onto Z
d_loss_D1 = dN_Z.F1_1 .* d_loss_N(:,:,1) + dN_Z.F1_2 .* d_loss_N(:,:,2) + dN_Z.F1_3 .* d_loss_N(:,:,3);
d_loss_D2 = dN_Z.F2_1 .* d_loss_N(:,:,1) + dN_Z.F2_2 .* d_loss_N(:,:,2) + dN_Z.F2_3 .* d_loss_N(:,:,3);
d_loss_Z = conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);

end
