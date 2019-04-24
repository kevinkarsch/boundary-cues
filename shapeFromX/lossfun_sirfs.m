function [loss, d_loss, state_return] =  lossfun_saifs(state, data, params)

Z = reconLpyr_simple(state.Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');

L_white = state.L_white / data.L_precond;

if params.WHITEN_LIGHT
  L = bsxfun(@plus, L_white * data.L_whiten_params.inverse, data.L_whiten_params.mean);
else
  L = L_white;
end

L = reshape(L, [9,3]);

[N, dN_Z] = getNormals_conv(Z);
N_vec = reshape(N, [], 3);

[loss_Z, d_loss_Z, losses_Z] = priorZ(Z, data, params, N, dN_Z);


if params.FORCE_MONOCHROME_LIGHT
  L0 = L;
  L = repmat(mean(L,2), 1, 3);
end

S = {};
dS_Z = {};
dS_L = {};
for c = 1:3
    
  [s, ds_n, ds_l] = renderSH_helper(N_vec, L(:,c));
    
  %   s = s + state.Lshift;
  S{c} = reshape(s, size(Z));
  dS_N = reshape(ds_n, [size(Z), 3]);
  
  dS_Z{c}.F1 = dN_Z.F1_1 .* dS_N(:,:,1) + dN_Z.F1_2 .* dS_N(:,:,2) + dN_Z.F1_3 .* dS_N(:,:,3);
  dS_Z{c}.F2 = dN_Z.F2_1 .* dS_N(:,:,1) + dN_Z.F2_2 .* dS_N(:,:,2) + dN_Z.F2_3 .* dS_N(:,:,3);
  
  dS_L{c} = ds_l;
end
S = cat(3, S{:});

A = data.true.log_im - S;
  
if strcmp(params.ALBEDO_MODEL, 'ours')
  [loss_A, d_loss_A, losses_A] = priorA(A, data, params);
elseif strcmp(params.ALBEDO_MODEL, 'rgb')
  [loss_A, d_loss_A, losses_A] = priorA_rgb(A, data, params);
elseif strcmp(params.ALBEDO_MODEL, 'yuv')
  [loss_A, d_loss_A, losses_A] = priorA_yuv(A, data, params);
end

d_loss_S = -d_loss_A;

d_loss_D1 = 0;
d_loss_D2 = 0;
for c = 1:3
  d_loss_D1 = d_loss_D1 + d_loss_S(:,:,c) .* dS_Z{c}.F1;
  d_loss_D2 = d_loss_D2 + d_loss_S(:,:,c) .* dS_Z{c}.F2;
end  
  
d_loss_ZA = conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);


if params.FORCE_MONOCHROME_LIGHT
  [loss_L, d_loss_L] = priorL(L0, data, params);
  d_loss_L = d_loss_L + repmat(mean(dS_L{c}' * reshape(d_loss_S, [], 3),2), [1, 3]);
else
  [loss_L, d_loss_L] = priorL(L, data, params);
  d_loss_L = d_loss_L + dS_L{c}' * reshape(d_loss_S, [], 3);
end

loss = loss_A + loss_Z + loss_L;
d_loss_Z = d_loss_ZA + d_loss_Z;

d_loss.Zpyr = buildGpyr_simple(d_loss_Z, length(state.Zpyr), params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');

if params.DISABLE_ZPYR > 0
  d_loss.Zpyr(1:params.DISABLE_ZPYR) = cellfun(@(x) zeros(size(x)), d_loss.Zpyr(1:params.DISABLE_ZPYR), 'UniformOutput', false);
end

if params.WHITEN_LIGHT
  d_loss_L_white = (d_loss_L(:)' * data.L_whiten_params.inverse);
else
  d_loss_L_white = d_loss_L(:);
end

d_loss.L_white = d_loss_L_white / data.L_precond;
% d_loss.Lvec = d_loss.Lvec / size(LCs,3);

if ~params.SOLVE_LIGHT
  d_loss.L_white = zeros(size(d_loss.L_white));
end


state_return.height = Z;
state_return.light = L;
state_return.normal = N;
state_return.shading = exp(S);
state_return.albedo = exp(A);


global global_loss_best
if isempty(global_loss_best)
  global_loss_best = inf;
end

global global_loss_best_video
if isempty(global_loss_best_video)
  global_loss_best_video = inf;
end

if (loss <= global_loss_best)
  global_loss_best = loss;
  
  if params.MAKE_VIDEO && (loss <= (global_loss_best_video - params.VIDEO_LOSS_EPS))
    global_loss_best_video = loss;
    
    S = {};
    for c = 1:3
      S{c} = reshape(renderSH_helper(N_vec, L(:,c)), size(Z));
    end
    S = cat(3, S{:});
    
    A = data.true.log_im - S;
    
    invalid = repmat(~data.valid, [1,1,3]);
    Z(~data.valid) = nan;
    S(invalid) = nan;
    A(invalid) = nan;
    
    I = data.true.im;
    
    Lv = visSH_color(L, [128, 192]);
    
    A = min(1, exp(A));
    S = exp(S);
    Lm = (1.25*max(Lv(:)));
    S = S ./ Lm;
    
    A(isnan(A)) = 1;
    S(isnan(S)) = 1;
    
    Zv = visualizeDEM(Z);
    
    Nv = visualizeNormals(Z);
    
    Lv = visSH_color(L, round([size(Nv,1), size(Nv,2)]/4));
    Lv(isnan(Lv)) = 1;
    Lv = min(1, Lv);
    
    Lv = [ones(floor((2*size(Nv,1) - size(Lv,1))/2), size(Lv,2), 3); Lv; ones(ceil((2*size(Nv,1) - size(Lv,1))/2), size(Lv,2), 3)];
    

%     Lv = visSH_color(L, round([size(Nv,1), size(Nv,2)]/2));
%     Lv(isnan(Lv)) = 1;
%     Lv = min(1, Lv);
%     
%     Lv = [ones(floor((size(Nv,1) - size(Lv,1))/2), size(Lv,2), 3); Lv; ones(ceil((size(Nv,1) - size(Lv,1))/2), size(Lv,2), 3)];
%     Lv = [ones(size(Lv,1), floor((size(Nv,2) - size(Lv,2))/2), 3), Lv, ones(size(Lv,1), ceil((size(Nv,2) - size(Lv,2))/2), 3)];
    
    mat1 = mod([1:size(I,1)], 16) <= 8;
    mat2 = mod([1:size(I,2)], 16) <= 8;
    checker = bsxfun(@xor, mat2, mat1');
    I_bg = repmat(1-checker*0.2, [1,1,3]);
    I(invalid) = I_bg(invalid);
    
    I(1:4,:,1) = 0.75;
    I(:,1:4,1) = 0.75;
    I(end-3:end,:,1) = 0.75;
    I(:,end-3:end,1) = 0.75;
    
    I(1:4,:,2) = 0.25;
    I(:,1:4,2) = 0.25;
    I(end-3:end,:,2) = 0.25;
    I(:,end-3:end,2) = 0.25;
    
    I(1:4,:,3) = 0.25;
    I(:,1:4,3) = 0.25;
    I(:,end-3:end,3) = 0.25;
    I(end-3:end,:,3) = 0.25;
    
%     I(invalid) = 1;
%     I1 = I(:,:,1);
%     I2 = I(:,:,2);
%     I3 = I(:,:,3);
%     I1(~data.valid) = 0.75;
%     I2(~data.valid) = 0.25;
%     I3(~data.valid) = 0.25;
%     I = cat(3, I1, I2, I3);
    
%     V1 = [I, Zv, Nv];
%     V2 = [A, S, Lv];
%     V = [V1; V2];
    
    V1 = [I, Zv];
    V2 = [A, S];
    V = [[V1; V2], Lv];
    
    global GLOBAL_VIDEO_COUNT;
    imwrite(uint8(round(255*V)), [data.video_folder, 'frame_', num2str(GLOBAL_VIDEO_COUNT, '%05d'), '.jpg'], 'Quality', 100)
%     if GLOBAL_VIDEO_COUNT == 1
%       
%       imwrite(I, [data.video_folder, 'image.png'])
%     end
    
%     fprintf(data.video_log, '%e ', [loss_A, loss_Z, loss_L]);
%     fprintf(data.video_log, '%e ', [loss]);
%     fprintf(data.video_log, '\n');

    GLOBAL_VIDEO_COUNT = GLOBAL_VIDEO_COUNT + 1;
    
  else
  
    if params.DO_DISPLAY
      
      global last_display
      if isempty(last_display) || (etime(clock, last_display) > params.DISPLAY_PERIOD)
        last_display = clock;
        
        S = {};
        for c = 1:3          
          S{c} = reshape(renderSH_helper(N_vec, L(:,c)), size(Z));
        end
        S = cat(3, S{:});
        
        A = data.true.log_im - S;
        
        invalid = repmat(~data.valid, [1,1,3]);
        Z(~data.valid) = nan;
        S(invalid) = nan;
        A(invalid) = nan;
        
        I = data.true.im;
        
        A = exp(A);
        S = exp(S);
        
        %       A = max(0, min(1, A));
        %       S = max(0, min(1, S));
        
        state = struct('normal', N, 'height', Z, 'albedo', A, 'shading', S, 'light', L);
        clear N Z A S L
                        
        if isfield(data.true, 'height')
          
          [err, junk] = getError(state, data.true);
          
          s = size(state.height,1);
          
          Lv = visSH_color(state.light, [s, 150]);
          Lv(isnan(Lv)) = 0;
          
          if size(Lv,1) > s
            Lv = max(0, min(1, imresize(Lv, [s, 150])));
          else
            pad1 = floor((s - size(Lv,1))/2);
            pad2 = ceil((s - size(Lv,1))/2);
            Lv = padarray(padarray(Lv, [pad2], 0, 'post'), [pad1], 0, 'pre');
          end
          
%           m = max(Lv(:));
%           Lv = Lv ./ m;
%           state.shading = state.shading ./ m;
%           state.albedo = state.albedo .* m;
          Lv = max(0, min(1, Lv));
          
          if size(data.true.light,2) == 1
            data.true.light = repmat(data.true.light, [1,3]);
          end
          Lv_true = visSH_color(data.true.light, [size(state.height,1), 150]);
          Lv_true(isnan(Lv_true)) = 0;
          Lv_true = max(0, min(1, Lv_true));
          
          I(invalid) = 0;
          data.true.shading(invalid) = 0;
          Lv(isnan(Lv)) = 0;
          state.albedo(invalid) = 0;
          state.shading(invalid) = 0;
          
%           mat1 = mod([1:size(I,1)], 16) <= 8;
%           mat2 = mod([1:size(I,2)], 16) <= 8;
%           checker = bsxfun(@xor, mat2, mat1');
%           I_bg = repmat(1-checker*0.2, [1,1,3]);
%           I(invalid) = I_bg(invalid);
          
          Z = state.height;
          shift = mean(Z(~isnan(data.true.height)) - data.true.height(~isnan(data.true.height)));
          Z = Z - shift;
          Zv = visualizeDEM([Z; data.true.height]);
%           Zv = [visualizeDEM(state.height); visualizeDEM(data.true.height)];
%           Zv = visualizeDEM([state.height; data.true.height]);
          Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
          Nv = visualizeNormals_color(state.normal);
          Ntv = visualizeNormals_color(data.true.normal);
          Ntv(isnan(data.true.normal)) = 0;
          Nv(invalid) = 0;
          I_pad = padarray(padarray(I, floor(size(I,1)/2), 0, 'pre'), ceil(size(I,1)/2), 0, 'post'); % hah, ipad.
          V = min(1, [I_pad, Zv, [Nv; Ntv], [state.albedo; data.true.albedo], [state.shading; data.true.shading], [Lv; Lv_true]]);
          
        else
          
          Lv = visSH_color(state.light, [size(state.height,1), 150]);
          Lv(isnan(Lv)) = 0;
          Lv = max(0, min(1, Lv));
          
          I(invalid) = 0;
          Lv(isnan(Lv)) = 0;
          state.albedo(invalid) = 0;
          state.shading(invalid) = 0;
          
          state.albedo = state.albedo ./ max(state.albedo(:));
          state.shading = state.shading ./ max(state.shading(:));
          Lv = Lv ./ max(Lv(:));
          
          Zv = visualizeDEM([state.height]);
          Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
          Nv = visualizeNormals_color(state.normal);
          Nv(invalid) = 0;
          V = min(1, [I, Zv, Nv, state.albedo, state.shading, Lv]);
          
        end
        
        figure(1);
        imagesc(V);
        axis image off;
        set(gca, 'PlotBoxAspectRatio', [1 1 1])
        set(gca, 'Position', [0 0 1 1]);
        set(gcf, 'Position', [1, 500, 1600, 1600/size(V,2)*size(V,1)])
        
        if isfield(data.true, 'height')
          %     figure(3)
          %     a = flattenPyr(Apyr);
          %     imagesc(min(1,exp(a))); imtight;
          %     set(gcf, 'Position', [1, 1, 1600, 1600/size(V,2)*size(V,1)])
          %     drawnow;
          
          global global_errors;
          global_errors{end+1} = err;
          light_errs = cellfun(@(x) x.light, global_errors);
%           curvature_errs = cellfun(@(x) x.curvature, global_errors);
          normal_errs = cellfun(@(x) x.normal, global_errors);
          %     survey_errs = cellfun(@(x) x.survey, global_errors);
          albedo_errs = cellfun(@(x) x.albedo, global_errors);
          shading_errs = cellfun(@(x) x.shading, global_errors);
          grosse_errs = cellfun(@(x) x.grosse, global_errors);
          avg_errs = cellfun(@(x) x.avg, global_errors);
          
          
          figure(2);
          plot(light_errs ./ max(eps,light_errs(1)), 'Color', [0.5, 0.5, 0.5]); hold on;
%           plot(curvature_errs ./ max(eps,curvature_errs(1)), 'Color', [0.5, 0.5, 0.5]); hold on;
          plot(normal_errs ./ max(eps,normal_errs (1)), 'k-');
          
          plot(grosse_errs ./ max(eps,grosse_errs(1)), 'b-');
          
          plot(albedo_errs ./ max(eps,albedo_errs(1)), 'r-');
          plot(shading_errs ./ max(eps,shading_errs(1)), 'g-');
          
          plot(avg_errs ./ max(eps,avg_errs(1)), 'm-');
          
          axis square
          set(gca, 'YLim', [.009, 2])
          set(gca, 'YScale', 'log');
          title(['gray = LightMSE, black = NMSE, b = LMSE, r = AMSE, g = SMSE, m = avg'])
%           title(['gray = KMSE, black = NMSE, b = LMSE, r = AMSE, g = SMSE, m = avg'])
          grid on;
          hold off;
        end
        drawnow;
        
      end
    end
  end
end










% function [loss, d_loss, Z_return, L_return, Ls] =  lossfun_saifs(state, data, params)
% 
% Z = reconLpyr_simple(state.Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% % Z = reconOLpyr(state.Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% Z_return = Z;
% 
% % figure(4); imagesc_root(getKup(Z), 0.25); colormap_usa; imtight; drawnow;
% 
% Ls_white = state.Ls_white / data.L_precond;
% Ls = bsxfun(@plus, Ls_white * data.L_whiten_params.inverse, data.L_whiten_params.mean);
% Ls = reshape(Ls', 9, 3, []);
% 
% % Ls_vec = reshape(Ls, 27, [])';
% 
% % [Ls_vec_unique, global_Ls_weights, dShrink] = shrinkTransform(Ls_vec, params.LIGHT_SHRINK_EPSILON{1}, params.LIGHT_SHRINK_RATIO{1});
% % Ls_unique = reshape(Ls_vec_unique', [9,3,size(Ls_vec_unique,1)]);
% 
% global global_Ls_weights global_Ls_dependency
% if isempty(global_Ls_weights)
%   global_Ls_weights = ones(size(Ls,3), 1);
% end
% if isempty(global_Ls_dependency)
%   global_Ls_dependency = [1:size(Ls,3)]';
% end
%   
% % Ls = Ls;
% 
% % Ls = Ls;
% % global_Ls_weights = ones(size(Ls,3),1);
% % dShrink = eye(size(Ls,3));
% 
% % epsilon = 0.005;
% % Dsq = distMat(Ls_vec, Ls_vec);
% % C = (Dsq < epsilon)^size(Dsq,2) > 0;
% % 
% % [S, CC] = graphconncomp(sparse(C));
% % [junk, heads] = unique(CC, 'first');
% % Ls_parents = heads(CC);
% % 
% % Ls_avg = nan(size(Ls));
% % for hi = 1:length(heads)
% %   Ls_avg(:,:,heads(hi)) = reshape(mean(Ls_vec(CC == hi,:),1), [9,3]);
% % end
% 
% [N, dN_Z] = getNormals_conv(Z);
% N_vec = reshape(N, [], 3);
% 
% [loss_Z, d_loss_Z, losses_Z] = priorZ(Z, data, params, N, dN_Z);
% 
% losses_A = 10^10*ones(size(Ls,3), 1);
% d_loss_D1s = {};
% d_loss_D2s = {};
% d_loss_Ls = {};
% 
% for L_i = 1:size(Ls,3)%find(global_Ls_weights > 0)'
%   
%   if global_Ls_weights(L_i) == 0
%     
%     d_loss_Ls{L_i} = zeros(9,3);
%     d_loss_D1s{L_i} = zeros(size(Z));
%     d_loss_D2s{L_i} = zeros(size(Z));
%     
%   else
%     
%     L = Ls(:,:,L_i);
%     
%     S = {};
%     dS_Z = {};
%     dS_L = {};
%     for c = 1:3
%       
%       [s, ds_n, ds_l] = renderSH_helper(N_vec, L(:,c));
%       
%       %   s = s + state.Lshift;
%       S{c} = reshape(s, size(Z));
%       dS_N = reshape(ds_n, [size(Z), 3]);
%       
%       dS_Z{c}.F1 = dN_Z.F1_1 .* dS_N(:,:,1) + dN_Z.F1_2 .* dS_N(:,:,2) + dN_Z.F1_3 .* dS_N(:,:,3);
%       dS_Z{c}.F2 = dN_Z.F2_1 .* dS_N(:,:,1) + dN_Z.F2_2 .* dS_N(:,:,2) + dN_Z.F2_3 .* dS_N(:,:,3);
%       
%       dS_L{c} = ds_l;
%     end
%     S = cat(3, S{:});
%     
%     A = data.true.log_im - S;
%     
%     if strcmp(params.ALBEDO_MODEL, 'ours')
%       [loss_A_sub, d_loss_A, losses_A_sub(L_i)] = priorA(A, data, params);
%     elseif strcmp(params.ALBEDO_MODEL, 'rgb')
%       [loss_A_sub, d_loss_A, losses_A_sub(L_i)] = priorA_rgb(A, data, params);
%     elseif strcmp(params.ALBEDO_MODEL, 'yuv')
%       [loss_A_sub, d_loss_A, losses_A_sub(L_i)] = priorA_yuv(A, data, params);
%     end
%     losses_A(L_i) = loss_A_sub;
%     %   loss_A = loss_A + loss_A_sub;
%     
%     d_loss_S = -d_loss_A;
%     
%     d_loss_D1s{L_i} = 0;
%     d_loss_D2s{L_i} = 0;
%     for c = 1:3
%       d_loss_D1s{L_i} = d_loss_D1s{L_i} + d_loss_S(:,:,c) .* dS_Z{c}.F1;
%       d_loss_D2s{L_i} = d_loss_D2s{L_i} + d_loss_S(:,:,c) .* dS_Z{c}.F2;
%     end
%     
%     d_loss_Ls{L_i} = {};
%     for c = 1:3
%       d_loss_Ls{L_i}{c} = (reshape(d_loss_S(:,:,c), 1, []) * dS_L{c})';
%     end
%     d_loss_Ls{L_i} = cat(2, d_loss_Ls{L_i}{:});
%     
%     
%     [loss_L_sub, d_loss_L_sub, losses_L_sub(L_i)] = priorL(L, data, params);
%     losses_A(L_i) = losses_A(L_i) + loss_L_sub;
%     d_loss_Ls{L_i} = d_loss_Ls{L_i} + d_loss_L_sub;
% 
% %     mult = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS){1}/numel(L);
% %     if mult ~= 0
% %       [l, dl] = lmvnpdf2(L(:)', data.L_gaussian.mu, data.L_gaussian.Sigma, 1);
% %       losses_A(L_i) = losses_A(L_i) + mult*l;
% %       d_loss_Ls{L_i} = d_loss_Ls{L_i} + mult * reshape(dl, size(L));
% %     end
%     
%     %   mult = params.MULT_OPTS.srfs.light.GSM{1}/numel(L);
%     %   if mult ~= 0
%     %     [l, dl] = GSM_mvn_pdf(data.L_GSM, L(:)', 2);
%     %     losses_A(L_i) = losses_A(L_i) + mult*l;
%     %     d_loss_Ls{L_i} = d_loss_Ls{L_i} + mult * reshape(dl, size(L));
%     %   end
%     
%   end
%   
% end
% 
% % Sigma = cov(Ls_vec');
% % -log(det(Sigma))
% 
% alpha = params.LIGHT_ALPHA{1};
% if params.DO_LIGHT_POSTERIOR{1}
%   PA = global_Ls_weights .* exp(-alpha * (losses_A - min(losses_A)));
% else
%   PA = global_Ls_weights .* ones(length(losses_A),1);
% end
% PA = PA ./ sum(PA);
% 
% 
% loss_A = sum(losses_A .* PA);
% 
% if params.DO_LIGHT_POSTERIOR{1}
%   dPA = PA .* ( alpha*(loss_A - losses_A) + 1 );
% else
%   dPA = PA;
% end
% 
% d_loss_D1s = cat(3,d_loss_D1s{:});
% d_loss_D2s = cat(3,d_loss_D2s{:});
% 
% d_loss_D1 = sum(bsxfun(@times, d_loss_D1s, permute(dPA, [2,3,1])),3);
% d_loss_D2 = sum(bsxfun(@times, d_loss_D2s, permute(dPA, [2,3,1])),3);
% 
% d_loss_Ls = cellfun(@(x) x(:), d_loss_Ls, 'UniformOutput', false);
% d_loss_Ls = cat(2, d_loss_Ls{:});
% d_loss_Ls = bsxfun(@times, d_loss_Ls, dPA');
% % d_loss_Ls = d_loss_Ls * dShrink;
% 
% L_return = sum(bsxfun(@times, Ls, permute(PA, [2,3,1])),3);
% L = L_return;
% 
% d_loss_ZA = conv3d( d_loss_D1, dN_Z.f1) + conv3d( d_loss_D2, dN_Z.f2);
% 
% loss = loss_A + loss_Z;
% d_loss_Z = d_loss_ZA + d_loss_Z;
% 
% d_loss.Zpyr = buildGpyr_simple(d_loss_Z, length(state.Zpyr), params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% % d_loss.Zpyr = buildOGpyr(d_loss_Z, length(state.Zpyr)/5, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% 
% if params.DISABLE_ZPYR > 0
%   d_loss.Zpyr(1:params.DISABLE_ZPYR) = cellfun(@(x) zeros(size(x)), d_loss.Zpyr(1:params.DISABLE_ZPYR), 'UniformOutput', false);
% end
% 
% % if ~params.WHITEN_LIGHT_ENTROPY
% %   mult = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).entropy_mult{1}*27;
% %   if mult ~= 0
% %     X = reshape(Ls, 27, [])';
% %     [V, dV] = renyiK_brute(X, params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).entropy_sigma{1}, 2);
% %     loss = loss - mult*V;
% %     d_loss_Ls = d_loss_Ls - mult*dV';
% %   end
% % end
% 
% 
% d_loss_Ls_white = (d_loss_Ls' * data.L_whiten_params.inverse);
% 
% % if params.WHITEN_LIGHT_ENTROPY
% %   mult = params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).entropy_mult{1}*27;
% %   if mult ~= 0
% %     [V, dV] = renyiK_brute(Ls_white, params.MULT_OPTS.saifs.light.(params.L_WHITEN_PARAMS).entropy_sigma{1}, 2);
% %     loss = loss - mult*V;
% %     d_loss_Ls_white = d_loss_Ls_white - mult*dV;
% %   end
% % end
% 
% d_loss.Ls_white = d_loss_Ls_white / data.L_precond;
% % d_loss.Lvec = d_loss.Lvec / size(LCs,3);
% 
% if ~params.SOLVE_LIGHT
%   d_loss.Ls_white = zeros(size(d_loss.Ls_white));
% end
% 
% 
% to_disable = (PA < params.PA_EPSILON{1}) & (global_Ls_weights > 0);
% if any(to_disable)
%   
%   global_Ls_weights(to_disable) = 0;
% 
% end
% 
% Ls_vec = reshape(Ls, 27, [])';
% Dsq = distMat(Ls_vec, Ls_vec);
% 
% is_close = triu(double(global_Ls_weights > 0) * double(global_Ls_weights > 0)' & (Dsq < (params.L_D_MIN{1}.^2)), 1);
% if any(is_close(:))
%   
%   [i_min,i_max] = find(is_close);
%   if losses_A(i_max) < losses_A(i_min)
%     j = i_min;
%     i_min = i_max;
%     i_max = j;
%   end
%   
%   global_Ls_dependency(i_max) = i_min;
%   
%   global_Ls_weights(i_min) = global_Ls_weights(i_min) + global_Ls_weights(i_max);
%   global_Ls_weights(i_max) = 0;
%   
% end
% 
% % keyboard
% 
% % dispStruct('losses.height', losses.height)
% % dispStruct('losses.albedo', losses.albedo)
% % dispStruct('losses.light', losses.light)
% 
% % v1 = struct2vector(d_loss.Zpyr);
% % v2 = struct2vector(d_loss.Lwhite);
% % figure(70)
% % plot(1:length(v1), v1, 'b-'); hold on;
% % plot(length(v1) + [1:length(v2)], v2, 'r-'); hold off;
% % drawnow;
% 
% global global_loss_best
% if isempty(global_loss_best)
%   global_loss_best = inf;
% end
% 
% if (loss <= global_loss_best)
%   global_loss_best = loss;
% 
%   
%   
%   [vec, template] = struct2vector(losses_L_sub(find(PA > 0, 1, 'first')));
%   vec = vec * 0;
%   for li = find(PA > 0)'
%     vec = vec + PA(li)*struct2vector(losses_L_sub(li));
%   end
%   losses_L = vector2struct(vec, template);
%   
%   [vec, template] = struct2vector(losses_A_sub(find(PA > 0, 1, 'first')));
%   vec = vec * 0;
%   for li = find(PA > 0)'
%     vec = vec + PA(li)*struct2vector(losses_A_sub(li));
%   end
%   losses_A = vector2struct(vec, template);
%   
%   losses = [];
%   losses.albedo_smooth = losses_A.smooth;
%   losses.albedo_entropy = losses_A.entropy;
%   losses.albedo_hist = losses_A.hist;
%   
%   losses.shape_smooth = losses_Z.smooth;
%   losses.shape_contour = losses_Z.contour;
%   losses.shape_slant = losses_Z.slant;
%   
%   losses.light_gaussian = losses_L.(params.L_WHITEN_PARAMS).gaussian;
%   
%   global global_losses_history;
%   if isempty(global_losses_history)
%     global_losses_history = losses;
%   else
%     global_losses_history(end+1) = losses;
%   end
%   
%   
%   
%   if params.MAKE_VIDEO
%     
%     S = {};
%     for c = 1:3
%       S{c} = reshape(renderSH_helper(N_vec, L(:,c)), size(Z));
%     end
%     S = cat(3, S{:});
%     
%     A = data.true.log_im - S;
%     
%     invalid = repmat(~data.valid, [1,1,3]);
%     Z(~data.valid) = nan;
%     S(invalid) = nan;
%     A(invalid) = nan;
%     
%     I = data.true.im;
%     
%     Lv = visSH_color(L, [128, 192]);
%     
%     A = min(1, exp(A));
%     S = exp(S);
%     Lm = (1.25*max(Lv(:)));
%     S = S ./ Lm;
%     
%     A(isnan(A)) = 1;
%     S(isnan(S)) = 1;
%     
%     Zv = visualizeDEM(Z);
%     
%     Nv = visualizeNormals(Z);
%     
% %     keyboard
%     
%     Lv = Lv ./ Lm;
%     
%     Lv(isnan(Lv)) = 1;
%     Lv = max(0, min(1, Lv));
%     
%     Lvs = {};
%     for li = 1:size(Ls,3)
%       v = visSH_color(Ls(:,:,global_Ls_dependency(li)), [64,64]);
%       v(isnan(v)) = 1;
%       v = max(0, min(1, v));
%       Lvs{li} = v;
%     end
%         
% %     m = ((PA./max(PA)).^(1/4));% ./ max(PA);
% %     for li = 1:size(Ls,3)
% %       if global_Ls_weights(li) == 0
% %         Lvs{li} = ones(64,64,3);
% %       else
% %         v = visSH_color(Ls(:,:,li), [64,64]);
% %         v = v ./ (1.25*max(v(:)));
% %         v = 1 - (1-v)*m(li);
% %         v = max(0, min(1, v));
% %         Lvs{li} = v;
% %       end
% %     end
%     
% %     keyboard
%     
% %     Lvs = {};
% %     m = PA ./ max(PA);
% %     for li = 1:size(Ls,3)
% %       if global_Ls_weights(li) == 0
% %         Lvs{li} = zeros(64,64,3);
% %       else
% %         v = visSH_color(Ls(:,:,li), [64,64]);
% %         v = v ./ max(v(:));
% %         v(isnan(v)) = 1;
% %         v = 1 - (1-v) * m(li);
% %         Lvs{li} = v;
% %       end
% %     end
% %     
% %     for li = (size(Lvs,3)+1):(ceil(sqrt(size(Lvs,3)))^2)
% %       Lvs{li} = ones(size(v));
% %     end
% %     d = sqrt(length(Lvs));
% %     for li = 1:d
% %       Lvs{li} = cat(2,Lvs{(li-1)*d + [1:d]});
% %     end
% %     Lvs = cat(1,Lvs{1:d});
% %     Lvs = {};
% %     m = (PA.^(1/4));% ./ max(PA);
% %     for li = 1:size(Ls,3)
% % %       if global_Ls_weights(li) == 0
% % %         Lvs{li} = ones(64,64,3);
% % %       else
% %         v = visSH_color(Ls(:,:,li), [64,64]);
% %         iv = any(isnan(v),3);
% %         v1 = v(:,:,1);
% %         v2 = v(:,:,2);
% %         v3 = v(:,:,3);
% %         v1(iv) = 1;
% %         v2(iv) = 1-m(li);
% %         v3(iv) = 1-m(li);
% %         v = cat(3, v1, v2, v3);
% %         v = max(0, min(1, v));
% %         Lvs{li} = v;
% % %       end
% %     end
% 
%     for li = (size(Lvs,3)+1):(ceil(sqrt(size(Lvs,3)))^2)
%       Lvs{li} = zeros(size(v));
%     end
%     d = sqrt(length(Lvs));
%     for li = 1:d
%       Lvs{li} = cat(2,Lvs{(li-1)*d + [1:d]});
%     end
%     Lvs = cat(1,Lvs{1:d});
%     Lv = [Lv; Lvs];
% 
%     Lv = max(0, min(1, imresize(Lv, min(size(Z) ./ [size(Lv,1), size(Lv,2)]))));
%     
%     if size(Lv,1) < size(Z,1)
%       Lv = padarray(padarray(Lv, [floor((size(Z,1) - size(Lv,1))/2), 0], 1, 'pre'), [ceil((size(Z,1) - size(Lv,1))/2), 0], 1, 'post');
%     end
% 
%     if size(Lv,2) < size(Z,2)
%       Lv = padarray(padarray(Lv, [0, floor((size(Z,2) - size(Lv,2))/2)], 1, 'pre'), [0, ceil((size(Z,2) - size(Lv,2))/2)], 1, 'post');
%     end
%     
%     mat1 = mod([1:size(I,1)], 16) <= 8;
%     mat2 = mod([1:size(I,2)], 16) <= 8;
%     checker = bsxfun(@xor, mat2, mat1');
%     I_bg = repmat(1-checker*0.2, [1,1,3]);
%     I(invalid) = I_bg(invalid);
%     
%     I(1:4,:,1) = 0.75;
%     I(:,1:4,1) = 0.75;
%     I(end-3:end,:,1) = 0.75;
%     I(:,end-3:end,1) = 0.75;
%     
%     I(1:4,:,2) = 0.25;
%     I(:,1:4,2) = 0.25;
%     I(end-3:end,:,2) = 0.25;
%     I(:,end-3:end,2) = 0.25;
%     
%     I(1:4,:,3) = 0.25;
%     I(:,1:4,3) = 0.25;
%     I(:,end-3:end,3) = 0.25;
%     I(end-3:end,:,3) = 0.25;
%     
% %     I(invalid) = 1;
% %     I1 = I(:,:,1);
% %     I2 = I(:,:,2);
% %     I3 = I(:,:,3);
% %     I1(~data.valid) = 0.75;
% %     I2(~data.valid) = 0.25;
% %     I3(~data.valid) = 0.25;
% %     I = cat(3, I1, I2, I3);
%     
%     V1 = [I, Zv, Nv];
%     V2 = [A, S, Lv];
%     V = [V1; V2];
%     
%     global GLOBAL_VIDEO_COUNT;
%     imwrite(uint8(round(255*V)), [data.video_folder, 'frame_', num2str(GLOBAL_VIDEO_COUNT, '%05d'), '.jpg'], 'Quality', 100)
% %     if GLOBAL_VIDEO_COUNT == 1
% %       
% %       imwrite(I, [data.video_folder, 'image.png'])
% %     end
%     
% %     fprintf(data.video_log, '%e ', [loss_A, loss_Z, loss_L]);
%     fprintf(data.video_log, '%e ', [loss]);
%     fprintf(data.video_log, '\n');
%     GLOBAL_VIDEO_COUNT = GLOBAL_VIDEO_COUNT + 1;
%     
%   else
%   
%     if params.DO_DISPLAY
%       
%       global last_display
%       if isempty(last_display) || (etime(clock, last_display) > params.DISPLAY_PERIOD)
%         last_display = clock;
%         
% %         keyboard
%         trace = {};
%         [trace{1}, template] = struct2vector(global_losses_history(1));
%         for i = 2:length(global_losses_history)
%           trace{i} = struct2vector(global_losses_history(i));
%         end
%         trace = cat(2,trace{:});
%         
%         for i = 1:length(template.fields)
%           template.fields{i}(template.fields{i} == '_') = '.';
%         end
%                 
% %         if length(global_losses_history) > 1
% %           figure(3);
% %           plot(max(trace', 0.01))
% %           set(gca, 'Yscale', 'log');
% %           legend(template.fields)
% %         end
%         
%         
%         S = {};
%         for c = 1:3          
%           S{c} = reshape(renderSH_helper(N_vec, L(:,c)), size(Z));
%         end
%         S = cat(3, S{:});
%         
%         A = data.true.log_im - S;
%         
%         invalid = repmat(~data.valid, [1,1,3]);
%         Z(~data.valid) = nan;
%         S(invalid) = nan;
%         A(invalid) = nan;
%         
%         I = data.true.im;
%         
%         A = exp(A);
%         S = exp(S);
%         
%         %       A = max(0, min(1, A));
%         %       S = max(0, min(1, S));
%         
%         state = struct('normal', N, 'height', Z, 'albedo', A, 'shading', S, 'light', L);
%         clear N Z A S L
%         
%         Lvs = {};
% %         PA
%         m = PA ./ max(PA);
%         
%         for li = 1:size(Ls,3)
%           if global_Ls_weights(li) == 0
%             Lvs{li} = zeros(64,64,3);
%           else
%             v = visSH_color(Ls(:,:,global_Ls_dependency(li)), [64,64]);
%             iv = any(isnan(v),3);
%             v(isnan(v)) = 0;
%             v1 = v(:,:,1);
%             v2 = v(:,:,2);
%             v3 = v(:,:,3);
%             v1(iv) = m(li);
%             v2(iv) = 0;
%             v3(iv) = 0;
%             v = cat(3, v1, v2, v3);
%             v = max(0, min(1, v));
%             Lvs{li} = v;
%           end
%         end
%         
%         for li = (length(Lvs)+1):9%(ceil(sqrt(size(Lvs,3)))^2)
%           Lvs{li} = zeros(size(v));
%         end
%         d = sqrt(length(Lvs));
%         for li = 1:d
%           Lvs{li} = cat(2,Lvs{(li-1)*d + [1:d]});
%         end
%         Lvs = cat(1,Lvs{1:d});
%         Lvs = max(0, min(1, imresize(Lvs, [nan, 150])));
%         
%         if isfield(data.true, 'height')
%           
% %           [err, state] = getError(state, data.true);
%           [err, junk] = getError(state, data.true);
%           
%           Lv = visSH_color(state.light, [150, 150]);
%           Lv(isnan(Lv)) = 0;
%           Lv = [Lv; Lvs];
%           
%           s = size(state.height,1);
%           if size(Lv,1) > s
%             Lv = max(0, min(1, imresize(Lv, [s, 150])));
%           else
%             pad1 = floor((s - size(Lv,1))/2);
%             pad2 = ceil((s - size(Lv,1))/2);
%             Lv = padarray(padarray(Lv, [pad2], 0, 'post'), [pad1], 0, 'pre');
%           end
%           
% %           m = max(Lv(:));
% %           Lv = Lv ./ m;
% %           state.shading = state.shading ./ m;
% %           state.albedo = state.albedo .* m;
%           Lv = max(0, min(1, Lv));
%           
%           if size(data.true.light,2) == 1
%             data.true.light = repmat(data.true.light, [1,3]);
%           end
%           Lv_true = visSH_color(data.true.light, [size(state.height,1), 150]);
%           Lv_true(isnan(Lv_true)) = 0;
%           Lv_true = max(0, min(1, Lv_true));
%           
%           I(invalid) = 0;
%           data.true.shading(invalid) = 0;
%           Lv(isnan(Lv)) = 0;
%           state.albedo(invalid) = 0;
%           state.shading(invalid) = 0;
%           
% %           mat1 = mod([1:size(I,1)], 16) <= 8;
% %           mat2 = mod([1:size(I,2)], 16) <= 8;
% %           checker = bsxfun(@xor, mat2, mat1');
% %           I_bg = repmat(1-checker*0.2, [1,1,3]);
% %           I(invalid) = I_bg(invalid);
%           
%           Z = state.height;
%           shift = mean(Z(~isnan(data.true.height)) - data.true.height(~isnan(data.true.height)));
%           Z = Z - shift;
%           Zv = visualizeDEM([Z; data.true.height]);
% %           Zv = [visualizeDEM(state.height); visualizeDEM(data.true.height)];
% %           Zv = visualizeDEM([state.height; data.true.height]);
%           Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
%           Nv = visualizeNormals_color(state.normal);
%           Ntv = visualizeNormals_color(data.true.normal);
%           Ntv(isnan(data.true.normal)) = 0;
%           Nv(invalid) = 0;
%           I_pad = padarray(padarray(I, floor(size(I,1)/2), 0, 'pre'), ceil(size(I,1)/2), 0, 'post'); % hah, ipad.
%           V = min(1, [I_pad, Zv, [Nv; Ntv], [state.albedo; data.true.albedo], [state.shading; data.true.shading], [Lv; Lv_true]]);
%           
%         else
%           
%           Lv = visSH_color(state.light, [size(state.height,1), 150]);
%           Lv(isnan(Lv)) = 0;
%           Lv = max(0, min(1, Lv));
%           
%           I(invalid) = 0;
%           Lv(isnan(Lv)) = 0;
%           state.albedo(invalid) = 0;
%           state.shading(invalid) = 0;
%           
%           state.albedo = state.albedo ./ max(state.albedo(:));
%           state.shading = state.shading ./ max(state.shading(:));
%           Lv = Lv ./ max(Lv(:));
%           
%           Zv = visualizeDEM([state.height]);
%           Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
%           Nv = visualizeNormals_color(state.normal);
%           Nv(invalid) = 0;
%           V = min(1, [I, Zv, Nv, state.albedo, state.shading, Lv]);
%           
%         end
%         
%         figure(1);
%         imagesc(V);
%         axis image off;
%         set(gca, 'PlotBoxAspectRatio', [1 1 1])
%         set(gca, 'Position', [0 0 1 1]);
%         set(gcf, 'Position', [1, 500, 1600, 1600/size(V,2)*size(V,1)])
%         
%         if isfield(data.true, 'height')
%           %     figure(3)
%           %     a = flattenPyr(Apyr);
%           %     imagesc(min(1,exp(a))); imtight;
%           %     set(gcf, 'Position', [1, 1, 1600, 1600/size(V,2)*size(V,1)])
%           %     drawnow;
%           
%           global global_errors;
%           global_errors{end+1} = err;
%           light_errs = cellfun(@(x) x.light, global_errors);
% %           curvature_errs = cellfun(@(x) x.curvature, global_errors);
%           normal_errs = cellfun(@(x) x.normal, global_errors);
%           %     survey_errs = cellfun(@(x) x.survey, global_errors);
%           albedo_errs = cellfun(@(x) x.albedo, global_errors);
%           shading_errs = cellfun(@(x) x.shading, global_errors);
%           grosse_errs = cellfun(@(x) x.grosse, global_errors);
%           avg_errs = cellfun(@(x) x.avg, global_errors);
%           
%           
%           figure(2);
%           plot(light_errs ./ max(eps,light_errs(1)), 'Color', [0.5, 0.5, 0.5]); hold on;
% %           plot(curvature_errs ./ max(eps,curvature_errs(1)), 'Color', [0.5, 0.5, 0.5]); hold on;
%           plot(normal_errs ./ max(eps,normal_errs (1)), 'k-');
%           
%           plot(grosse_errs ./ max(eps,grosse_errs(1)), 'b-');
%           
%           plot(albedo_errs ./ max(eps,albedo_errs(1)), 'r-');
%           plot(shading_errs ./ max(eps,shading_errs(1)), 'g-');
%           
%           plot(avg_errs ./ max(eps,avg_errs(1)), 'm-');
%           
%           axis square
%           set(gca, 'YLim', [.009, 2])
%           set(gca, 'YScale', 'log');
%           title(['gray = LightMSE, black = NMSE, b = LMSE, r = AMSE, g = SMSE, m = avg'])
% %           title(['gray = KMSE, black = NMSE, b = LMSE, r = AMSE, g = SMSE, m = avg'])
%           grid on;
%           hold off;
%         end
%         drawnow;
%         
%       end
%     end
%   end
% end




% function [loss, d_loss, L_expected, L_max, L_prob] =  lossfun_srfs_color(state, data, params)
% 
% state.height = upsample(state.Zdown, data.pind);
% 
% alpha = params.ALPHA_LIGHT{1};
% 
% Z = state.height;
% 
% [N, dN_Z] = getNormals_conv(Z);
% N_vec = reshape(N, [], 3);
% 
% [loss_Z, d_loss_Z] = priorZ(Z, data, params);
% [loss_Z_simple, d_loss_Z_simple] = priorZ_simple(Z, data, params, N, dN_Z);
% 
% loss_Z = alpha * loss_Z;
% d_loss_Z = alpha * d_loss_Z;
% loss_Z_simple = alpha * loss_Z_simple;
% d_loss_Z_simple = alpha * d_loss_Z_simple;
% 
% L_meta = data.L_meta;
% 
% Ls = data.true.light;
% 
% global global_L_active_idx;
% 
% d_loss_N = 0;
% d_losses_N = {};
% 
% % keyboard
% 
% if params.USE_GREEDY
% 
%   global global_in_line_search
%   if isempty(global_in_line_search)
%     global_in_line_search = false;
%   end
%   
%   losses_A = (10^10) * ones(1, size(Ls,3));
%   
% %   if ~global_in_line_search
% % %     fprintf('I''m just starting out, fingers crossed!\n');
%     queue = [L_meta.head, -inf];
%     ADD_CHILDREN = true;
% %   else
% % %     fprintf('Im in a line search!\n');
% %     queue = [global_L_active_idx', [1:length(global_L_active_idx)]'];
% %     ADD_CHILDREN = false;
% %   end
%   
%   
%   loss_A_min = [];
%   
%   closed_set = false(size(Ls,3),1);
%   for i_search = 1:params.GREEDY_N_SEARCH
%     
%     if isempty(queue)
%       break
%     end
%     
%     [junk, min_i] = min(queue(:,2));
%     i = queue(min_i,1);
%     queue = queue([[1:(min_i-1)], [(min_i+1):size(queue,1)]],:);
%     closed_set(i) = true;
%     L = Ls(:,:,i);
%     
%     S = {};
%     dS_N = {};
%     for c = 1:3
%       [s, ds_n] = renderSH_helper(N_vec, L(:,c));
%       S{c} = reshape(s, size(Z));
%       dS_N{c} = reshape(ds_n, [size(Z), 3]);
%     end
%     S = cat(3, S{:});
%     
%     A = data.true.log_im - S;
%     
%     [loss_A, d_loss_A] = priorA(A, data, params);
%     loss_A = alpha * loss_A;
%     d_loss_A = alpha * d_loss_A;
% 
%     losses_A(i) = loss_A;
%     
%     d_loss_Ni = 0;
%     for c = 1:3
%       d_loss_Ni = d_loss_Ni + repmat(-d_loss_A(:,:,c), [1,1,3]) .* dS_N{c};
%     end
%     d_losses_N{i} = d_loss_Ni;
% 
%     if ~isempty(loss_A_min)
%       
%       add_children = loss_A < (loss_A_min + params.GREEDY_LL_MIN);
%       loss_A_min = min(loss_A_min, loss_A);
%       
%     else
%       
%       loss_A_min = loss_A;
%       add_children = true;
%       
%     end
%     
%     if ADD_CHILDREN && add_children
%       kids = find((L_meta.children(:,i) | L_meta.children(i,:)') & ~closed_set);
%       %     kids = find(L_meta.children(:,i));
%       kids_depth = L_meta.depth(kids);
%       queue = [queue; [kids, loss_A + kids_depth*params.GREEDY_DEPTH_COST*alpha]];
%     end
%     
%   end
%   
%   L_active_idx = find(losses_A  <  10^10);
%   global_L_active_idx = L_active_idx;
%   
%   LL = -losses_A;
%   LL = LL - max(LL);
%   P = exp(LL);
%   P = P ./ sum(P);
%   
%   loss_A = sum(losses_A .* P);
%   L_expected = sum(Ls .* repmat(reshape(P, [1,1, size(P,2)]), size(Ls,1), size(Ls,2)),3);
%   [junk, max_idx] = max(P);
%   L_max = Ls(:,:,max_idx);
%   L_prob = P;
%   
%   
%   d_loss_N = 0;
%   
%   dLs = P .* ( (loss_A - losses_A) + 1 );
%   for i = L_active_idx
%     d_loss_N = d_loss_N + dLs(i) * d_losses_N{i};
%   end
% 
% %   for i = L_active_idx
% %     d_loss_N = d_loss_N + P(i) * d_losses_N{i};
% %   end
% 
% else
%   
%   L_active_idx = find(data.L_active);
% 
%   losses_A = (10^10) * ones(1, size(Ls,3));
%   LL_online = -inf(1, size(Ls,3));
%   loss_online = 0;
%   for i = L_active_idx
%     L = Ls(:,:,i);
%     
%     S = {};
%     dS_N = {};
%     for c = 1:3
%       [s, ds_n] = renderSH_helper(N_vec, L(:,c));
%       S{c} = reshape(s, size(Z));
%       dS_N{c} = reshape(ds_n, [size(Z), 3]);
%     end
%     S = cat(3, S{:});
%     
%     A = data.true.log_im - S;
%     
%     [loss_A, d_loss_A] = priorA(A, data, params);
%     losses_A(i) = loss_A;
%     
%     d_loss_Ni = 0;
%     for c = 1:3
%       d_loss_Ni = d_loss_Ni + repmat(-d_loss_A(:,:,c), [1,1,3]) .* dS_N{c};
%     end
%     
%     if params.USE_APPROX_LIGHT_GRADIENT
%       
%       LL_online(i) = -alpha*loss_A;
%       P_online = exp(LL_online - max(LL_online));
%       P_online = P_online ./ sum(P_online);
%       
%       loss_online = P_online(i).* loss_A + (1-P_online(i)) * loss_online;
%       
%       d_loss_N = P_online(i).* d_loss_Ni + (1-P_online(i)) * d_loss_N;
%     else
%       d_losses_N{i} = d_loss_Ni;
%     end
%     
%   end
%   
%   LL = -alpha * losses_A;
%   LL = LL - max(LL);
%   P = exp(LL);
%   P = P ./ sum(P);
%   
%   loss_A = sum(losses_A .* P);
%   L_expected = sum(Ls .* repmat(reshape(P, [1,1, size(P,2)]), size(Ls,1), size(Ls,2)),3);
%   L_prob = P;
%   
%   
%   if ~params.USE_APPROX_LIGHT_GRADIENT
%     dLs = P .* ( alpha*(loss_A - losses_A) + 1 );
%     for i = L_active_idx
%       d_loss_N = d_loss_N + dLs(i) * d_losses_N{i};
%     end
%   end
% 
% end
% 
% d_loss_ZA = ...
%   + conv3d( dN_Z.F1_1 .* d_loss_N(:,:,1) + dN_Z.F1_2 .* d_loss_N(:,:,2) + dN_Z.F1_3 .* d_loss_N(:,:,3), dN_Z.f1) ...
%   + conv3d( dN_Z.F2_1 .* d_loss_N(:,:,1) + dN_Z.F2_2 .* d_loss_N(:,:,2) + dN_Z.F2_3 .* d_loss_N(:,:,3), dN_Z.f2);
% 
% 
% loss = loss_A + loss_Z + loss_Z_simple;
% d_loss_Z = d_loss_Z + d_loss_ZA + d_loss_Z_simple;
% 
% loss = loss / alpha;
% d_loss_Z = d_loss_Z / alpha;
% 
% d_loss.Zdown = downsample(d_loss_Z, size(data.pind,1), namedFilter(PYR_FILTER_NAME)*2);
% % d_loss.Zdown = downsample(d_loss_Z, size(data.pind,1), [1;4;6;4;1]/8);
% 
% 
% global global_loss_best
% if isempty(global_loss_best)
%   global_loss_best = inf;
% end
% 
% if (loss <= global_loss_best)
%   global_loss_best = loss;
%   
% %   global global_greedy_head
% %   [m, m_idx] = max(LL);
% %   if m > (log(8) + LL(global_greedy_head))
% %     global_greedy_head = m_idx;
% %     if params.GREEDY_MOVE_HEAD
% %       fprintf('New Head: %d\n', global_greedy_head);
% %     end
% %   end
% 
%   global global_P global_losses_A
%   global_P = P;
%   global_losses_A = losses_A;
%   
%   losses_A(losses_A == 10^10) = nan;
%   
%   
%   global global_light_fid;
%   fprintf(global_light_fid, '%e ', losses_A);
%   fprintf(global_light_fid, '\n');
% 
%   
%   if params.DO_DISPLAY
%     
%     global last_display
%     if isempty(last_display) || (etime(clock, last_display) > params.DISPLAY_PERIOD)
%       last_display = clock;
%       
%       invalid = repmat(~data.valid, [1,1,3]);
%       Z(~data.valid) = nan;
%       
%       I = data.true.im;
%       
%       L = visSH_color(L_expected, size(Z));
%       %     L = visSH_color(data.true.light, size(Z));
%       L = L ./ max(L(:));
% 
%       S = {};
%       for c = 1:3
%         S{c} = reshape(renderSH_helper(N_vec, L_expected(:,c)), size(Z));
%       end
%       S = cat(3,S{:});
% 
%       A = data.true.log_im - S;
%       
%       A = exp(A);
%       S = exp(S);
%       
%       A = max(0, min(1, A));
%       S = max(0, min(1, S));
%       
%       S(invalid) = nan;
%       A(invalid) = nan;
% 
%       
%       state = struct('normal', N, 'height', Z, 'albedo', A, 'shading', S);
%       
%       [err, state] = getError(state, data.true);
%       
%       I(invalid) = 0;
%       data.true.shading(invalid) = 0;
%       L(isnan(L)) = 0;
%       state.albedo(invalid) = 0;
%       state.shading(invalid) = 0;
%       
%       Zv = visualizeDEM([state.height; data.true.height]);
%       Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
%       Nv = visualizeNormals_color(N);
%       Ntv = visualizeNormals_color(data.true.normal);
%       Ntv(isnan(data.true.normal)) = 0;
%       Nv(isnan(data.true.normal)) = 0;
%       V = min(1, [Zv, [Nv; Ntv], [state.shading; data.true.shading], [state.albedo; data.true.albedo], [L; I]]);
%       
%       figure(1);
%       imagesc(V);
%       axis image off;
%       set(gca, 'PlotBoxAspectRatio', [1 1 1])
%       set(gca, 'Position', [0 0 1 1]);
%       set(gcf, 'Position', [1, 500, 1600, 1600/size(V,2)*size(V,1)])
%       %       set(gcf, 'Position', get(gcf, 'Position') .* [1 1 0 0] + [0 0 size(V,2) size(V,1)]);
%       
%       
%       global global_errors;
%       global_errors{end+1} = err;
%       height_errs = cellfun(@(x) x.height, global_errors);
%       %     affine_height_errs = cellfun(@(x) x.affine_height, global_errors);
%       survey_errs = cellfun(@(x) x.survey, global_errors);
%       albedo_errs = cellfun(@(x) x.albedo, global_errors);
%       shading_errs = cellfun(@(x) x.shading, global_errors);
%       grosse_errs = cellfun(@(x) x.grosse, global_errors);
%       
%       figure(2);
%       plot(height_errs ./ max(eps,height_errs(1)), 'b-'); hold on;
%       %     plot(affine_height_errs ./ max(eps,affine_height_errs(1)), 'Color', [0.5, 0.5, 1]);
%       plot(survey_errs ./ max(eps,survey_errs(1)), 'r-');
%       plot(grosse_errs ./ max(eps,grosse_errs(1)), 'k-');
%       plot(albedo_errs ./ max(eps,albedo_errs(1)), 'g-');
%       plot(shading_errs ./ max(eps,shading_errs(1)), 'm-');
%       axis square
%       ylim = get(gca, 'YLim');
%       set(gca, 'YLim', [.09, 1.1])
%       set(gca, 'YScale', 'log');
%       title(['b = Z(', num2str(height_errs(end)),'), r = I(', num2str(survey_errs(end)),'), k = grosse(', num2str(grosse_errs(end)),'), g = A, m = S, m = loss'])
%       hold off;
%       
%       X = L_meta.X;
%       Y = L_meta.Y;
%       children = L_meta.children;
%       V = L_meta.vis;
%       
%       children = children(L_active_idx, L_active_idx);
%       X = X(L_active_idx);
%       Y = Y(L_active_idx);
%       V = V(L_active_idx);
%       
%       [i,j] = find(children);
%       
%       figure(4);
%       clf
% 
%       line(X([i,j])', Y([i,j])', 'color', 'k')
%       axis equal off
%       P_active = P(L_active_idx);
% %       sort(P_active, 'descend')
%       Q = log(P_active ./ max(P_active));
%       Q = max(0, min(1, (Q + 5)/5));
%       
%       hold on;
%       for vi = 1:length(V)
%         v = V{vi};
%         v = cat(3, padarray(v(:,:,1), [3,3], 1), padarray(v(:,:,2:3), [3,3], 1-Q(vi)));
%         imagesc(X(vi) - size(V{vi},1)/2, Y(vi) - size(V{vi},1)/2, v );
%       end
%       hold off;
%       
%       drawnow;
%       
%     end
%   end
%   
% end


% function [loss, d_loss, L_expected, L_prob] =  lossfun_srfs_color(state, data, params)
% 
% state.height = upsample(state.Zdown, data.pind);
% 
% Z = state.height;
% 
% [N, dN_Z] = getNormals_conv(Z);
% N_vec = reshape(N, [], 3);
% 
% [loss_Z, d_loss_Z] = priorZ(Z, data, params);
% [loss_Z_simple, d_loss_Z_simple] = priorZ_simple(Z, data, params, N, dN_Z);
% 
% 
% Ls = data.true.light;
% global global_L_active;
% global global_losses_A
% if isempty(global_L_active)
%   global_L_active = true(1,size(Ls,3));
%   global_losses_A = inf(1, size(Ls,3));
%   
% %   global_L_active = rand(size(global_L_active))>0.5;
% %   global_losses_A = randn(size(global_losses_A))*4 + 62;
%   
% end
% 
% L_active_idx = find(global_L_active);
% 
% losses_A = global_losses_A;
% 
% alpha = params.ALPHA_LIGHT{1};
% 
% d_loss_N = 0;
% d_losses_N = {};
% 
% LL_online = -inf(1, length(L_active_idx));
% loss_online = 0;
% for i = L_active_idx
%   L = Ls(:,:,i);
%   
%   S = {};
%   dS_N = {};
%   for c = 1:3
%     [s, ds_n] = renderSH_helper(N_vec, L(:,c));
%     S{c} = reshape(s, size(Z));
%     dS_N{c} = reshape(ds_n, [size(Z), 3]);
%   end
%   S = cat(3, S{:});
%   
%   A = data.true.log_im - S;
% 
%   [loss_A, d_loss_A] = priorA(A, data, params);
%   losses_A(i) = loss_A;
%   
%   d_loss_Ni = 0;
%   for c = 1:3
%     d_loss_Ni = d_loss_Ni + repmat(-d_loss_A(:,:,c), [1,1,3]) .* dS_N{c};
%   end
%   
%   if params.USE_APPROX_LIGHT_GRADIENT
%     
%     LL_online(i) = -alpha*loss_A;
%     P_online = exp(LL_online - max(LL_online));
%     P_online = P_online ./ sum(P_online);
%     
%     loss_online = P_online(i).* loss_A + (1-P_online(i)) * loss_online;
%     
%     d_loss_N = P_online(i).* d_loss_Ni + (1-P_online(i)) * d_loss_N;
%   else
%     d_losses_N{i} = d_loss_Ni;
%   end
%   
% end
% 
% LL = -alpha * losses_A;
% LL = LL - max(LL);
% P = exp(LL);
% P = P ./ sum(P);
% 
% loss_A = sum(losses_A .* P);
% L_expected = sum(Ls .* repmat(reshape(P, [1,1, size(P,2)]), size(Ls,1), size(Ls,2)),3);
% L_prob = P;
% 
% 
% if ~params.USE_APPROX_LIGHT_GRADIENT
%   dLs = P .* ( alpha*(loss_A - losses_A) + 1 );
%   for i = L_active_idx
%     d_loss_N = d_loss_N + dLs(i) * d_losses_N{i};
%   end
% end
% 
% 
% d_loss_ZA = ...
%   + conv3d( dN_Z.F1_1 .* d_loss_N(:,:,1) + dN_Z.F1_2 .* d_loss_N(:,:,2) + dN_Z.F1_3 .* d_loss_N(:,:,3), dN_Z.f1) ...
%   + conv3d( dN_Z.F2_1 .* d_loss_N(:,:,1) + dN_Z.F2_2 .* d_loss_N(:,:,2) + dN_Z.F2_3 .* d_loss_N(:,:,3), dN_Z.f2);
% 
% 
% loss = loss_A + loss_Z + loss_Z_simple;
% d_loss_Z = d_loss_Z + d_loss_ZA + d_loss_Z_simple;
% 
% d_loss.Zdown = downsample(d_loss_Z, size(data.pind,1), [1;4;6;4;1]/8);
% 
% 
% global global_loss_best
% if isempty(global_loss_best)
%   global_loss_best = inf;
% end
% 
% 
% if (loss < global_loss_best)
%   global_loss_best = loss;
% 
%   below = LL < -params.LIGHT_LL_CUTOFF{1};
%   global_L_active(below) = false;
%   global_losses_A(below) = losses_A(below);
%   
%   global global_light_fid;
%   fprintf(global_light_fid, '%e ', losses_A);
%   fprintf(global_light_fid, '\n');
% end
% 
% 
% if params.DO_DISPLAY
%   
%   global last_display
%   if isempty(last_display) || (etime(clock, last_display) > params.DISPLAY_PERIOD)
%     last_display = clock;
%     
%     invalid = repmat(~data.valid, [1,1,3]);
%     Z(~data.valid) = nan;
%     S(invalid) = nan;
%     A(invalid) = nan;
%     
%     I = data.true.im;
%     
%     L = visSH_color(L_expected, size(Z));
% %     L = visSH_color(data.true.light, size(Z));
%     L = L ./ max(L(:));
%     
%     A = exp(A);
%     S = exp(S);
%     A = A ./ max(A(:));
%     S = S ./ max(S(:));
%     
%     state = struct('normal', N, 'height', Z, 'albedo', A, 'shading', S);
%     
%     [err, state] = getError(state, data.true);
%     
%     I(invalid) = 0;
%     data.true.shading(invalid) = 0;
%     L(isnan(L)) = 0;
%     state.albedo(invalid) = 0;
%     state.shading(invalid) = 0;
%     
%     Zv = visualizeDEM([state.height; data.true.height]);
%     Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
%     Nv = visualizeNormals_color(N);
%     Ntv = visualizeNormals_color(data.true.normal);
%     Ntv(isnan(data.true.normal)) = 0;
%     Nv(isnan(data.true.normal)) = 0;
%     V = min(1, [Zv, [Nv; Ntv], [state.shading; data.true.shading], [state.albedo; data.true.albedo], [L; I]]);
%     
%     figure(1);
%     imagesc(V);
%     axis image off;
%     set(gca, 'PlotBoxAspectRatio', [1 1 1])
%     set(gca, 'Position', [0 0 1 1]);
%     set(gcf, 'Position', [1, 500, 1600, 1600/size(V,2)*size(V,1)])
%     %       set(gcf, 'Position', get(gcf, 'Position') .* [1 1 0 0] + [0 0 size(V,2) size(V,1)]);
%     
%     
%     global global_errors;
%     global_errors{end+1} = err;
%     height_errs = cellfun(@(x) x.height, global_errors);
% %     affine_height_errs = cellfun(@(x) x.affine_height, global_errors);
%     survey_errs = cellfun(@(x) x.survey, global_errors);
%     albedo_errs = cellfun(@(x) x.albedo, global_errors);
%     shading_errs = cellfun(@(x) x.shading, global_errors);
%     grosse_errs = cellfun(@(x) x.grosse, global_errors);
%     
%     figure(2);
%     plot(height_errs ./ max(eps,height_errs(1)), 'b-'); hold on;
% %     plot(affine_height_errs ./ max(eps,affine_height_errs(1)), 'Color', [0.5, 0.5, 1]);
%     plot(survey_errs ./ max(eps,survey_errs(1)), 'r-');
%     plot(grosse_errs ./ max(eps,grosse_errs(1)), 'k-');
%     plot(albedo_errs ./ max(eps,albedo_errs(1)), 'g-');
%     plot(shading_errs ./ max(eps,shading_errs(1)), 'm-');
%     axis square
%     ylim = get(gca, 'YLim');
%     set(gca, 'YLim', [.09, 1.1])
%     set(gca, 'YScale', 'log');
%     title(['b = Z(', num2str(height_errs(end)),'), r = I(', num2str(survey_errs(end)),'), k = grosse(', num2str(grosse_errs(end)),'), g = A, m = S, m = loss, ', num2str(sum(global_L_active)), ' L active'])
%     hold off;
%     
%     drawnow;
%     
%   end
% end
% 
% 
% 
% % function [loss, d_loss, L_expected, L_prob] =  lossfun_srfs_color(state, data, params)
% % 
% % state.height = upsample(state.Zdown, data.pind);
% % 
% % Z = state.height;
% % 
% % [N, dN_Z] = getNormals_conv(Z);
% % N_vec = reshape(N, [], 3);
% % 
% % [loss_Z, d_loss_Z] = priorZ(Z, data, params);
% % [loss_Z_simple, d_loss_Z_simple] = priorZ_simple(Z, data, params, N, dN_Z);
% % 
% % 
% % Ls = data.true.light;
% % global global_L_active;
% % global global_losses_A
% % if isempty(global_L_active)
% %   global_L_active = true(1,size(Ls,3));
% %   global_losses_A = inf(1, size(Ls,3));
% %   
% % %   global_L_active = rand(size(global_L_active))>0.5;
% % %   global_losses_A = randn(size(global_losses_A))*4 + 62;
% %   
% % end
% % 
% % L_active_idx = find(global_L_active);
% % 
% % 
% % losses_A = global_losses_A;
% % d_losses_N = {};
% % 
% % alpha = params.ALPHA_LIGHT{1};
% % 
% % % LL_max = -inf;
% % % LL_online = -inf(1, length(L_active_idx));
% % % L_online = 0;
% % for i = L_active_idx
% %   L = Ls(:,:,i);
% %   
% %   S = {};
% %   dS_N = {};
% %   for c = 1:3
% %     [s, ds_n] = renderSH_helper(N_vec, L(:,c));
% %     S{c} = reshape(s, size(Z));
% %     dS_N{c} = reshape(ds_n, [size(Z), 3]);
% %   end
% %   S = cat(3, S{:});
% %   
% %   A = data.true.log_im - S;
% % 
% %   [loss_A, d_loss_A] = priorA(A, data, params);
% %   losses_A(i) = loss_A;
% %   
% % %   LL_new = -loss_A;
% % %   
% % %   if i == 1
% % %     LL_online(i) = 0;
% % %     LL_max = LL_new;
% % %   else
% % %     if LL_new > LL_max
% % %       LL_online(1:(i-1)) = LL_online(1:(i-1)) + LL_max - LL_new;
% % %       LL_online(i) = 0;
% % %       LL_max = LL_new;
% % %     else
% % %       LL_online(i) = LL_new - LL_max;
% % %     end
% % %   end
% % 
% % %   P = exp(LL_online);
% % %   P_sum = sum(P_online);
% % %   P = P ./ P_sum;
% % %   
% % %   L_online = L_online + 
% %   
% %   d_losses_N{i} = 0;
% %   for c = 1:3
% %     d_losses_N{i} = d_losses_N{i} + repmat(-d_loss_A(:,:,c), [1,1,3]) .* dS_N{c};
% %   end
% %   
% % end
% % 
% % LL = -alpha * losses_A;
% % LL = LL - max(LL);
% % P = exp(LL);
% % P = P ./ sum(P);
% % 
% % loss_A = sum(losses_A .* P);
% % L_expected = sum(Ls .* repmat(reshape(P, [1,1, size(P,2)]), size(Ls,1), size(Ls,2)),3);
% % L_prob = P;
% % 
% % 
% % dLs = P .* ( alpha*(loss_A - losses_A) + 1 );
% % d_loss_N = 0;
% % for i = L_active_idx
% %   d_loss_N = d_loss_N + dLs(i) * d_losses_N{i};
% % end
% % 
% % d_loss_ZA = ...
% %   + conv3d( dN_Z.F1_1 .* d_loss_N(:,:,1) + dN_Z.F1_2 .* d_loss_N(:,:,2) + dN_Z.F1_3 .* d_loss_N(:,:,3), dN_Z.f1) ...
% %   + conv3d( dN_Z.F2_1 .* d_loss_N(:,:,1) + dN_Z.F2_2 .* d_loss_N(:,:,2) + dN_Z.F2_3 .* d_loss_N(:,:,3), dN_Z.f2);
% % 
% % 
% % loss = loss_A + loss_Z + loss_Z_simple;
% % d_loss_Z = d_loss_Z + d_loss_ZA + d_loss_Z_simple;
% % 
% % d_loss.Zdown = downsample(d_loss_Z, size(data.pind,1), [1;4;6;4;1]/8);
% % 
% % 
% % global global_loss_best
% % if isempty(global_loss_best)
% %   global_loss_best = inf;
% % end
% % 
% % 
% % if (loss < global_loss_best)
% %   global_loss_best = loss;
% % 
% %   L = alpha*LL;
% %   below = L < -params.LIGHT_LL_CUTOFF{1};
% %   global_L_active(below) = false;
% %   global_losses_A(below) = losses_A(below);
% %   
% %   global global_light_fid;
% %   fprintf(global_light_fid, '%e ', losses_A);
% %   fprintf(global_light_fid, '\n');
% % end
% % 
% % 
% % if params.DO_DISPLAY
% %   
% %   global last_display
% %   if isempty(last_display) || (etime(clock, last_display) > params.DISPLAY_PERIOD)
% %     last_display = clock;
% %     
% %     invalid = repmat(~data.valid, [1,1,3]);
% %     Z(~data.valid) = nan;
% %     S(invalid) = nan;
% %     A(invalid) = nan;
% %     
% %     I = data.true.im;
% %     
% %     L = visSH_color(L_expected, size(Z));
% % %     L = visSH_color(data.true.light, size(Z));
% %     L = L ./ max(L(:));
% %     
% %     A = exp(A);
% %     S = exp(S);
% %     A = A ./ max(A(:));
% %     S = S ./ max(S(:));
% %     
% %     state = struct('normal', N, 'height', Z, 'albedo', A, 'shading', S);
% %     
% %     [err, state] = getError(state, data.true);
% %     
% %     I(invalid) = 0;
% %     data.true.shading(invalid) = 0;
% %     L(isnan(L)) = 0;
% %     state.albedo(invalid) = 0;
% %     state.shading(invalid) = 0;
% %     
% %     Zv = visualizeDEM([state.height; data.true.height]);
% %     Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
% %     Nv = visualizeNormals_color(N);
% %     Ntv = visualizeNormals_color(data.true.normal);
% %     Ntv(isnan(data.true.normal)) = 0;
% %     Nv(isnan(data.true.normal)) = 0;
% %     V = min(1, [Zv, [Nv; Ntv], [state.shading; data.true.shading], [state.albedo; data.true.albedo], [L; I]]);
% %     
% %     figure(1);
% %     imagesc(V);
% %     axis image off;
% %     set(gca, 'PlotBoxAspectRatio', [1 1 1])
% %     set(gca, 'Position', [0 0 1 1]);
% %     set(gcf, 'Position', [1, 500, 1600, 1600/size(V,2)*size(V,1)])
% %     %       set(gcf, 'Position', get(gcf, 'Position') .* [1 1 0 0] + [0 0 size(V,2) size(V,1)]);
% %     
% %     
% %     global global_errors;
% %     global_errors{end+1} = err;
% %     height_errs = cellfun(@(x) x.height, global_errors);
% % %     affine_height_errs = cellfun(@(x) x.affine_height, global_errors);
% %     survey_errs = cellfun(@(x) x.survey, global_errors);
% %     albedo_errs = cellfun(@(x) x.albedo, global_errors);
% %     shading_errs = cellfun(@(x) x.shading, global_errors);
% %     grosse_errs = cellfun(@(x) x.grosse, global_errors);
% %     
% %     figure(2);
% %     plot(height_errs ./ max(eps,height_errs(1)), 'b-'); hold on;
% % %     plot(affine_height_errs ./ max(eps,affine_height_errs(1)), 'Color', [0.5, 0.5, 1]);
% %     plot(survey_errs ./ max(eps,survey_errs(1)), 'r-');
% %     plot(grosse_errs ./ max(eps,grosse_errs(1)), 'k-');
% %     plot(albedo_errs ./ max(eps,albedo_errs(1)), 'g-');
% %     plot(shading_errs ./ max(eps,shading_errs(1)), 'm-');
% %     axis square
% %     ylim = get(gca, 'YLim');
% %     set(gca, 'YLim', [.09, 1.1])
% %     set(gca, 'YScale', 'log');
% %     title(['b = Z(', num2str(height_errs(end)),'), r = I(', num2str(survey_errs(end)),'), k = grosse(', num2str(grosse_errs(end)),'), g = A, m = S, m = loss, ', num2str(sum(global_L_active)), ' L active'])
% %     hold off;
% %     
% %     drawnow;
% %     
% %   end
% % end
% % 
% % 
% % % 
% % % beta = params.ADAPTIVE_BETA{1};
% % % 
% % % if beta == 0
% % % 
% % %   loss = loss_A + loss_Z;
% % %   d_loss_Z = d_loss_Z + d_loss_ZA;
% % % 
% % % else
% % %   loss_A = loss_A * beta;
% % %   loss_Z = loss_Z * beta;
% % %   
% % %   d_loss_ZA = d_loss_ZA * beta;
% % %   d_loss_Z = d_loss_Z * beta;
% % %   
% % %   LL_A = -loss_A;
% % %   LL_Z = -loss_Z;
% % %   m = max(LL_A, LL_Z);
% % %   P_A = exp(LL_A - m);
% % %   P_Z = exp(LL_Z - m);
% % %   s = P_A + P_Z;
% % %   P_A = P_A ./ s;
% % %   P_Z = P_Z ./ s;
% % %   
% % %   loss = loss_A * P_A + loss_Z * P_Z;
% % %   
% % %   dL_A = P_A * ( loss - loss_A + 1 );
% % %   dL_Z = P_Z * ( loss - loss_Z + 1 );
% % %   
% % %   d_loss_Z = dL_Z * d_loss_Z + dL_A * d_loss_ZA;
% % %   
% % %   d_loss_Z = d_loss_Z / (0.5*beta);
% % %   loss = loss / (0.5*beta);
% % % 
% % % end
% % % 
% % % 








% function [loss, d_loss, L_expected, L_max, L_prob] =  lossfun_srfs_color(state, data, params)
% 
% state.height = upsample(state.Zdown, data.pind);
% 
% Z = state.height;
% 
% [N, dN_Z] = getNormals_conv(Z);
% N_vec = reshape(N, [], 3);
% 
% [loss_Z, d_loss_Z] = priorZ(Z, data, params);
% [loss_Z_simple, d_loss_Z_simple] = priorZ_simple(Z, data, params, N, dN_Z);
% 
% L_meta = data.L_meta;
% 
% Ls = data.true.light;
% 
% L_active_idx = find(data.L_active);
% 
% alpha = params.ALPHA_LIGHT{1};
% 
% d_loss_N = 0;
% d_losses_N = {};
% 
% % keyboard
% 
% if params.USE_GREEDY
% 
%   losses_A = (10^10) * ones(1, size(Ls,3));
%   LL_online = -inf(1, size(Ls,3));
%   loss_online = 0;
%   
%   
% %   if params.GREEDY_MOVE_HEAD
% %     global global_greedy_head
% %     if isempty(global_greedy_head)
% %       global_greedy_head = L_meta.head;
% %     end
% %     head = global_greedy_head;
% %   else
% %     head= L_meta.head;
% %   end
%   
%   head= L_meta.head;
%   queue = [head, -inf];
%   
%   closed_set = false(size(Ls,3),1);
%   for i_search = 1:params.GREEDY_N_SEARCH
%     
% %     queue
%     
%     [junk, min_i] = min(queue(:,2));
%     i = queue(min_i,1);
%     queue = queue([[1:(min_i-1)], [(min_i+1):size(queue,1)]],:);
%     closed_set(i) = true;
%     L = Ls(:,:,i);
%     
%     S = {};
%     dS_N = {};
%     for c = 1:3
%       [s, ds_n] = renderSH_helper(N_vec, L(:,c));
%       S{c} = reshape(s, size(Z));
%       dS_N{c} = reshape(ds_n, [size(Z), 3]);
%     end
%     S = cat(3, S{:});
%     
%     A = data.true.log_im - S;
%     
%     [loss_A, d_loss_A] = priorA(A, data, params);
%     losses_A(i) = loss_A;
%     
%     d_loss_Ni = 0;
%     for c = 1:3
%       d_loss_Ni = d_loss_Ni + repmat(-d_loss_A(:,:,c), [1,1,3]) .* dS_N{c};
%     end
%     d_losses_N{i} = d_loss_Ni;
%     
%     kids = find((L_meta.children(:,i) | L_meta.children(i,:)') & ~closed_set);
% %     kids = find(L_meta.children(:,i));
%     kids_depth = L_meta.depth(kids);
%     queue = [queue; [kids, loss_A + kids_depth*params.GREEDY_DEPTH_COST]];
%     
%   end
%   
%   L_active_idx = find(losses_A  <  10^10);
%   
%   LL = -alpha * losses_A;
%   LL = LL - max(LL);
%   P = exp(LL);
%   P = P ./ sum(P);
%   
%   loss_A = sum(losses_A .* P);
%   L_expected = sum(Ls .* repmat(reshape(P, [1,1, size(P,2)]), size(Ls,1), size(Ls,2)),3);
%   [junk, max_idx] = max(P);
%   L_max = Ls(:,:,max_idx);
%   L_prob = P;
%   
%   
%   d_loss_N = 0;
%   
%   dLs = P .* ( alpha*(loss_A - losses_A) + 1 );
%   for i = L_active_idx
%     d_loss_N = d_loss_N + dLs(i) * d_losses_N{i};
%   end
% 
% %   for i = L_active_idx
% %     d_loss_N = d_loss_N + P(i) * d_losses_N{i};
% %   end
% 
% else
%   
%   losses_A = (10^10) * ones(1, size(Ls,3));
%   LL_online = -inf(1, size(Ls,3));
%   loss_online = 0;
%   for i = L_active_idx
%     L = Ls(:,:,i);
%     
%     S = {};
%     dS_N = {};
%     for c = 1:3
%       [s, ds_n] = renderSH_helper(N_vec, L(:,c));
%       S{c} = reshape(s, size(Z));
%       dS_N{c} = reshape(ds_n, [size(Z), 3]);
%     end
%     S = cat(3, S{:});
%     
%     A = data.true.log_im - S;
%     
%     [loss_A, d_loss_A] = priorA(A, data, params);
%     losses_A(i) = loss_A;
%     
%     d_loss_Ni = 0;
%     for c = 1:3
%       d_loss_Ni = d_loss_Ni + repmat(-d_loss_A(:,:,c), [1,1,3]) .* dS_N{c};
%     end
%     
%     if params.USE_APPROX_LIGHT_GRADIENT
%       
%       LL_online(i) = -alpha*loss_A;
%       P_online = exp(LL_online - max(LL_online));
%       P_online = P_online ./ sum(P_online);
%       
%       loss_online = P_online(i).* loss_A + (1-P_online(i)) * loss_online;
%       
%       d_loss_N = P_online(i).* d_loss_Ni + (1-P_online(i)) * d_loss_N;
%     else
%       d_losses_N{i} = d_loss_Ni;
%     end
%     
%   end
%   
%   LL = -alpha * losses_A;
%   LL = LL - max(LL);
%   P = exp(LL);
%   P = P ./ sum(P);
%   
%   loss_A = sum(losses_A .* P);
%   L_expected = sum(Ls .* repmat(reshape(P, [1,1, size(P,2)]), size(Ls,1), size(Ls,2)),3);
%   L_prob = P;
%   
%   
%   if ~params.USE_APPROX_LIGHT_GRADIENT
%     dLs = P .* ( alpha*(loss_A - losses_A) + 1 );
%     for i = L_active_idx
%       d_loss_N = d_loss_N + dLs(i) * d_losses_N{i};
%     end
%   end
% 
% end
% 
% d_loss_ZA = ...
%   + conv3d( dN_Z.F1_1 .* d_loss_N(:,:,1) + dN_Z.F1_2 .* d_loss_N(:,:,2) + dN_Z.F1_3 .* d_loss_N(:,:,3), dN_Z.f1) ...
%   + conv3d( dN_Z.F2_1 .* d_loss_N(:,:,1) + dN_Z.F2_2 .* d_loss_N(:,:,2) + dN_Z.F2_3 .* d_loss_N(:,:,3), dN_Z.f2);
% 
% 
% loss = loss_A + loss_Z + loss_Z_simple;
% d_loss_Z = d_loss_Z + d_loss_ZA + d_loss_Z_simple;
% 
% d_loss.Zdown = downsample(d_loss_Z, size(data.pind,1), namedFilter(PYR_FILTER_NAME)*2);
% % d_loss.Zdown = downsample(d_loss_Z, size(data.pind,1), [1;4;6;4;1]/8);
% 
% 
% global global_loss_best
% if isempty(global_loss_best)
%   global_loss_best = inf;
% end
% 
% if (loss <= global_loss_best)
%   global_loss_best = loss;
%   
% %   global global_greedy_head
% %   [m, m_idx] = max(LL);
% %   if m > (log(8) + LL(global_greedy_head))
% %     global_greedy_head = m_idx;
% %     if params.GREEDY_MOVE_HEAD
% %       fprintf('New Head: %d\n', global_greedy_head);
% %     end
% %   end
% 
%   global global_P global_losses_A
%   global_P = P;
%   global_losses_A = losses_A;
%   
%   losses_A(losses_A == 10^10) = nan;
%   
%   
%   global global_light_fid;
%   fprintf(global_light_fid, '%e ', losses_A);
%   fprintf(global_light_fid, '\n');
% 
%   
%   if params.DO_DISPLAY
%     
%     global last_display
%     if isempty(last_display) || (etime(clock, last_display) > params.DISPLAY_PERIOD)
%       last_display = clock;
%       
%       invalid = repmat(~data.valid, [1,1,3]);
%       Z(~data.valid) = nan;
%       
%       I = data.true.im;
%       
%       L = visSH_color(L_expected, size(Z));
%       %     L = visSH_color(data.true.light, size(Z));
%       L = L ./ max(L(:));
% 
%       S = {};
%       for c = 1:3
%         S{c} = reshape(renderSH_helper(N_vec, L_expected(:,c)), size(Z));
%       end
%       S = cat(3,S{:});
% 
%       A = data.true.log_im - S;
%       
%       A = exp(A);
%       S = exp(S);
%       
%       A = max(0, min(1, A));
%       S = max(0, min(1, S));
%       
%       S(invalid) = nan;
%       A(invalid) = nan;
% 
%       
%       state = struct('normal', N, 'height', Z, 'albedo', A, 'shading', S);
%       
%       [err, state] = getError(state, data.true);
%       
%       I(invalid) = 0;
%       data.true.shading(invalid) = 0;
%       L(isnan(L)) = 0;
%       state.albedo(invalid) = 0;
%       state.shading(invalid) = 0;
%       
%       Zv = visualizeDEM([state.height; data.true.height]);
%       Zv(repmat(all(Zv == 1,3), [1,1,3])) = 0;
%       Nv = visualizeNormals_color(N);
%       Ntv = visualizeNormals_color(data.true.normal);
%       Ntv(isnan(data.true.normal)) = 0;
%       Nv(isnan(data.true.normal)) = 0;
%       V = min(1, [Zv, [Nv; Ntv], [state.shading; data.true.shading], [state.albedo; data.true.albedo], [L; I]]);
%       
%       figure(1);
%       imagesc(V);
%       axis image off;
%       set(gca, 'PlotBoxAspectRatio', [1 1 1])
%       set(gca, 'Position', [0 0 1 1]);
%       set(gcf, 'Position', [1, 500, 1600, 1600/size(V,2)*size(V,1)])
%       %       set(gcf, 'Position', get(gcf, 'Position') .* [1 1 0 0] + [0 0 size(V,2) size(V,1)]);
%       
%       
%       global global_errors;
%       global_errors{end+1} = err;
%       height_errs = cellfun(@(x) x.height, global_errors);
%       %     affine_height_errs = cellfun(@(x) x.affine_height, global_errors);
%       survey_errs = cellfun(@(x) x.survey, global_errors);
%       albedo_errs = cellfun(@(x) x.albedo, global_errors);
%       shading_errs = cellfun(@(x) x.shading, global_errors);
%       grosse_errs = cellfun(@(x) x.grosse, global_errors);
%       
%       figure(2);
%       plot(height_errs ./ max(eps,height_errs(1)), 'b-'); hold on;
%       %     plot(affine_height_errs ./ max(eps,affine_height_errs(1)), 'Color', [0.5, 0.5, 1]);
%       plot(survey_errs ./ max(eps,survey_errs(1)), 'r-');
%       plot(grosse_errs ./ max(eps,grosse_errs(1)), 'k-');
%       plot(albedo_errs ./ max(eps,albedo_errs(1)), 'g-');
%       plot(shading_errs ./ max(eps,shading_errs(1)), 'm-');
%       axis square
%       ylim = get(gca, 'YLim');
%       set(gca, 'YLim', [.09, 1.1])
%       set(gca, 'YScale', 'log');
%       title(['b = Z(', num2str(height_errs(end)),'), r = I(', num2str(survey_errs(end)),'), k = grosse(', num2str(grosse_errs(end)),'), g = A, m = S, m = loss'])
%       hold off;
%       
%       X = L_meta.X;
%       Y = L_meta.Y;
%       children = L_meta.children;
%       V = L_meta.vis;
%       
%       children = children(L_active_idx, L_active_idx);
%       X = X(L_active_idx);
%       Y = Y(L_active_idx);
%       V = V(L_active_idx);
%       
%       [i,j] = find(children);
%       
%       figure(4);
%       clf
% 
%       line(X([i,j])', Y([i,j])', 'color', 'k')
%       axis equal off
%       P_active = P(L_active_idx);
% %       sort(P_active, 'descend')
%       Q = log(P_active ./ max(P_active));
%       Q = max(0, min(1, (Q + 5)/5));
%       
%       hold on;
%       for vi = 1:length(V)
%         v = V{vi};
%         v = cat(3, padarray(v(:,:,1), [3,3], 1), padarray(v(:,:,2:3), [3,3], 1-Q(vi)));
%         imagesc(X(vi) - size(V{vi},1)/2, Y(vi) - size(V{vi},1)/2, v );
%       end
%       hold off;
%       
%       drawnow;
%       
%     end
%   end
%   
% end
