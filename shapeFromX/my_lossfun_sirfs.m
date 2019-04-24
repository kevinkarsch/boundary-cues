function [loss, d_loss, state_return] =  my_lossfun_sirfs(state, data, params)

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

[loss_Z, d_loss_Z, losses_Z] = myPriorZ(Z, data, params, N, dN_Z);

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
