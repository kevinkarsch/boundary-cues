function [results, state, data, params] = go_Solve(eval_string)

addpath(genpath('./IPPConvsToolbox'))
addpath(genpath('./matlabPyrTools/'));
addpath('./binaryops/');
addpath('./bmorph/');
addpath('./fminlbfgs_version2c');
addpath('./fftconv')

rand('twister',5489)
randn('state',0)

curdir = pwd;
if strcmp(curdir(1:10), '/Users/jon')
%   compile
else
  fprintf('not compiling\n');
  
  try
    maxNumCompThreads(1)
  catch ME
    fprintf('ERR: %s\n', ME.message);
  end
  
end

if nargin == 0
  eval_string = '';
end

params.EVAL_STRING = eval_string;
eval(params.EVAL_STRING);

CONSTANTS;

fprintf('%s\n', params.EVAL_STRING);
eval(params.EVAL_STRING);

BASE_FOLDER = 'data/MITC-intrinsic/data/';
FALSE_COLOR_FOLDER = 'data/MIT-colorized/data/';

load('prior.mat')
% load('prior_good.mat')

fprintf('params.MULT_OPTS.srfs.albedo = \n');
disp(params.MULT_OPTS.saifs.albedo)

fprintf('params.MULT_OPTS.srfs.height = \n');
disp(params.MULT_OPTS.saifs.height)

fprintf('params.MULT_OPTS.srfs.light = \n');
disp(params.MULT_OPTS.saifs.light)

if isfield(params, 'EVAL_NAMES')
  names = params.EVAL_NAMES;
else
%   names = {'apple'};
%   names = {'frog1'};
%   names = {'paper1'};
%   names = {'frog2'};
%   names = {'turtle'};
%   names = {'panther'};
%   names = {'raccoon'};
%   names = {'dinosaur'};
  names = {'Lpaper2'};
%   names = {'cup2'};
%   names = {'box'};
%   names = {'sun'};
end

for name_i = 1:length(names)
  name = names{name_i};
    
  if name(1) == 'C'
    FALSE_COLOR_LIGHT = 1;
  elseif name(1) == 'L'
    FALSE_COLOR_LIGHT = 0;
  else
    assert(1==0);
  end
  subname = name(2:end);
  
  
  load([BASE_FOLDER, subname, '/Z.mat']);
  clear lights
  data.true.height = depth;
  
  data.true.albedo = imread([BASE_FOLDER, subname, '/reflectance.png']);
  data.true.albedo = double(data.true.albedo) ./ double(intmax('uint16'));

  
  if FALSE_COLOR_LIGHT
  
    data.true.shading = imread([FALSE_COLOR_FOLDER, subname, '/shading_color.png']);
    data.true.shading = double(data.true.shading) ./ double(intmax('uint16'));
    
    data.true.im = imread([FALSE_COLOR_FOLDER, subname, '/diffuse.png']);
    data.true.im = double(data.true.im) ./ double(intmax('uint16'));
    
    load([FALSE_COLOR_FOLDER, subname, '/L.mat']);
    data.true.light = L;
    
    params.L_WHITEN_PARAMS = 'natural_color';
    
  else
    
    loaded = load([BASE_FOLDER, subname, '/shading_corrected_color.mat']);
    data.true.shading = loaded.shading_correct;
    
    data.true.im = imread([BASE_FOLDER, subname, '/diffuse.png']);
    data.true.im = double(data.true.im) ./ double(intmax('uint16'));
    
    load([BASE_FOLDER, subname, '/L.mat']);
    data.true.light = color_lights.diffuse;

    params.L_WHITEN_PARAMS = 'lab_color';
    
  end

%   if params.USE_BOTH_LIGHTS
%     params.L_WHITEN_PARAMS = 'both_color';
%   end
  
  data.true.mask = all(imread([BASE_FOLDER, subname, '/mask.png']) > 0,3);
  

  if params.FORCE_MONOCHROME_LIGHT
    data.true.light = repmat(mean(data.true.light,2), [1,3]);
  end

  
%   data.true.shading = min(1, data.true.im ./ data.true.albedo);
  
%   if params.FALSE_COLOR_LIGHT
% %     load(['data/MITC-intrinsic/data/panther/log_lights.mat']);
% %     L = [color_lights.light{1}(:,1), color_lights.light{2}(:,2), color_lights.light{3}(:,3)];
%     
% %     L = data.true.light;
% %     rand('twister',5489)
% %     randn('state',0)
% %     L = L + normrnd(0, 1, size(L))/10;
% 
%     load(['blob_lights/light_', num2str(params.FALSE_COLOR_LIGHT, '%02d'), '.mat'])
%     
%     % L = [color_lights.light{1}(:,1), color_lights.light{1}(:,2), color_lights.light{1}(:,3)];
%     log_shading = renderSH(depth, L_color);
%     log_albedo = log(data.true.albedo);
%     
%     data.true.shading = exp(log_shading);
%     data.true.log_im = log_shading + log_albedo;
%     data.true.im = exp(data.true.log_im);
%     
%     if params.FORCE_MONOCHROME_LIGHT
%       data.true.light = L_gray;
%     else
%       data.true.light = L_color;
%     end
%     
%     
%   end
    
  
  load([BASE_FOLDER, subname, '/crop_idx.mat'])
    
  for field = {'height', 'albedo', 'im', 'shading', 'mask'}
    field = field{1};
    data.true.(field) = data.true.(field)(crop_idx1, crop_idx2,:);
  end
  
  
  if params.RESIZE_INPUT ~= 1
  
    sz_before = size(data.true.height);
    
    Z = inpaintZ(data.true.height, 0, 1);
    
    S = data.true.shading;
    S(S==0) = nan;
    S = inpaintZ(S, 1, 0);
    
    A = data.true.albedo;
    A(A==0) = nan;
    A = inpaintZ(A, 1, 0);
    
    I = data.true.im;
    I(I==0) = nan;
    I = inpaintZ(I, 1, 0);
    
    M = double(data.true.mask);
    
    I = max(0, min(1, imresize(I, params.RESIZE_INPUT)));
    S = max(0, min(1, imresize(S, params.RESIZE_INPUT)));
    A = max(0, min(1, imresize(A, params.RESIZE_INPUT)));
    Z = imresize(Z, params.RESIZE_INPUT)*params.RESIZE_INPUT;
    %     plot(Z2(end/2,:)); hold on; plot(Z(end/2,:)); hold on;
    M = imresize(M, params.RESIZE_INPUT) > 0.5;
    
    Z(~M) = nan;
    I(repmat(~M, [1,1,3])) = 0;
    S(repmat(~M, [1,1,3])) = 0;
    A(repmat(~M, [1,1,3])) = 0;
    
    data.true.shading = S;
    data.true.albedo = A;
    data.true.im = I;
    data.true.mask = M;
    data.true.height = Z;
    
  end
  
  
  im = data.true.im;
  
%   invalid = any(im == 0,3);
%   im(im == 0) = nan;
%   im = inpaintZ(im, 1, 0);
%   im = max(0, min(1, blur3(im)));
%   im(repmat(invalid, [1,1,3])) = 0;
%   data.true.im = im;
  
  valid = all(im > 0,3);
  im(repmat(~valid, [1,1,3])) = nan;  
  log_im = log(im);
%   log_im = inpaintZ(log_im, 1, 0);
%   im = exp(log_im);
  
  data.true.im = im;
  data.true.log_im = log_im;
    
  data.valid = valid;
%   data.valid_pyr = buildGpyr_simple(double(data.valid), params.APYR_DEPTH{1}, [1;4;6;4;1]/16, 'repeat');
%   data.valid_pyr = cellfun(@(x) x > 0.5, data.valid_pyr, 'UniformOutput', false);

%   data.median_filter_mats = {};
%   for s = 1:length(data.valid_pyr)
%     data.median_filter_mats{s} = medianFilterMat_mask(~data.valid_pyr{s}, params.Z_MEDIAN_HALFWIDTH);
%     [i, j] = find((data.median_filter_mats{s}~=0)');
%     data.median_filter_pairs{s} = [i(1:2:end), i(2:2:end)];
% %     data.median_filter_mats_weighted{s} = medianFilterMat_mask_weighted(~data.valid_pyr{s}, params.Z_MEDIAN_HALFWIDTH);
%   end

  if params.Z_MEDIAN_HALFWIDTH == 0
    data.Z_median_filter_mat = medianFilterMat_mask_4conn(~data.valid);
  else
    data.Z_median_filter_mat = medianFilterMat_mask(~data.valid, params.Z_MEDIAN_HALFWIDTH);
  end
  
  if params.A_MEDIAN_HALFWIDTH == 0
    data.A_median_filter_mat = medianFilterMat_mask_4conn(~data.valid);
  else
    data.A_median_filter_mat = medianFilterMat_mask(~data.valid, params.A_MEDIAN_HALFWIDTH);
  end
  
  data.Z_median_filter_mat_T = data.Z_median_filter_mat';
  data.A_median_filter_mat_T = data.A_median_filter_mat';
  
  data.border = [];
  [data.border.mask,data.border.normal] = getBorderNormals(data.true.mask);
  data.true.normal = getNormals_conv(data.true.height);
  data.true.curvature = getK(data.true.height);
  
  data.prior = prior;

  for v = params.GLOBAL_VARS
    eval(['global ', v{1}, ';']);
    eval([v{1}, ' = [];']);
  end
  
  load Ldata
  
  if params.DISABLE_LIGHT_PRIOR
%     data.L_whiten_params = Ldata.whiten.(params.L_WHITEN_PARAMS);
    data.L_whiten_params.mean = zeros(1,27);
    data.L_whiten_params.map = eye(27,27);
    data.L_whiten_params.inverse = eye(27,27);
    data.L_gaussian.mu = zeros(1,27);
  else  
    data.L_whiten_params = Ldata.whiten.(params.L_WHITEN_PARAMS);
    data.L_gaussian = Ldata.gaussian.(params.L_WHITEN_PARAMS);
  end
    
  
%   data.L_GSM = Ldata.GSM.(params.L_WHITEN_PARAMS);

  if params.MAKE_VIDEO
    mkdir video_frames;
    data.video_folder = ['video_frames/', num2str(params.EXPERIMENT), '_', num2str(params.SOLVE_LIGHT), '_', name, '/'];
    mkdir(data.video_folder)
    system(['rm -rf ', data.video_folder, '*']);
%     data.video_log = fopen([data.video_folder, 'log.txt'], 'w');
    global GLOBAL_VIDEO_COUNT;
    GLOBAL_VIDEO_COUNT = 1;
  end
  
  start_time = clock;
    
%   params.LOSSFUN = 'lossfun_srfs_color_light2';
%   [state] = do_Solve_Light2(data, params);

    
%   params.LOSSFUN = 'lossfun_srfs_color_light';
%   [state] = do_Solve_Light(data, params);
  
%   params.LOSSFUN = 'lossfun_srfs_color_light_pyr';
%   
%   if params.DO_CHEAT
%     Z_true = inpaintZ(data.true.height, 0, 1);
%     L_true = data.true.light;
%     [state] = do_Solve_Light_Pyr(data, params, Z_true, L_true);
%   else
%     
%     if params.SOLVE_LIGHT == 1
%       params.LOSSFUN = 'lossfun_srfs_color_light_pyr_multi';  
%       [state] = do_Solve_Light_Pyr_Multi(data, params);
% %       params.LOSSFUN = 'lossfun_srfs_color_light_pyr';  
% %       [state] = do_Solve_Light_Pyr(data, params);
%     else
%       [state] = do_Solve_Light_Pyr(data, params);
%     end
%   end
  

%   if params.USE_INCEPTION
%     params.LOSSFUN = 'lossfun_saifs_inception';
%   else
  params.LOSSFUN = 'lossfun_sirfs';
%   end

  
  if params.DO_CHEAT
    Z_true = inpaintZ(data.true.height, 0, 1);
    L_true = data.true.light;
    [state] = do_Solve(data, params, Z_true, L_true);
  else
    
%     if params.USE_INCEPTION
%       [state] = do_Solve_Inception(data, params);
%     else
    [state] = do_Solve(data, params);
%     end
    
  end

  if params.MAKE_VIDEO
%     fclose(data.video_log)
  end
  
  state.solve_time = etime(clock, start_time);
  
%   [S, junk, N] = renderSH(state.height, state.light);
%   state.normal = N;
%   A = data.true.log_im - S;
% 
%   state.shading = exp(S);
%   state.albedo = exp(A);
  
  losses = [];
  
  if params.DO_DISPLAY;
    Ztrue = data.true.height;
    Ztrue = inpaintZ(Ztrue, 0, 1);
%     [Strue] = renderSH(Ztrue, data.true.light);
%     Atrue = data.true.log_im - Strue;
    Atrue = data.true.albedo;
    Atrue(Atrue == 0) = nan;
    Atrue = log(Atrue);
%     Atrue = inpaintZ(Atrue, 1, 0);
    

    [junk, junk, losses.gt.height] = priorZ(Ztrue, data, params);
    [junk, junk, losses.gt.albedo] = priorA(Atrue, data, params);
    [junk, junk, losses.gt.light] = priorL(data.true.light, data, params);
    
    [junk, junk, losses.us.height] = priorZ(state.height, data, params);
    [junk, junk, losses.us.albedo] = priorA(log(state.albedo), data, params);
    [junk, junk, losses.us.light] = priorL(state.light, data, params);
    
    dispStruct('losses.gt.height', losses.gt.height)
    dispStruct('losses.gt.albedo', losses.gt.albedo)
    dispStruct('losses.gt.light', losses.gt.light)
    fprintf('Total gt loss: %f\n', sum(struct2vector(losses.gt)))
    fprintf('\n');
    dispStruct('losses.us.height', losses.us.height)
    dispStruct('losses.us.albedo', losses.us.albedo)
    dispStruct('losses.us.light', losses.us.light)
    fprintf('Total us loss: %f\n', sum(struct2vector(losses.us)))
    
%     [junk, template] = struct2vector(losses.us);
% 
%     pow = 1/8;
% %     pow = 1/32;
%     
%     Ms_new = [];
%     for i1 = 1:length(template.fields)
%       f1 = template.fields{i1};
%       Ms_new.(f1) = [];
%       for i2 = 1:length(template.shapes{i1}.fields)
%         f2 = template.shapes{i1}.fields{i2};
%         if iscell(params.MULT_OPTS.saifs.(f1).(f2))
%           Ms_new.(f1).(f2) = {params.MULT_OPTS.saifs.(f1).(f2){1} * ((losses.us.(f1).(f2) ./ max(eps,losses.gt.(f1).(f2))) .^ (pow))};
%         else
%           for i3 = 1:length(template.shapes{i1}.shapes{i2}.fields)
%             f3 = template.shapes{i1}.shapes{i2}.fields{i3};
%             Ms_new.(f1).(f2).(f3) = {params.MULT_OPTS.saifs.(f1).(f2).(f3){1} * ((losses.us.(f1).(f2).(f3) ./ max(eps,losses.gt.(f1).(f2).(f3))) .^ (pow))};
%           end
%         end
%       end
%     end
% %     
%     [Ms_vec, template] = struct2vector(Ms_new);
%     Ms_vec = Ms_vec ./ sum(Ms_vec);
%     Ms_norm = vector2struct(Ms_vec, template);
%     
%     dispStruct('params.MULT_OPTS.saifs.height', Ms_norm.height);
%     dispStruct('params.MULT_OPTS.saifs.albedo', Ms_norm.albedo);
%     dispStruct('params.MULT_OPTS.saifs.light', Ms_norm.light);
    
%     pow = 0.125;
%     
%     params.MULT_OPTS.saifs.height.smooth{1} = params.MULT_OPTS.saifs.height.smooth{1} * (losses.us.height.smooth ./ losses.gt.height.smooth)^pow;
%     params.MULT_OPTS.saifs.height.contour{1} = params.MULT_OPTS.saifs.height.contour{1} * (losses.us.height.contour ./ losses.gt.height.contour)^pow;
%     params.MULT_OPTS.saifs.height.slant{1} = params.MULT_OPTS.saifs.height.slant{1} * (losses.us.height.slant ./ losses.gt.height.slant)^pow;
% 
%     params.MULT_OPTS.saifs.albedo.smooth{1} = params.MULT_OPTS.saifs.albedo.smooth{1} * (losses.us.albedo.smooth ./ losses.gt.albedo.smooth)^pow;
%     params.MULT_OPTS.saifs.albedo.entropy{1} = params.MULT_OPTS.saifs.albedo.entropy{1} * (losses.us.albedo.entropy ./ losses.gt.albedo.entropy)^pow;
%     params.MULT_OPTS.saifs.albedo.hist{1} = params.MULT_OPTS.saifs.albedo.hist{1} * (losses.us.albedo.hist ./ losses.gt.albedo.hist)^pow;
%     
%     
%     total = params.MULT_OPTS.saifs.height.smooth{1} + params.MULT_OPTS.saifs.height.contour{1} + params.MULT_OPTS.saifs.height.slant{1} + params.MULT_OPTS.saifs.albedo.smooth{1} + params.MULT_OPTS.saifs.albedo.entropy{1} + params.MULT_OPTS.saifs.albedo.hist{1};
%     
%     params.MULT_OPTS.saifs.height.smooth{1} = params.MULT_OPTS.saifs.height.smooth{1} / total;
%     params.MULT_OPTS.saifs.height.contour{1} = params.MULT_OPTS.saifs.height.contour{1} / total;
%     params.MULT_OPTS.saifs.height.slant{1} = params.MULT_OPTS.saifs.height.slant{1} / total;
%     params.MULT_OPTS.saifs.albedo.smooth{1} = params.MULT_OPTS.saifs.albedo.smooth{1} / total;
%     params.MULT_OPTS.saifs.albedo.entropy{1} = params.MULT_OPTS.saifs.albedo.entropy{1} / total;
%     params.MULT_OPTS.saifs.albedo.hist{1} = params.MULT_OPTS.saifs.albedo.hist{1} / total;
%     
%     fprintf('If you liked that, you''ll love this\n');
%     fprintf('params.MULT_OPTS.saifs.height.smooth = { %f };\n', params.MULT_OPTS.saifs.height.smooth{1});
%     fprintf('params.MULT_OPTS.saifs.height.contour = { %f };\n', params.MULT_OPTS.saifs.height.contour{1});
%     fprintf('params.MULT_OPTS.saifs.height.slant = { %f };\n', params.MULT_OPTS.saifs.height.slant{1});
%     fprintf('\n');
%     fprintf('params.MULT_OPTS.saifs.albedo.smooth = { %f };\n', params.MULT_OPTS.saifs.albedo.smooth{1});
%     fprintf('params.MULT_OPTS.saifs.albedo.entropy = { %f };\n', params.MULT_OPTS.saifs.albedo.entropy{1});
%     fprintf('params.MULT_OPTS.saifs.albedo.hist = { %f };\n', params.MULT_OPTS.saifs.albedo.hist{1});
    
  end
  
  [err] = getError(state, data.true);
    
  state_init = [];
  state_init.height = zeros(size(state.height));
  [S, junk, N] = renderSH(state_init.height, data.true.light);
  A = data.true.log_im - S;       
  state_init.shading = exp(S);
  state_init.albedo = exp(A);
  state_init.normal = N;

  state_init.light = reshape(data.L_gaussian.mu, [9,3]);
  err_init = getError(state_init, data.true);
  
  [v, template] = struct2vector(err);
  [vi, template] = struct2vector(err_init);
  err_percent = vector2struct(v ./ max(eps,vi) * 100, template)

  results{name_i}.err = err;
  results{name_i}.err_init = err_init;
  results{name_i}.err_percent = err_percent;
  results{name_i}.mult = params.MULT_OPTS;
  results{name_i}.losses = losses;
  results{name_i}.solve_time = state.solve_time;
  
  state_pad = state;
  state_pad.height = nan(crop_init_size);
  state_pad.albedo = nan([crop_init_size,3]);
  state_pad.albedo_exp = nan([crop_init_size,3]);
  state_pad.albedo_max = nan([crop_init_size,3]);
  state_pad.shading = nan([crop_init_size,3]);
  
  if params.RESIZE_INPUT ~= 1
  
    state.shading = max(0, min(1, imresize(state.shading, sz_before)));
    state.albedo = max(0, min(1, imresize(inpaintZ(state.albedo, 1, 0), sz_before)));
    state.height = imresize(inpaintZ(state.height, 0, 1), sz_before)/params.RESIZE_INPUT;
    data.valid = imresize(double(data.valid), sz_before)>0;

  end
  
  
  state_pad.height(crop_idx1, crop_idx2) = state.height;
  state_pad.albedo(crop_idx1, crop_idx2,:) = state.albedo;
  state_pad.shading(crop_idx1, crop_idx2,:) = state.shading;
  
  invalid = true(crop_init_size);
  invalid(crop_idx1, crop_idx2) = ~data.valid;
  state_pad.height(invalid) = nan;
  state_pad.albedo(repmat(invalid,[1,1,3])) = nan;
  state_pad.shading(repmat(invalid,[1,1,3])) = nan;
  
  state = state_pad;
    
  if ~isempty(params.DUMP_OUTPUT)
    system(['mkdir ', params.DUMP_OUTPUT]);
    save([params.DUMP_OUTPUT, '/', name, '.mat'], 'state');
  end
  
end

if ~isempty(params.OUTPUT_FILENAME)

  fprintf('Saving results and params to %s... ', params.OUTPUT_FILENAME);
  save(params.OUTPUT_FILENAME, 'results', 'params');
  fprintf('done.\n');

end

save last.mat state

% data = rmfield(data, 'prior');
% save output state state_im results data

