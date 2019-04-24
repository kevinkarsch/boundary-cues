params.OUTPUT_FILENAME = [];

params.SOLVE_LIGHT = 0;
params.SOLVE_SHAPE = 1;

params.SFS = 0;

params.GAMUT_POWER = 2;

params.APYR_ROOT = {4};
params.ZPYR_ROOT = {4};

params.APYR_DEPTH = 5;
params.ZPYR_DEPTH = 5;
    
params.ONE_ALBEDO_MODEL = 0;
params.ONE_SHAPE_MODEL = 0;

params.DISABLE_ZPYR = 0;

params.VIDEO_LOSS_EPS = 0.005;

params.FAST_LIGHT_SOLVE = 1;

% params.ALBEDO_RANGE_LOW = -5;
% params.ALBEDO_RANGE_HIGH = 0;
% params.ALBEDO_RANGE_POWER = 2;

params.ENTROPY_METHOD = 'renyi';
% params.ENTROPY_METHOD = 'parzen';

params.YC_MODEL = 'model1';

params.ADAPTIVE_BETA = {0};

% params.CURVATURE_MODE = 'K';

params.MULT_OPTS.sfs = { .01 };

params.FALSE_COLOR_LIGHT = 0;
params.FORCE_MONOCHROME_LIGHT = 0;

params.WHITEN_ALBEDO_ENTROPY = 1;
params.ALBEDO_ENTROPY_MINSCALE = 1;
params.ENTROPY_PROJECTIONS_ORDER = 4;

params.USE_ORACLE = 0;

params.USE_APPROX_RENYI = 0;

params.PIXEL_PAD = 1;

params.USE_KZUP = 1;

params.RESTART_LIGHTS = 0;
params.PA_EPSILON = {0.001};
params.L_D_MIN = {0.01};

params.ZPYR_MULT = {2};
params.L_MULT = {1};

params.MA_MODEL = 'MA';

params.PRIOR_PDF = 'GSM';
% params.PRIOR_PDF = 'myGaussian';

params.INIT_SHAPE_FROM_CONTOUR = 0;

params.ALBEDO_BIN_LOW =  [-7, -7, -7];
params.ALBEDO_BIN_HIGH = [ 4,  4,  4];

% params.FREE_LIGHT_SHIFT = 0;

params.GAMUT_POWER = 1;

params.L_WHITEN_PARAMS = 'natural_color'; %'natural_gray', 'lab_gray', 'lab_color'
    
params.RESIZE_INPUT = 1;

params.LIGHT_ALPHA = {0};

params.N_LIGHTS = {1};
params.DO_LIGHT_POSTERIOR = {1};

params.Z_MEDIAN_HALFWIDTH = 2;
params.A_MEDIAN_HALFWIDTH = 2;

params.A_SMOOTH_EPSILON = {0};
params.Z_SMOOTH_EPSILON = {0};
% params.A_SMOOTH_EPSILON = {10^-3};
% params.Z_SMOOTH_EPSILON = {10^-4};


params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };%{ 12.3377 };
params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

params.MULT_OPTS.saifs.light.lab_color.GSM =       { 0 };
params.MULT_OPTS.saifs.light.natural_color.GSM =   { 0 };

params.MULT_OPTS.saifs.light.lab_color.MOG =       { 0 };%{ 12.3377 };
params.MULT_OPTS.saifs.light.natural_color.MOG =   { 0 };

params.ALBEDO_SMOOTH_MODEL = 0;


params.WHITEN_LIGHT = 1;


params.ALBEDO_MODEL = 'ours';

params.MULT_OPTS.saifs.height.contour =         { 0.7071 };
params.MULT_OPTS.saifs.height.slant =   { 1.0905 };
params.MULT_OPTS.saifs.height.smooth =  { 3.668 };

params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
params.MULT_OPTS.saifs.albedo.hist =    { 9.5137 };
params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };

params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 11.3137 };
params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.181 };

params.DISABLE_LIGHT_PRIOR = 0;



if isfield(params, 'EXPERIMENT')
    
  if params.EXPERIMENT == 101
    
    
  elseif params.EXPERIMENT == 102 % No Entropy
    
    params.MULT_OPTS.saifs.albedo.entropy = { 0 };
    
  elseif params.EXPERIMENT == 103 % No Smoothness
    
    params.MULT_OPTS.saifs.albedo.smooth =  { 0 };
    
  elseif params.EXPERIMENT == 104 % No histogram
    
    params.MULT_OPTS.saifs.albedo.hist =    { 0 };

  elseif params.EXPERIMENT == 105 % No lights
    
    params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
    params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };
    params.DISABLE_LIGHT_PRIOR = 1;
    
  elseif params.EXPERIMENT == 106 % Monochrome lights
    
    params.FORCE_MONOCHROME_LIGHT = 1;
    
  elseif params.EXPERIMENT == 150 % Do nothing
    
    params.MULT_OPTS.saifs.height.contour =         { 0 };
    params.MULT_OPTS.saifs.height.slant =   { 0 };
    params.MULT_OPTS.saifs.height.smooth =  { 0 };
    
    params.MULT_OPTS.saifs.albedo.entropy =         { 0 };
    params.MULT_OPTS.saifs.albedo.hist =    { 0 };
    params.MULT_OPTS.saifs.albedo.smooth =  { 0 };
    
    params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
    params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

    
  elseif params.EXPERIMENT == 151 % Shape from Contour
    
    params.MULT_OPTS.saifs.height.contour =         { 3.3636 };
    params.MULT_OPTS.saifs.height.slant =   { 0.2973 };
    params.MULT_OPTS.saifs.height.smooth =  { 1024 };

    params.MULT_OPTS.saifs.albedo.entropy =         { 0 };
    params.MULT_OPTS.saifs.albedo.hist =    { 0 };
    params.MULT_OPTS.saifs.albedo.smooth =  { 0 };
    
    params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
    params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

    
    
  elseif params.EXPERIMENT == 121 % RGB
    
    params.ALBEDO_MODEL = 'rgb';
    params.MULT_OPTS.saifs.albedo.entropy =         { 0.5946 };
    params.MULT_OPTS.saifs.albedo.smooth =  { 4.362 };
    
  elseif params.EXPERIMENT == 122 % YUV
    
    params.ALBEDO_MODEL = 'yuv';
    params.MULT_OPTS.saifs.albedo.entropy =         { 1 };
    params.MULT_OPTS.saifs.albedo.smooth =  { 8.7241 };
        
    
  end
  
end


params.DUMP_OUTPUT = '';
params.DUMP_FINAL = false;


params.LMSE_WINDOW_SIZE = 20;

params.DEBUG_GRADIENT = 0;
params.USE_NUMERICAL = 0;
params.DO_DISPLAY = 0;
params.DISPLAY_NFUN = 1;
params.DISPLAY_PERIOD = 30;

params.LBFGS_NCORR = 100;
params.N_ITERS_OPTIMIZE = 1000;
params.N_SUB_ITERS = 1;
params.F_STREAK = 5;
params.F_PERCENT = .005;
params.PROG_TOL = 10^-12;
params.OPT_TOL = 10^-5;

% params.F_PERCENT = 0;
% params.PROG_TOL = 0;
% params.OPT_TOL = 0;


% params.N_ITERS_CONJUGATE = 50;
% params.PRE_SOLVE = 1;


params.MAX_N_EVAL = inf;
params.N_EVAL_SKIP = 0;

params.DO_CHEAT = 0;

params.GLOBAL_VARS = {'Z_last_global', 'Ls_white_global', 'global_Ls_dependency', 'global_Ls_weights', 'global_in_line_search', 'global_L_active_idx', 'global_greedy_head', 'global_loss_best', 'global_P', 'global_light_fid', 'last_display', 'num_display', 'display_figure', 'global_Z_last', 'global_L_best', 'global_L', 'global_state', 'global_errors', 'num_display', 'last_Zprior', 'global_losses', 'global_losses_history', 'GLOBAL_VIDEO_COUNT', 'global_loss_best_video'};



params.BLUR_STEP = 1;
params.MAX_PYR_DEPTH = inf;

params.MAKE_VIDEO = 0;

params.USE_SELF_OCCLUSION = 1; %Use self occlusion boundary information
