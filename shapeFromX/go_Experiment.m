function go_ExperimentFace(EXPERIMENT_NUM)

system('rm -rf *.err');
system('rm -rf *.log');

curdir = pwd;
if strcmp(curdir(1:10), '/Users/jon')
  fprintf('not compiling\n');
else
%   compile
end

params = [];
params.DO_DISPLAY = 0;
params.MAX_N_EVAL = 1;

OBJECTIVES = {'light', 'normal', 'grosse', 'shading', 'albedo'};
% OBJECTIVES = {'curvature', 'normal', 'grosse', 'shading', 'albedo'};

N_ROUNDS = 1;


if EXPERIMENT_NUM == 101
  
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  params.MULT_OPTS.saifs.height.smooth =  { 3.668 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 12.3377 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 3.0844 };

  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.height.smooth',
    'params.MULT_OPTS.saifs.height.slant',
    'params.MULT_OPTS.saifs.height.contour',
    'params.MULT_OPTS.saifs.albedo.hist',
    'params.MULT_OPTS.saifs.albedo.entropy',
    'params.MULT_OPTS.saifs.albedo.smooth',
    'params.MULT_OPTS.saifs.light.lab_color.gaussian'
    'params.MULT_OPTS.saifs.light.natural_color.gaussian',
    };
  
  SEARCH_VALS = [0, [2.^[-2:.125:4]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
elseif EXPERIMENT_NUM == 102
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.125 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 0 };
  params.MULT_OPTS.saifs.albedo.hist =            { 10.3747 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 4.362 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 6.7272 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.8284 };

  N_ROUNDS = 2;
  params.SOLVE_LIGHT = 1;

  params.ALBEDO_MODEL = 'rgb';
  
  search_strings = {
%     'params.MULT_OPTS.saifs.albedo.smooth',
    'params.MULT_OPTS.saifs.albedo.entropy',
    'params.MULT_OPTS.saifs.albedo.hist',
    'params.MULT_OPTS.saifs.height.smooth',
    'params.MULT_OPTS.saifs.height.slant',
    'params.MULT_OPTS.saifs.height.contour',
    'params.MULT_OPTS.saifs.light.lab_color.gaussian'
    'params.MULT_OPTS.saifs.light.natural_color.gaussian',
    };
  
  SEARCH_VALS = [0, [2.^[-2:.125:4]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  

elseif EXPERIMENT_NUM == 103
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.125 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 0 };
  params.MULT_OPTS.saifs.albedo.hist =            { 10.3747 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 8.7241 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 6.7272 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.8284 };

  N_ROUNDS = 2;
  params.SOLVE_LIGHT = 1;

  params.ALBEDO_MODEL = 'yuv';
  
  search_strings = {
%     'params.MULT_OPTS.saifs.albedo.smooth',
    'params.MULT_OPTS.saifs.albedo.entropy',
    'params.MULT_OPTS.saifs.albedo.hist',
    'params.MULT_OPTS.saifs.height.smooth',
    'params.MULT_OPTS.saifs.height.slant',
    'params.MULT_OPTS.saifs.height.contour',
    'params.MULT_OPTS.saifs.light.lab_color.gaussian'
    'params.MULT_OPTS.saifs.light.natural_color.gaussian',
    };
  
  SEARCH_VALS = [0, [2.^[-2:.125:4]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
elseif EXPERIMENT_NUM == 106
  
  params.MULT_OPTS.saifs.height.contour =         { 3.3636 };
  params.MULT_OPTS.saifs.height.smooth =  { 1024 };
  params.MULT_OPTS.saifs.height.slant =   { 0.2973 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.125 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 0 };
  params.MULT_OPTS.saifs.albedo.hist =            { 0 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 0 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

  N_ROUNDS = 2;
  params.SOLVE_LIGHT = 1;

  OBJECTIVES = {'normal'};
  
  params.ALBEDO_MODEL = 'ours';
  
  search_strings = {
    'params.MULT_OPTS.saifs.height.contour',
    'params.MULT_OPTS.saifs.height.smooth',
    'params.MULT_OPTS.saifs.height.slant',
    };
  
  SEARCH_VALS = [0, [2.^[-4:.25:12]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
  
elseif EXPERIMENT_NUM >= 201 && EXPERIMENT_NUM <= 210
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy = { 0 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 6.1688 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 6.7272 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.8284 };

  params.ALBEDO_HIST_MODEL = EXPERIMENT_NUM - 200;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.albedo.hist',
    };
  
  SEARCH_VALS = [2.^[1.5:.125:4.5]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
elseif EXPERIMENT_NUM >= 301 && EXPERIMENT_NUM <= 310
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy = { 0 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 6.7272 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.8284 };

  params.MULT_OPTS.saifs.albedo.entropy_sigma = { 2.^((EXPERIMENT_NUM - 323)/4) };
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.albedo.entropy',
    };
  
  SEARCH_VALS = [[2.^[-2:.125:0.5]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
elseif EXPERIMENT_NUM == 401
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 6.7272 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.8284 };

  params.ALBEDO_SMOOTH_MODEL = 0;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.albedo.smooth',
    };
  
  SEARCH_VALS = [0, [2.^[2:.125:5]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
elseif EXPERIMENT_NUM == 402
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };
  
  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 6.7272 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.8284 };

  params.ALBEDO_SMOOTH_MODEL = 1;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.albedo.smooth',
    };
  
  SEARCH_VALS = [0, [2.^[2:.125:5]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];

elseif EXPERIMENT_NUM == 501
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };

  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

  params.MULT_OPTS.saifs.light.lab_color.MOG = {0};
  params.MULT_OPTS.saifs.light.natural_color.MOG = {0};
  
  params.ALBEDO_SMOOTH_MODEL = 1;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.light.lab_color.gaussian',
    };
  
  SEARCH_VALS = [0, [2.^[1:.125:4]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'})];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
  
elseif EXPERIMENT_NUM == 502
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };

  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

  params.MULT_OPTS.saifs.light.lab_color.MOG = {0};
  params.MULT_OPTS.saifs.light.natural_color.MOG = {0};
  
  params.ALBEDO_SMOOTH_MODEL = 1;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.light.lab_color.MOG',
    };
  
  SEARCH_VALS = [0, [2.^[0:.125:5]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'})];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
elseif EXPERIMENT_NUM == 503
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };

  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

  params.MULT_OPTS.saifs.light.lab_color.MOG = {0};
  params.MULT_OPTS.saifs.light.natural_color.MOG = {0};
  
  params.ALBEDO_SMOOTH_MODEL = 1;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.light.natural_color.gaussian',
    };
  
  SEARCH_VALS = [0, [2.^[1:.125:4]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'C'})];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
  
elseif EXPERIMENT_NUM == 504
  
  params.MULT_OPTS.saifs.height.contour = { 0.7071 };
  params.MULT_OPTS.saifs.height.smooth =  { 2.5937 };
  params.MULT_OPTS.saifs.height.slant =   { 1 };
  
  params.MULT_OPTS.saifs.albedo.entropy_sigma =   { 0.0625 };
  params.MULT_OPTS.saifs.albedo.entropy =         { 1.1892 };
  params.MULT_OPTS.saifs.albedo.smooth =  { 6.1688 };
  params.MULT_OPTS.saifs.albedo.hist =    { 10.3747 };

  params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
  params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };

  params.MULT_OPTS.saifs.light.lab_color.MOG = {0};
  params.MULT_OPTS.saifs.light.natural_color.MOG = {0};
  
  params.ALBEDO_SMOOTH_MODEL = 1;
  
  N_ROUNDS = 1;
  params.SOLVE_LIGHT = 1;
  
  search_strings = {
    'params.MULT_OPTS.saifs.light.natural_color.MOG',
    };
  
  SEARCH_VALS = [0, [2.^[-2:.125:0]]];
  params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'C'})];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
  
  
  
  
else
  assert(1==0)
end
  
  
NAME = ['exp', num2str(EXPERIMENT_NUM), '.'];


search_strings_expand = {};
for i = 1:length(search_strings)
  x = eval(search_strings{i});
  for j = length(x):-1:1
    s = [search_strings{i}, '{', num2str(j), '}'];
    search_strings_expand{length(search_strings_expand)+1} = s;
  end
end
search_strings = search_strings_expand;

search_strings_sorted = {};
for j = length(x):-1:1
  take = cellfun(@(x) ~isempty(strfind(x, num2str(j))), search_strings, 'UniformOutput', false);
  take = cat(2,take{:});
  search_strings_sorted = [search_strings_sorted, search_strings(take)];
end
search_strings = search_strings_sorted;

cellfun(@(x) fprintf([x, '\n']), search_strings)

[boosted_params, history_err, history_param] = coordExperiment(NAME, params, search_strings, SEARCH_VALS, N_ROUNDS, OBJECTIVES);
































% function go_ExperimentFace(EXPERIMENT_NUM)
% 
% system('rm -rf *.err');
% system('rm -rf *.log');
% 
% curdir = pwd;
% if strcmp(curdir(1:10), '/Users/jon')
%   fprintf('not compiling\n');
% else
%   compile
% end
% 
% params = [];
% params.DO_DISPLAY = 0;
% params.MAX_N_EVAL = 1;
% 
% OBJECTIVES = {'light', 'normal', 'grosse', 'shading', 'albedo'};
% % OBJECTIVES = {'curvature', 'normal', 'grosse', 'shading', 'albedo'};
% 
% N_ROUNDS = 1;
% 
% 
% params.MULT_OPTS.srfs.height.lowflat = { 0 };
% params.MULT_OPTS.srfs.height.flat = { 0 };
% params.MULT_OPTS.srfs.height.slant = { 0 };
% params.MULT_OPTS.srfs.height.multislant = { 0 };
% params.MULT_OPTS.srfs.height.multifront = { 0 };
% params.MULT_OPTS.srfs.height.frontness = { 0 };
% params.MULT_OPTS.srfs.height.contour = { 0 };
% params.MULT_OPTS.srfs.height.GKZ = { 0 };
% params.MULT_OPTS.srfs.height.GIKZ = { 0 };
% params.MULT_OPTS.srfs.height.MKZ = { 0 };
% params.MULT_OPTS.srfs.height.NKZ = { 0 };
% params.MULT_OPTS.srfs.height.NIKZ = { 0 };
% params.MULT_OPTS.srfs.height.KZ2 = { 0 };
% 
% params.MULT_OPTS.srfs.light.L2 = { 0 };
% params.MULT_OPTS.srfs.light.L1 = { 0 };
% 
% params.MULT_OPTS.srfs.albedo.range = { 0 };
% params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% params.MULT_OPTS.srfs.albedo.YC = { 0 };
% params.MULT_OPTS.srfs.albedo.MYC = { 0 };
% params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% params.MULT_OPTS.srfs.albedo.multi_entropy = {0};
% params.MULT_OPTS.srfs.albedo.gamut = {0};
% 
% params.MULT_OPTS.srfs.albedo.white1 = {0};
% params.MULT_OPTS.srfs.albedo.white2 = {0};
% params.MULT_OPTS.srfs.albedo.gaussian = {0};
% params.MULT_OPTS.srfs.albedo.hist = {0};
% params.MULT_OPTS.srfs.albedo.W_hist = {0};
% 
% 
% params.MULT_OPTS.srfs.light.L2 = { 0 };
% params.MULT_OPTS.srfs.light.L1 = { 0 };
% 
% 
% if EXPERIMENT_NUM == 101 % Best Model
%   
%   params.MULT_OPTS.srfs.height.slant = { 1 };
%   params.MULT_OPTS.srfs.height.contour = { 4 };
%   params.MULT_OPTS.srfs.height.GKZ = { 3 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 3 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
%   params.RENYI_SIGMA = {0.136310};
% 
%   params.ALPHA_LIGHT = {0.5}; % Correct
% 
%   N_ROUNDS = 2;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.albedo.multi_entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-1:.125:3]]];
%   
%   
% elseif EXPERIMENT_NUM == 601 % Best Model
%     
%   
%   params.MULT_OPTS.srfs.height.slant = { 1.5422 };
%   params.MULT_OPTS.srfs.height.contour = { 1.6818 };
%   params.MULT_OPTS.srfs.height.GKZ = { 1.6818 };
%   
%   params.MULT_OPTS.srfs.albedo.MYC = { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
%   params.MULT_OPTS.srfs.albedo.gamut = { 2.5937 };
%   params.MULT_OPTS.srfs.albedo.white1 = { 2.5937 };
%   
%   params.MULT_OPTS.srfs.light.L2 = { 1 };
% 
% 
%   params.RENYI_SIGMA = {0.136310};
%     
%   params.FALSE_COLOR_LIGHT = 1;
%   params.FORCE_MONOCHROME_LIGHT = 0;
%   params.LIGHT_COLLECTION = 'natural_train_color';
%     
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.L2',
%     'params.MULT_OPTS.srfs.albedo.white1',
% %     'params.MULT_OPTS.srfs.albedo.gamut',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.height.GKZ',
%     'params.MULT_OPTS.srfs.albedo.multi_entropy',
%     'params.MULT_OPTS.srfs.albedo.MYC',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-2:.125:4]]];
%   
% elseif EXPERIMENT_NUM == 611 % Best Model
%     
%   
% %   params.MULT_OPTS.srfs.height.slant = { 1.5422 };
% %   params.MULT_OPTS.srfs.height.contour = { 1.6818 };
% %   params.MULT_OPTS.srfs.height.MKZ = { 2.5937 };
% %   
% %   params.MULT_OPTS.srfs.albedo.MYC = { 6.1688 };
% %   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
% %   params.MULT_OPTS.srfs.albedo.gamut = { 2.5937 };
% %   params.MULT_OPTS.srfs.albedo.white1 = { 2.5937 };
% % 
% %   params.MULT_OPTS.srfs.light.L2 = { 1 };
% 
% params.MULT_OPTS.srfs.height.slant = { 1.5422 };
% params.MULT_OPTS.srfs.height.contour = { 1.6818 };
% params.MULT_OPTS.srfs.height.MKZ = { 2.5937 };
% 
% params.MULT_OPTS.srfs.albedo.MYC = { 5.6569 };
% params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
% params.MULT_OPTS.srfs.albedo.gamut = { 2.5937 };
% params.MULT_OPTS.srfs.albedo.white1 = { 2.8284 };
% 
% params.MULT_OPTS.srfs.light.L2 = { 1.2968 };
% 
% 
%   params.RENYI_SIGMA = {0.136310};
%     
%   params.FALSE_COLOR_LIGHT = 1;
%   params.FORCE_MONOCHROME_LIGHT = 0;
%   params.LIGHT_COLLECTION = 'natural_train_color';
%     
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 3;
%     
%   search_strings = {
%   'params.MULT_OPTS.srfs.height.MKZ',
%   'params.MULT_OPTS.srfs.height.slant',
%   'params.MULT_OPTS.srfs.albedo.MYC',
%   'params.MULT_OPTS.srfs.albedo.multi_entropy',
%   };
%   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.light.L2',
% %     'params.MULT_OPTS.srfs.albedo.white1',
% % %     'params.MULT_OPTS.srfs.albedo.gamut',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.height.GKZ',
% %     'params.MULT_OPTS.srfs.albedo.multi_entropy',
% %     'params.MULT_OPTS.srfs.albedo.MYC',
% %     };
%   
%   SEARCH_VALS = [0, [2.^[-2:.125:4]]];
%   
% elseif EXPERIMENT_NUM == 612 % Best Model
% 
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.32421 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 2.5937 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.3636 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
% 
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0.35355 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 3.3636 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 9.5137 };
% 
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 0.77111 };
% 
%   params.APYR_DEPTH = 1;
%   params.ZPYR_DEPTH = 1;
% 
%   params.RENYI_SIGMA = {0.136310};
%     
%   params.FORCE_MONOCHROME_LIGHT = 0;
%   params.LIGHT_COLLECTION = 'natural_train_color';
%     
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 2;
%     
%   search_strings = {
%   'params.MULT_OPTS.srfs.light.L2',
%   'params.MULT_OPTS.srfs.height.slant',
%   'params.MULT_OPTS.srfs.height.contour',
%   'params.MULT_OPTS.srfs.height.MKZ',
%   'params.MULT_OPTS.srfs.albedo.white1',
%   'params.MULT_OPTS.srfs.albedo.gray2',
%   'params.MULT_OPTS.srfs.albedo.gamut',
%   'params.MULT_OPTS.srfs.albedo.MA3',
%   'params.MULT_OPTS.srfs.albedo.multi_entropy',
%   };
%   
%   SEARCH_VALS = [0, [2.^[-4:.125:4]]];
%   
%   
% elseif EXPERIMENT_NUM == 613
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 1.1892 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 1.6818 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.3636 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 13.4543 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 0.77111 };
%   
%   params.APYR_DEPTH = {5};
%   params.ZPYR_DEPTH = {1};
%   
%   params.APYR_DEPTH_SMOOTH = {1};
%   params.APYR_DEPTH_ENTROPY = {1};
%   params.APYR_DEPTH_WHITE = {1};
% 
%   params.RENYI_SIGMA = {0.136310};
%     
%   params.FORCE_MONOCHROME_LIGHT = 0;
%   params.LIGHT_COLLECTION = 'natural_train_color';
%     
%   params.SOLVE_LIGHT = 0;
%   
%   N_ROUNDS = 2;
%     
%   search_strings = {
%   'params.MULT_OPTS.srfs.height.slant',
%   'params.MULT_OPTS.srfs.height.contour',
%   'params.MULT_OPTS.srfs.height.MKZ',
%   'params.MULT_OPTS.srfs.albedo.white1',
%   'params.MULT_OPTS.srfs.albedo.gray2',
%   'params.MULT_OPTS.srfs.albedo.gamut',
%   'params.MULT_OPTS.srfs.albedo.MA3',
%   'params.MULT_OPTS.srfs.albedo.multi_entropy',
%   };
%   
%   SEARCH_VALS = [0, [2.^[-4:.125:4]]];
%    
%   
% elseif EXPERIMENT_NUM == 614
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.3636 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 0.77111 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%  
%   params.N_LIGHTS = {1};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%     
%   search_strings = {
%     'params.N_LIGHTS'
%   };
%   
%   SEARCH_VALS = 1:15;
%   
% elseif EXPERIMENT_NUM == 621
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.3636 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 0.77111 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%  
%   params.N_LIGHTS = {1};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.L2'
%   };
%   
%   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% 
% elseif EXPERIMENT_NUM == 622
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.3636 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 0.77111 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%  
%   params.N_LIGHTS = {9};
%   params.SOLVE_LIGHT = 1;
%   params.LIGHT_ALPHA = {1};
%   N_ROUNDS = 2;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.L2'
%     'params.LIGHT_ALPHA'
%   };
%   
%   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
%   
% elseif EXPERIMENT_NUM == 623
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.3636 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 0.77111 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%  
%   params.N_LIGHTS = {15};
%   params.SOLVE_LIGHT = 1;
%   params.LIGHT_ALPHA = {1};
%   N_ROUNDS = 2;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.L2'
%     'params.LIGHT_ALPHA'
%   };
%   
% SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% 
% 
% 
% elseif EXPERIMENT_NUM == 630
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 1.0905 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.01};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {1.090500};
%   
%   params.N_LIGHTS = {9};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-6:.5:6]]];
%   
%   
% elseif EXPERIMENT_NUM == 631
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 1.0905 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.03};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {1.090500};
%   
%   params.N_LIGHTS = {9};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-6:.5:6]]];
%   
% 
% elseif EXPERIMENT_NUM == 632
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 1.0905 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.1};
% 
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {1.090500};
%   
%   params.N_LIGHTS = {9};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%   };
%   
%   SEARCH_VALS = [0, [2.^[-6:.5:6]]];
%   
% elseif EXPERIMENT_NUM == 633
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 1.0905 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.3};
% 
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {1.090500};
%   
%   params.N_LIGHTS = {9};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%   };
%   
%   SEARCH_VALS = [0, [2.^[-6:.5:6]]];
%   
% elseif EXPERIMENT_NUM == 634
%   
%   params.MULT_OPTS.srfs.height.lowflat =  { 0 };
%   params.MULT_OPTS.srfs.height.flat =     { 0 };
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.multislant =       { 0 };
%   params.MULT_OPTS.srfs.height.multifront =       { 0 };
%   params.MULT_OPTS.srfs.height.frontness =        { 0 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.GKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.GIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   params.MULT_OPTS.srfs.height.MIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.NKZ =      { 0 };
%   params.MULT_OPTS.srfs.height.NIKZ =     { 0 };
%   params.MULT_OPTS.srfs.height.KZ2 =      { 0 };
%   params.MULT_OPTS.srfs.height.dN =       { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range =    { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy =  { 0 };
%   params.MULT_OPTS.srfs.albedo.YC =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MYC =      { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB =      { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.MA =       { 0 };
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.white2 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian =         { 0 };
%   params.MULT_OPTS.srfs.albedo.gray =     { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 10.3747 };
%   
%   params.MULT_OPTS.srfs.light.L1 =        { 0 };
%   params.MULT_OPTS.srfs.light.L2 =        { 1.0905 };
%   params.MULT_OPTS.srfs.light.GSM =        { 0 };
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {1};
% 
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {1.090500};
%   
%   params.N_LIGHTS = {9};
%   params.SOLVE_LIGHT = 1;
%   
%   N_ROUNDS = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%   };
%   
%   SEARCH_VALS = [0, [2.^[-6:.5:6]]];
%   
%   
%   elseif EXPERIMENT_NUM == 701
%   
%     params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%     params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%     params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%     
%     params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%     params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%     params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%     params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%     params.MULT_OPTS.srfs.albedo.gray2 =    { 2.5937 };
%     
%     params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%     params.MULT_OPTS.srfs.light.lab_color.gaussian = { 7.336 };
%     
%     params.RENYI_SIGMA = {0.136310};
%     params.LIGHT_ALPHA = {1};
%     
%     params.N_LIGHTS = {9};
%     params.DO_LIGHT_POSTERIOR = {1};
%     
%     params.PA_EPSILON = {0.001};
%     params.L_D_MIN = {0.01};
%     
%     N_ROUNDS = 1;
%     params.SOLVE_LIGHT = 1;
%     
%     search_strings = {
%       'params.PA_EPSILON',
%       'params.L_D_MIN'
%       };
%     
%     SEARCH_VALS = [0, [2.^[-16:.5:1]]];
%     params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   elseif EXPERIMENT_NUM == 702
%   
%     params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%     params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%     params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%     
%     params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%     params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%     params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%     params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%     params.MULT_OPTS.srfs.albedo.gray2 =    { 2.5937 };
%     
%     params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%     params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%     
%     params.RENYI_SIGMA = {0.136310};
%     params.LIGHT_ALPHA = {0};
%     
%     params.N_LIGHTS = {9};
%     params.DO_LIGHT_POSTERIOR = {1};
%     
%     params.PA_EPSILON = {0.001};
%     params.L_D_MIN = {0.01};
%     
%     params.ZPYR_MULT = {2};
%     params.L_MULT = {1};
%     
%     N_ROUNDS = 2;
%     params.SOLVE_LIGHT = 1;
%       
%     search_strings = {
%       'params.MULT_OPTS.srfs.albedo.multi_entropy',
%       'params.MULT_OPTS.srfs.albedo.white1',
%       'params.MULT_OPTS.srfs.albedo.gray2',
%       'params.MULT_OPTS.srfs.albedo.gamut',
%       'params.MULT_OPTS.srfs.albedo.MA3',
%       'params.MULT_OPTS.srfs.height.slant',
%       'params.MULT_OPTS.srfs.height.contour',
%       'params.MULT_OPTS.srfs.height.MKZ',
%       'params.MULT_OPTS.srfs.light.natural_color.gaussian',
%       'params.MULT_OPTS.srfs.light.lab_color.gaussian',
%       };
%         
%     SEARCH_VALS = [0, [2.^[-2:.125:4]]];
%     params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%     
% elseif EXPERIMENT_NUM == 703
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 0 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 2;
%   
%   N_ROUNDS = 2;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.hist',
%     'params.MULT_OPTS.srfs.albedo.gamut',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-4:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
% elseif EXPERIMENT_NUM == 711
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   N_ROUNDS = 2;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     'params.MULT_OPTS.srfs.height.MKZ',
%     'params.MULT_OPTS.srfs.light.natural_color.gaussian',
%     'params.MULT_OPTS.srfs.light.lab_color.gaussian',
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.albedo.multi_entropy',
%     'params.MULT_OPTS.srfs.albedo.hist',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-2:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
%   
%     
% elseif EXPERIMENT_NUM == 712
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.L_D_MIN'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-16:.5:1]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   elseif EXPERIMENT_NUM == 713
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-2:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
%   
% elseif EXPERIMENT_NUM == 714
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.5};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.5:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
% elseif EXPERIMENT_NUM == 715
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.3};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.5:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 716
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.5:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
% elseif EXPERIMENT_NUM == 721
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.5};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     };
%   
%   SEARCH_VALS = [0, [2.^[2:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
% elseif EXPERIMENT_NUM == 722
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.3};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   params.A_MEDIAN_HALFWIDTH = 2;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     };
%   
%   SEARCH_VALS = [0, [2.^[2:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 723
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   params.A_MEDIAN_HALFWIDTH = 3;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     };
%   
%   SEARCH_VALS = [0, [2.^[1:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 724
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   params.A_MEDIAN_HALFWIDTH = 4;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     };
%   
%   SEARCH_VALS = [0, [2.^[1:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
%   
%   
%   
% elseif EXPERIMENT_NUM == 731
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.5};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.MKZ',
%     };
%   
%   SEARCH_VALS = [0, [2.^[1:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
% elseif EXPERIMENT_NUM == 732
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.3};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 2;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.MKZ',
%     };
%   
%   SEARCH_VALS = [0, [2.^[1:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 733
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.A_MEDIAN_HALFWIDTH = 4;
%   params.Z_MEDIAN_HALFWIDTH = 3;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.MKZ',
%     };
%   
%   SEARCH_VALS = [0, [2.^[0:.125:2]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 734
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 2.181 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.A_MEDIAN_HALFWIDTH = 4;
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%     
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.MKZ',
%     };
%   
%   SEARCH_VALS = [0, [2.^[0:.125:2]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
% elseif EXPERIMENT_NUM == 751
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 0 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 0 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   params.A_MEDIAN_HALFWIDTH = 4;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.W_hist',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-1:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
%     
% elseif EXPERIMENT_NUM == 761
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 0 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 4.362 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 1;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.W_hist',
% %     'params.MULT_OPTS.srfs.height.MKZ',
% %     'params.MULT_OPTS.srfs.albedo.MA3',
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     'params.MULT_OPTS.srfs.albedo.multi_entropy',
% %     'params.MULT_OPTS.srfs.light.natural_color.gaussian'
% %     'params.MULT_OPTS.srfs.light.lab_color.gaussian'
%     };
%   
%   SEARCH_VALS = [[2.^[1:.125:3]]];
% %   SEARCH_VALS = [0, [2.^[-2:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 762
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 0 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 4.362 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 2;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.W_hist',
%     };
%   
%   SEARCH_VALS = [[2.^[1:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 763
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 0 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 4.362 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 3;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.W_hist',
%     };
%   
%   SEARCH_VALS = [[2.^[1:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
% elseif EXPERIMENT_NUM == 764
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 0 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 0 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 0 };
%   params.MULT_OPTS.srfs.albedo.hist =   { 0 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 4.362 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy = {0};
%   params.MULT_OPTS.srfs.light.entropy_sigma = {.8};
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   params.ALBEDO_HIST_MODEL = 4;
%   
%   params.Z_MEDIAN_HALFWIDTH = 1;
%   params.A_MEDIAN_HALFWIDTH = 1;
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.W_hist',
%     };
%   
%   SEARCH_VALS = [[2.^[1:.125:3]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 771
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 8 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy_sigma =     { 0.5 };
%   params.MULT_OPTS.srfs.light.entropy =   { 0 };%{ 64 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.Z_MEDIAN_HALFWIDTH = 2;
%   params.A_MEDIAN_HALFWIDTH = 2;
%   
%   params.ALBEDO_HIST_MODEL = 2;
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.MKZ',
%     'params.MULT_OPTS.srfs.albedo.MA3',
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.light.natural_color.gaussian',
%     'params.MULT_OPTS.srfs.light.lab_color.gaussian'
%     'params.MULT_OPTS.srfs.albedo.W_hist',
%     'params.MULT_OPTS.srfs.albedo.multi_entropy',
%     };
%   
%   SEARCH_VALS = [[2.^[-2:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 781
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 5.6569 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy_sigma =     { 0.3 };
%   params.MULT_OPTS.srfs.light.entropy =   { 0 };%{ 64 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   params.A_MEDIAN_HALFWIDTH = 4;
%   
%   params.ALBEDO_HIST_MODEL = 2;
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%    
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.5:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
% elseif EXPERIMENT_NUM == 782
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 5.6569 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy_sigma =     { 0.5 };
%   params.MULT_OPTS.srfs.light.entropy =   { 0 };%{ 64 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   params.A_MEDIAN_HALFWIDTH = 4;
%   
%   params.ALBEDO_HIST_MODEL = 2;
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%    
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.5:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
% elseif EXPERIMENT_NUM == 783
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 2.8284 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.W_hist =   { 5.6569 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 12.3377 };
%   
%   params.MULT_OPTS.srfs.light.entropy_sigma =     { 0.7 };
%   params.MULT_OPTS.srfs.light.entropy =   { 0 };%{ 64 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {0};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.Z_MEDIAN_HALFWIDTH = 4;
%   params.A_MEDIAN_HALFWIDTH = 4;
%   
%   params.ALBEDO_HIST_MODEL = 2;
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   N_ROUNDS = 1;
%   params.SOLVE_LIGHT = 1;
%    
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.5:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
% elseif EXPERIMENT_NUM == 704
%   
%   params.MULT_OPTS.srfs.height.slant =    { 0.64842 };
%   params.MULT_OPTS.srfs.height.contour =  { 0.54525 };
%   params.MULT_OPTS.srfs.height.MKZ =      { 3.668 };
%   
%   params.MULT_OPTS.srfs.albedo.MA3 =      { 6.1688 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy =    { 0.77111 };
%   params.MULT_OPTS.srfs.albedo.gamut =    { 1.1892 };
%   params.MULT_OPTS.srfs.albedo.white1 =   { 1.6818 };
%   params.MULT_OPTS.srfs.albedo.gray2 =    { 2.5937 };
%   
%   params.MULT_OPTS.srfs.light.natural_color.gaussian = { 2.8284 };
%   params.MULT_OPTS.srfs.light.lab_color.gaussian = { 7.336 };
%   
%   params.RENYI_SIGMA = {0.136310};
%   params.LIGHT_ALPHA = {1};
%   
%   params.N_LIGHTS = {9};
%   params.DO_LIGHT_POSTERIOR = {1};
%   
%   params.PA_EPSILON = {0.001};
%   params.L_D_MIN = {0.01};
%   
%   params.ZPYR_MULT = {2};
%   params.L_MULT = {1};
%   
%   N_ROUNDS = 2;
%   params.SOLVE_LIGHT = 1;
%   
%   params.MULT_OPTS.srfs.light.entropy_sigma =     { 0.3 };
%   params.MULT_OPTS.srfs.light.entropy =   { 0 };%{ 64 };
% 
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.entropy',
%     'params.MULT_OPTS.srfs.light.entropy_sigma'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.125:4]]];
%   params.EVAL_NAMES = [MIT_TRAIN_EXPAND({'L'}), MIT_TRAIN_EXPAND({'C'})];
%   
%   
%   
% elseif EXPERIMENT_NUM == 102 % Best Model
%   
%   params.MULT_OPTS.srfs.height.slant = { 1 };
%   params.MULT_OPTS.srfs.height.contour = { 4 };
%   params.MULT_OPTS.srfs.height.GKZ = { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 3 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
%   params.RENYI_SIGMA = {0.136310};
%   
%   params.ALPHA_LIGHT = {0.5}; % Correct
%   
%   N_ROUNDS = 2;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.albedo.multi_entropy'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-1:.125:3]]];
%   
%     
% elseif EXPERIMENT_NUM == 103 % Best Model
%   
%   params.MULT_OPTS.srfs.height.slant = { 2.8284 };
%   params.MULT_OPTS.srfs.height.contour = { 2 };
%   params.MULT_OPTS.srfs.height.GKZ = { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 3 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
%   params.MULT_OPTS.srfs.albedo.gamut = { 0 };
%   params.RENYI_SIGMA = {0.136310};
%   
%   params.ALPHA_LIGHT = {0.5}; % Correct
%   
%   N_ROUNDS = 2;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.albedo.multi_entropy'
%     };
%   
%     SEARCH_VALS = [0, [2.^[-1:.125:3]]];
%   
%     
%     
% elseif EXPERIMENT_NUM == 105 
%   
%   params.MULT_OPTS.srfs.height.slant = { 3.3636 };
%   params.MULT_OPTS.srfs.height.contour = { 4 };
%   params.MULT_OPTS.srfs.height.GKZ = { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 3 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
%   params.MULT_OPTS.srfs.albedo.gamut = { 0 };
%   params.RENYI_SIGMA = {0.136310};
%   
%   params.GAMUT_POWER = 1;
%   
%   N_ROUNDS = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.gamut',
%     };
%   
%     SEARCH_VALS = [0, [2.^[-6:.25:6]]];
%     
% elseif EXPERIMENT_NUM == 106
%   
%   params.MULT_OPTS.srfs.height.slant = { 3.3636 };
%   params.MULT_OPTS.srfs.height.contour = { 4 };
%   params.MULT_OPTS.srfs.height.GKZ = { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 3 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 1 };
%   params.MULT_OPTS.srfs.albedo.gamut = { 0 };
%   params.RENYI_SIGMA = {0.136310};
%   
%   params.ALPHA_LIGHT = {0.5}; % Correct
%   
%   params.GAMUT_POWER = 1;
%   
%   N_ROUNDS = 1;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.gamut',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-12:.5:12]]];
%     
% elseif EXPERIMENT_NUM == 111 % Best Model
%   
%   params.MULT_OPTS.srfs.height.slant = { 2.8284 };
%   params.MULT_OPTS.srfs.height.contour = { 4.7568 };
%   params.MULT_OPTS.srfs.height.GKZ = { 1.4142 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 4 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 0 };
% 
%   params.RENYI_SIGMA = {0.3};
%   
% 
%   N_ROUNDS = 3;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.albedo.multi_entropy'
%     'params.RENYI_SIGMA'
%     };
%   
%   SEARCH_VALS = [0, [2.^[-3:.125:2]]];
%  
%   
%   
% elseif EXPERIMENT_NUM == 104 % Simple color model
%   
%   params.MULT_OPTS.srfs.height.slant = { 1 };
%   params.MULT_OPTS.srfs.height.contour = { 4 };
%   params.MULT_OPTS.srfs.height.GKZ = { 2 };
%   
%   params.MULT_OPTS.srfs.albedo.YC = { 0 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB = { 3 };
%   params.MULT_OPTS.srfs.albedo.entropy = { 1 };
%   
% %   params.MULT_OPTS.srfs.albedo.gamut = { 1 };
% %   
% %   params.MULT_OPTS.srfs.albedo.white1 = {1};
% %   params.MULT_OPTS.srfs.albedo.white2 = {0};
% %   params.MULT_OPTS.srfs.albedo.gaussian = {0};
% %   
% %   params.MULT_OPTS.srfs.light.L2 = { 1 };
% %   params.MULT_OPTS.srfs.light.L1 = { 0 };
% 
%   
% params.MULT_OPTS.srfs.albedo.gamut = { 0 };
% 
% params.MULT_OPTS.srfs.albedo.white1 = {0};
% params.MULT_OPTS.srfs.albedo.white2 = {0};
% params.MULT_OPTS.srfs.albedo.gaussian = {0};
% 
% params.MULT_OPTS.srfs.light.L2 = { 1 };
% params.MULT_OPTS.srfs.light.L1 = { 0 };
% 
% params.SOLVE_LIGHT = 0;
% 
%   params.RENYI_SIGMA = {0.3};
% 
%   params.ENTROPY_PROJECTIONS_ORDER = 1;
%   params.USE_APPROX_RENYI = 1;
%   
%   N_ROUNDS = 2;
%   
%   search_strings = {
% %     'params.MULT_OPTS.srfs.light.L2',
% %     'params.MULT_OPTS.srfs.albedo.white1',
% %     'params.MULT_OPTS.srfs.albedo.gamut',
%     'params.MULT_OPTS.srfs.height.contour',
%     'params.MULT_OPTS.srfs.height.slant',
%     'params.MULT_OPTS.srfs.albedo.entropy'
% %     'params.MULT_OPTS.srfs.height.GKZ',
% %     'params.MULT_OPTS.srfs.albedo.RGB',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-1:.125:4]]];
%   
%   
% elseif EXPERIMENT_NUM == 204 % Simple color model
%   
%   
%   params.MULT_OPTS.srfs.height.lowflat = { 0 };
%   params.MULT_OPTS.srfs.height.flat = { 0 };
%   params.MULT_OPTS.srfs.height.slant = { 1.4142 };
%   params.MULT_OPTS.srfs.height.multislant = { 0 };
%   params.MULT_OPTS.srfs.height.multifront = { 0 };
%   params.MULT_OPTS.srfs.height.frontness = { 0 };
%   params.MULT_OPTS.srfs.height.contour = { 0.70711 };
%   params.MULT_OPTS.srfs.height.GKZ = { 2 };
%   params.MULT_OPTS.srfs.height.GIKZ = { 0 };
%   
%   params.MULT_OPTS.srfs.albedo.range = { 0 };
%   params.MULT_OPTS.srfs.albedo.entropy = { 0.917 };
%   params.MULT_OPTS.srfs.albedo.YC = { 0 };
%   params.MULT_OPTS.srfs.albedo.RGB = { 3 };
%   params.MULT_OPTS.srfs.albedo.multi_entropy = { 0 };
%   params.MULT_OPTS.srfs.albedo.gamut = { 1 };
%   params.MULT_OPTS.srfs.albedo.white1 = { 1 };
%   params.MULT_OPTS.srfs.albedo.white2 = { 0 };
%   params.MULT_OPTS.srfs.albedo.gaussian = { 0 };
%   
%   params.MULT_OPTS.srfs.light.L1 = { 0 };
%   params.MULT_OPTS.srfs.light.L2 = { 1 };
%   
%   params.SOLVE_LIGHT = 0;
% 
%   params.RENYI_SIGMA = {0.3};
% 
%   params.ENTROPY_PROJECTIONS_ORDER = 1;
%   params.USE_APPROX_RENYI = 1;
%   
%   N_ROUNDS = 2;
%   
%   search_strings = {
%     'params.MULT_OPTS.srfs.light.L2',
%     'params.MULT_OPTS.srfs.albedo.gamut',
%     'params.MULT_OPTS.srfs.albedo.white1',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     'params.MULT_OPTS.srfs.height.GKZ',
% %     'params.MULT_OPTS.srfs.albedo.RGB',
%     };
%   
%   SEARCH_VALS = [0, [2.^[-1:.125:4]]];
%   
%   
% else 
%   assert(1==0);
% end
% 
% 
% % if EXPERIMENT_NUM == 141
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 7.336 };
% %   params.MULT_OPTS.srfs.height.contour = { 5.6569 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 2.181 };
% %   
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.range = { 0 };
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %     
% %   params.RENYI_SIGMA = {0.75};
% %   params.ADAPTIVE_BETA = {0};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.ADAPTIVE_BETA'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-16:.5:0]]];
% %   
% % elseif EXPERIMENT_NUM == 142
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 6 };
% %     
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.range = { 0 };
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %     
% %   params.RENYI_SIGMA = {0.75};
% %   params.ADAPTIVE_BETA = {0};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.ADAPTIVE_BETA'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-16:.5:0]]];
% %   
% % else
% %   assert(1==0)
% % end
% 
% 
% 
% % if EXPERIMENT_NUM == 101
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 121
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 2;
% %   
% %   params.CURVATURE_MODE = 'IK';
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 122
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 2;
% %   
% %   params.CURVATURE_MODE = 'K';
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:4]]];
% %   
% % elseif EXPERIMENT_NUM == 123
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   params.CURVATURE_MODE = 'IK';
% %   
% %   params.F_PERCENT = .001;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 124
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   params.CURVATURE_MODE = 'K';
% %   
% %   params.F_PERCENT = .001;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:4]]];
% %   
% % elseif EXPERIMENT_NUM == 125
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 2;
% %   
% %   params.CURVATURE_MODE = 'IK';
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.height.GKZ',
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 126
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 2;
% %   
% %   params.CURVATURE_MODE = 'K';
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.height.GKZ',
% %     'params.MULT_OPTS.srfs.height.slant',
% %     'params.MULT_OPTS.srfs.height.contour',
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:4]]];
% %   
% %   
% % elseif EXPERIMENT_NUM == 102
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model2';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% %   
% % elseif EXPERIMENT_NUM == 103
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model3';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 104
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model4';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 105
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model4';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.RGB'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 106
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model5';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 107
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model6';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% % elseif EXPERIMENT_NUM == 108
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model7';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% %   
% % elseif EXPERIMENT_NUM == 109
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 0 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model8';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.YC'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-1:.125:5]]];
% %   
% %   
% % elseif EXPERIMENT_NUM == 111
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.75};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% %   
% % elseif EXPERIMENT_NUM == 112
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.5};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% %   
% % elseif EXPERIMENT_NUM == 113
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.25};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% %   
% %   
% % elseif EXPERIMENT_NUM == 114
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.1};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% %   
% % elseif EXPERIMENT_NUM == 115
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.05};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% %   
% % elseif EXPERIMENT_NUM == 116
% %   
% %   
% %   params.MULT_OPTS.srfs.height.slant = { 4.362 };
% %   params.MULT_OPTS.srfs.height.contour = { 7.336 };
% %   params.MULT_OPTS.srfs.height.GKZ = { 4 };
% %   
% %   params.MULT_OPTS.srfs.albedo.entropy = { 0 };
% %   params.MULT_OPTS.srfs.albedo.YC = { 8.7241 };
% %   params.MULT_OPTS.srfs.albedo.RGB = { 0 };
% %   
% %   params.RENYI_SIGMA = {0.025};
% %   
% %   params.YC_MODEL = 'model1';
% %   N_ROUNDS = 1;
% %   
% %   search_strings = {
% %     'params.MULT_OPTS.srfs.albedo.entropy'
% %     };
% %   
% %   SEARCH_VALS = [0, [2.^[-2:.125:3]]];
% %   
% % else
% %   assert(1==0)
% % end
% 
% 
% 
% NAME = ['exp', num2str(EXPERIMENT_NUM), '.'];
% 
% 
% search_strings_expand = {};
% for i = 1:length(search_strings)
%   x = eval(search_strings{i});
%   for j = length(x):-1:1
%     s = [search_strings{i}, '{', num2str(j), '}'];
%     search_strings_expand{length(search_strings_expand)+1} = s;
%   end
% end
% search_strings = search_strings_expand;
% 
% search_strings_sorted = {};
% for j = length(x):-1:1
%   take = cellfun(@(x) ~isempty(strfind(x, num2str(j))), search_strings, 'UniformOutput', false);
%   take = cat(2,take{:});
%   search_strings_sorted = [search_strings_sorted, search_strings(take)];
% end
% search_strings = search_strings_sorted;
% 
% cellfun(@(x) fprintf([x, '\n']), search_strings)
% 
% [boosted_params, history_err, history_param] = coordExperiment(NAME, params, search_strings, SEARCH_VALS, N_ROUNDS, OBJECTIVES);
% 
