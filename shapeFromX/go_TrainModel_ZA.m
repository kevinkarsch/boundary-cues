clear all;

rand('twister',5489)
randn('state',0)

addpath(genpath('./matlabPyrTools/'));
addpath('./bmorph/');

CONSTANTS;

DEPTH = 7;
N_GAUSSIANS = 40;
N_TRAIN = 100000;

in_directory = ['data/MITC-Intrinsic/data/'];

% MODE = 'test'
% MODE = 'train';

out_file = ['./prior.mat'];
names = MIT_TRAIN;
% names = [MIT_TRAIN, MIT_TEST];

try
  load(out_file, 'prior');
  fprintf('Loaded existing prior file.\n');
catch
  prior = [];
  fprintf('No existing prior file, making a new one.\n');
end

% load texton_fb.mat
% fb = fb(:);

KZpyrs = {};
GApyrs = {};
NKZpyrs2 = {};
KZpyrs2 = {};
GApyrs2 = {};
Vs = {};

Ipyrs = {};
Apyrs = {};
Spyrs = {};
Zpyrs = {};
Vpyrs = {};
Tpyrs = {};
pinds = {};

As = {};
Is = {};
Zs = {};
Ns = {};
GZs = {};
Zsigmas = [];
props = [];
MKZs = {};
MAs = {};
Nzs = {};

I_Ygs = {};
I_Cgs = {};

HALF_WIDTH = 2;

for name_i = 1:length(names)
  
  name = names{name_i};
%   name = 'frog1';
  
  fprintf('loading %s\n', name)
  
  A = imread([in_directory, name, '/reflectance.png']);
  A = double(A) ./ double(intmax('uint16'));
  background = any(A == 0,3);
  A = log(max(eps,A));
  A(repmat(background, [1,1,3])) = nan;
    
  I = imread([in_directory, name, '/diffuse.png']);
  I = double(I) ./ double(intmax('uint16'));
  background = any(I == 0,3);
  I = log(max(eps,I));
  I(repmat(background, [1,1,3])) = nan;
  
  load([in_directory, name, '/Z.mat']);
  Z = depth; clear depth;
  
  KZ = getK_fast(Z);
  
  V = ~isnan(KZ) & ~background;
  
  A(repmat(~V, [1,1,3])) = nan;
  Av = reshape(A, [], 3);

  As{name_i} = Av(~any(isnan(Av),2),:);
  
  imagesc([visualizeDEM(Z), min(1,exp(A))]); imtight; colormap('gray'); drawnow;
  
%   YUV = rgb2yuv(A);
%   
%   YG = getGradNorm(YUV(:,:,1));
%   UG = getGradNorm(YUV(:,:,2));
%   VG = getGradNorm(YUV(:,:,3));
%   YUVs{name_i} = [YG(:), UG(:), VG(:)];
%   YUVs{name_i} = YUVs{name_i}(all(~isnan(YUVs{name_i}),2),:);
%   
%   RG = getGradNorm(A(:,:,1));
%   GG = getGradNorm(A(:,:,2));
%   BG = getGradNorm(A(:,:,3));
%   
%   RGBs{name_i} = [RG(:), GG(:), BG(:)];
%   RGBs{name_i} = RGBs{name_i}(all(~isnan(RGBs{name_i}),2),:);
  
  M = medianFilterMat_mask(~V, HALF_WIDTH);
  
  MAs{name_i} = M * Av;
  MKZs{name_i} = M * KZ(:);
  
end




As = cellfun(@(x) randomlySelect(x, 30000), As, 'UniformOutput', false);
A = cat(1, As{:});

% YUV = cellfun(@(x) randomlySelect(x, 30000), YUVs, 'UniformOutput', false);
% YUV = cat(1, YUV{:});
% for s = 1:3
%   prior.albedo.YUV.GSMs{s} = GSM_fit(YUV(:,s), 40, 1);
% end
% 
% 
% RGB = cellfun(@(x) randomlySelect(x, 30000), RGBs, 'UniformOutput', false);
% RGB = cat(1, RGB{:});
% prior.albedo.RGB.GSM = GSM_fit(RGB(:), 40, 1);

[A_white, whiten_params] = whiten(A, 0, 0);
prior.albedo.A_whiten = whiten_params.map;

As_white = {};
for i = 1:length(As)
  As_white{i} = As{i} * whiten_params.map;
end
A_white = A * whiten_params.map;



bin_low =  params.ALBEDO_BIN_LOW;
bin_high = params.ALBEDO_BIN_HIGH;

lambdas = 2.^[4:.5:10];
robust_costs = 1;%[0,1];
LLs = [];
for ri = 1:length(robust_costs)
  for li = 1:length(lambdas)
    for j = 1:length(As_white)
      A_train = cat(1,As_white{setdiff(1:length(As_white), j)});
      A_test = As_white{j};
      [junk, LL] = smoothHist3_fit(A_train, bin_low, bin_high, lambdas(li), robust_costs(ri));
      
      s = splat3_fast_wrapper(A_test, bin_low, bin_high);
      LLs(ri,li,j) = sum(sum(sum(LL .* s.N)));
    end
  end
end

cost = -sum(LLs,3);
[ri, li] = find(cost == min(cost(:)));
lambda = lambdas(li)
robust_cost = robust_costs(ri)

A_train = cat(1,As_white{:});
[prior.albedo.Aw_hist, LL] = smoothHist3_fit(A_train, bin_low, bin_high, lambda, robust_cost);


MAs = cellfun(@(x) randomlySelect(x, 30000), MAs, 'UniformOutput', false);
MA = cat(1, MAs{:});

prior.albedo.MA.GSM_mvn = GSM_mvn_fitsmart(MA, 40);


MKZ = MKZs(:);
MKZ = cellfun(@(x) randomlySelect(x, 30000), MKZ, 'UniformOutput', false);
MKZ = cat(1, MKZ{:});
prior.height.MKZ.GSM = GSM_fit(MKZ, 40, 0);


save(out_file, 'prior');

