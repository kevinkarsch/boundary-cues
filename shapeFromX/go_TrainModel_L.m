
clear all;

REGENERATE_DATA = 1;

rand('twister',5489)
randn('state',0)

if REGENERATE_DATA
  
  names = MIT_TRAIN;
  
  Ls = {};
  for name_i = 1:length(names)
    
    name = names{name_i};
    load(['data/MITC-intrinsic/data/', name, '/L.mat'])
    L = cat(3,color_lights.light{:}, color_lights.diffuse, color_lights.shading);
    Ls{name_i} = L;
  end
  Ls = cat(3,Ls{:});
  
  % Ls = Ls(:,Ls(3,:) > 0);
  
  valid = false(1, size(Ls,3));
  for i = 1:size(Ls,3)
    valid(i) = validSH(Ls(:,:,i));
  end
  
  % figure; visLights(Ls(:,:,valid)); title('keep'); drawnow;
  % figure; visLights(Ls(:,:,~valid)); title('discard'); drawnow;
  
  Ls = Ls(:,:,valid);
  
  Ldata.train_color = Ls;
  Ldata.train_gray = mean(Ls,2);
  
  
  load natural_lights
  
  
  
  Ls_flip = Ls;
  Ls_flip([2,5,6],:) = -Ls_flip([2,5,6],:);
  
  Ls = cat(3, Ls, Ls_flip);
  
  Ls_yuv_base = L_rgb2yuv(Ls);
  
  Ls_yuv = {};
  for c = [2, 2.5, 3]
    Ls_yuv{end+1} = Ls_yuv_base .* repmat([1, 2.^-(c-1), 2.^-(c-1)], [9, 1, size(Ls_yuv_base,3)]);
  end
  Ls_yuv = cat(3,Ls_yuv{:});
  
  Ls = L_yuv2rgb(Ls_yuv);
  
  mults = [2.5,3,3.5];%[1, 1.5, 2, 2.5];
  Ls2 = {};
  for m = mults
    Ls2{end+1} = Ls*m;
  end
  Ls = cat(3,Ls2{:});
  
  
  Ls2 = {};
  for i = 1:size(Ls,3)
    L = Ls(:,:,i);
    if validSH(L)
      V = log(visSH_color(Ls(:,:,i)));
      c4 = 0.886227;
      L(1,:) = L(1,:) - max(V(:))/c4;
      V2 = log(visSH_color(L));
      Ls2{end+1} = L;
    end
  end
  Ls2 = cat(3,Ls2{:});
  Ls = Ls2;
  
  
  X = {};
  for i = 1:size(Ls,3)
    v = log(visSH_color(Ls(:,:,i), [32,32]));
    v = v - max(v(:));
    X{i} = v(~isnan(v));
  end
  X = cat(2,X{:});
  Dsq = distMat(X', X');
  
  too_close = triu(Dsq < 2,1);
  keep = ~any(too_close,1);
  Ls = Ls(:,:,keep);
  
  
  
  rand('twister',5489)
  randn('state',0)
  ridx = randperm(size(Ls,3));
  idx_train = ridx(1:round(size(Ls,3)/2));
  idx_test = ridx((round(size(Ls,3)/2)+1):end);
  
  Ls_train = Ls(:,:,idx_train);
  Ls_test = Ls(:,:,idx_test);
  
  
  Ldata.natural_train_color = Ls_train;
  Ldata.natural_test_color = Ls_test;
  
  save Ldata
else
  load Ldata
end


Ls = Ldata.natural_train_color;
X = reshape(Ls, [], size(Ls,3))';
[Xw, L_whiten_params] = whiten(X, 1, 0);
Ldata.whiten.natural_color = L_whiten_params;
Ldata.gaussian.natural_color.mu = mean(X,1);
Ldata.gaussian.natural_color.Sigma = cov(X);
Ldata.X.natural_color = X;


Ls = Ldata.train_color;
X = reshape(Ls, [], size(Ls,3))';

mu = mean(X);
Sigma = cov(X);

LL = lmvnpdf(X, mu, Sigma);
LL = LL - max(LL(:));
X = X(LL > prctile(LL, 50),:);
Ls = Ls(:,:,LL > prctile(LL, 50),:);

p = perms(1:3);
Ls2 = {};
for i = 1:size(p,1)
  Ls2{i} = Ls(:,p(i,:),:);
end
Ls = cat(3,Ls2{:});

Lf = Ls;
Lf([2,5,6],:,:) = -Lf([2,5,6],:,:);
Ls = cat(3, Ls, Lf);

V = {};
for i = 1:size(Ls,3)
  v = visSH_color(Ls(:,:,i), [150,150]);
  v = v ./ max(v(:));
  v(isnan(v)) = 0;
  V{i} = v;
end
V = cat(4,V{:});
montage(V)

X = reshape(Ls, [], size(Ls,3))';
[Xw, L_whiten_params] = whiten(X, 1, 0);
Ldata.whiten.lab_color = L_whiten_params;
Ldata.gaussian.lab_color.mu = mean(X,1);
Ldata.gaussian.lab_color.Sigma = cov(X);
Ldata.X.lab_color = X;



save Ldata.mat Ldata

