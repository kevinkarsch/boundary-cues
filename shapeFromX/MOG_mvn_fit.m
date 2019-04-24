function [model, score] = MOG_mvn_fit(data, K, DISPLAY)
if nargin < 3
  DISPLAY = 0;
end

LOSS_FRACTION_THRESH = 1*10^-8;
N_ITERS = 500;
N_VIS_BINS = min(100, sqrt(size(data,1))/2);
N_DATA_BINS = 1000;
MIN_COV = 10^-13;
MIN_ITERS = 20;

% if all(data(:,2) == 0)
%   data(:,2) = data(:,2) + normrnd(0, 10^-8, size(data(:,2)));
% end
% 
% if all(data(:,1) == 0)
%   data(:,1) = data(:,1) + normrnd(0, 10^-8, size(data(:,1)));
% end


r = rand(size(data,1), 1);
idx = ceil(r * K);

model.pis = nan(1,K);
model.mus = nan(K, size(data,2));
model.Sigmas = nan(size(data,2), size(data,2), K);

for k = 1:K
  model.pis(k) = mean(idx == k);
  model.mus(k,:) = mean(data(idx == k,:));
  model.Sigmas(:,:,k) = cov(data(idx == k,:));
end

% model.pis = ones(1,K) / K;
% model.mus = repmat((mean(data) + median(data))/2, K, 1);
% model.Sigmas = repmat(cov(data)/2, [1, 1, K]) .* repmat(reshape((2.^([0:(K-1)] - (K-1)/2)), [1, 1, K]), 2, 2);
% model.Sigmas = max(model.Sigmas, repmat(MIN_COV*eye(2,2), [1, 1, K]));

score_bak = inf;
start_time = clock;


% data_step = 100/N_DATA_BINS;
% data_X = prctile(data(:,1), [data_step/2:data_step:(100-data_step/2)])';
% data_Y = prctile(data(:,2), [data_step/2:data_step:(100-data_step/2)])';
% data_N = hist3(data, {data_X, data_Y});
% [data_XY_X, data_XY_Y] = meshgrid(data_X, data_Y);
% 
% 
% valid = data_N>0;
% data_XY = [data_XY_X(valid), data_XY_Y(valid)];
% data_N = data_N(data_N>0);


extents = { [mean(data(:,1)) + [-2.5, 2.5]*std(data(:,1))], [mean(data(:,2)) + [-2.5, 2.5]*std(data(:,2))]};
steps = [((extents{1}(2) - extents{1}(1))/N_VIS_BINS), ((extents{2}(2) - extents{2}(1))/N_VIS_BINS)];
steps = max(steps, eps);
ranges = {...
  (extents{1}(1)-eps):steps(1):(extents{1}(2)+eps), ...
  (extents{2}(1)-eps):steps(2):(extents{2}(2)+eps)};
data_pixel_area = steps(1)*steps(2);

% [vis_N, vis_X] = hist3(data, ranges);
% vis_N = vis_N./sum(vis_N(:));
% 
% [vis_i, vis_j] = find(ones(size(vis_N)));
% vis_XY = [vis_X{1}(vis_i)', vis_X{2}(vis_j)'];


% Ps = nan(size(data,1), K);
Ps = nan(size(data,1), K);
% Ps_vis = nan([size(vis_N), K]);
DO_BREAK = false;
for iter = 1:N_ITERS
  
  % E-step
  for k = 1:K
    
%     data_diff = [data(:,1) - model.mus(k,1), data(:,2) - model.mus(k,2)];
    data_diff = data - repmat(model.mus(k,:), size(data,1), 1);
    Sigma_pad = model.Sigmas(:,:,k);
    iSigma = inv(Sigma_pad);
    Ps(:,k) = exp(-0.5 * sum((data_diff * iSigma) .* data_diff,2)) * (model.pis(k) ./  ((2*pi).^(size(data,2)/2) * sqrt(det(Sigma_pad))));
    
    if any(imag(Ps(:,k)))
      keyboard
    end

  
%     data_diff = [data(:,1) - model.mus(k,1), data(:,2) - model.mus(k,2)];
%     Ps(:,k) = exp(-0.5 * sum((data_diff * inv(model.Sigmas(:,:,k))) .* data_diff,2)) * (model.pis(k) ./  ((2*pi).^(size(data,2)/2) * sqrt(det(model.Sigmas(:,:,k)))));
%     Ps(:,k) = mvnpdf(data, model.mus(k,:), model.Sigmas(:,:,k)) * model.pis(k);
    
%     Ps_vis(:,:,k) = reshape(mvnpdf(vis_XY, model.mus(k,:), model.Sigmas(:,:,k)) * model.pis(k), size(Ps_vis(:,:,k)));
  end
%   P_vis = sum(Ps_vis,3);
  
  Ps = max(eps, Ps);
  
  Ps_sum = sum(Ps,2);
  data_tau = Ps ./ repmat(Ps_sum, 1, size(Ps,2));
  
  score = sum(log(Ps_sum));
%   score = sum(sum(Ps .* repmat(data_N, 1, K),1),2);
%   score = sum(log(Ps_sum));
  if DISPLAY >= 0
    fprintf('%3d) \t LL = %f \t dLL = %f \t %s\n', iter, score, (score - score_bak), time2human(etime(clock, start_time)));
  end

  
  % M-step

  model.pis = sum(data_tau) ./ size(data,1);
%   model.pis = sum(data_tau .* repmat(data_N, 1, K), 1);

  for k = 1:K

    mult = data_tau(:,k);
    sum_mult = sum(mult,1);
    model.mus(k,:) = sum(repmat(mult, 1, size(data,2)) .* data,1) ./ sum_mult;
    
  end
  
  % Fix the mus to be the same
%   model.mus = repmat(sum(model.mus .* repmat(model.pis', 1, size(data,2)),1), size(model.mus,1), 1);

  for k = 1:K
    mult = data_tau(:,k);
    sum_mult = sum(mult,1);

    data_off = (data - repmat(model.mus(k,:), size(data,1), 1));
    model.Sigmas(:,:,k) = ((repmat(mult, 1, size(data,2)) .* data_off)' * data_off) ./ (sum(mult));
    model.Sigmas(:,:,k) = (model.Sigmas(:,:,k) + model.Sigmas(:,:,k)')/2;
    
    model.Sigmas(:,:,k) = model.Sigmas(:,:,k) + 10^-6 * eye(size(model.Sigmas(:,:,k)));
    
%     mult = 10^-9;
%     while det(model.Sigmas(:,:,k)) <= 10^-9
%       model.Sigmas(:,:,k) = model.Sigmas(:,:,k) + mult*eye(size(model.Sigmas(:,:,k)));
%       mult = mult*2;
%     end
    
    
%     if any(any(imag(model.Sigmas(:,:,k))))
%       keyboard
%     end
    
%     mult = data_tau(:,k);
%     sum_mult = sum(mult,1);
%     model.mus(k,:) = sum(repmat(mult, 1, 2) .* data,1) ./ sum_mult;
%     data_off = (data - repmat(model.mus(k,:), size(data,1), 1));
%     model.Sigmas(:,:,k) = ((repmat(mult, 1, 2) .* data_off)' * data_off) ./ (sum(mult));
  end
%   model.Sigmas = max(model.Sigmas, repmat(MIN_COV*eye(size(data,2),size(data,2)), [1, 1, K]));

%   model.mus
    
  if abs(((score - score_bak) ./ score)) < LOSS_FRACTION_THRESH
    DO_BREAK = true;
  end
  
  if iter < MIN_ITERS
    DO_BREAK = false;
  end
  
  score_bak = score;

%   if nargin < 3 || DISPLAY || DO_BREAK || (iter == N_ITERS)
%     
%     clf;
%     
%     range = [min(log(vis_N(vis_N>0))), max(log(vis_N(:)))];
%     subplot(2,4,1);
%     imagesc(log(vis_N), range);
%     title('data');
%     axis image off;
%     
%     subplot(2,4,2);
%     imagesc(log(P_vis*data_pixel_area), range);
%     title('model');
%     axis image off;
%     colormap('hot');
% 
%     subplot(2,4,3);
%     bar(sum(P_vis*data_pixel_area,1));
%     hold on;
%     plot(sum(vis_N,1), 'r-');
%     axis tight off;
% 
%     subplot(2,4,4);
%     bar(sum(P_vis*data_pixel_area,2));
%     hold on;
%     plot(sum(vis_N,2), 'r-');
%     axis tight off;
% 
%     
%     for k = 1:K
%       subplot(2,K,K+k)
% %       imagesc(log(Ps_vis(:,:,k)*data_pixel_area), range);
%       imagesc(Ps_vis(:,:,k));
%       title(['mixture ', num2str(k)]);
%       axis image off;
%     end
%     drawnow;
%     
%   end
%   

  if DO_BREAK
    break
  end

end  


LL = [];
for i = 1:size(model.mus,1)
  LL(i) = MOG_mvn_pdf(model.mus(i,:), model);
end
model.LL_max = max(LL);

% edges = {[0 : 3/401 : 3],[-1.5 : 3/401 : 1.5]};
% [jj,ii] = meshgrid(edges{2}, edges{1});
% 
% model.F_prob = -reshape(MOG_mvn_pdf(model, [ii(:), jj(:)]), size(ii));
% model.F_cost = model.F_prob - min(model.F_prob(:));
% 
% imagesc(mod(model.F_cost, 10)); colormap('hsv'); drawnow;
% 
% edges1 = edges{1};
% edges2 = edges{2};
% 
% edges1_start = edges1(1);
% edges1_step = edges1(2) - edges1(1);
% edges1_end = edges1(end);
% 
% edges2_start = edges2(1);
% edges2_step = edges2(2) - edges2(1);
% edges2_end = edges2(end);
% 
% model.edges = {edges1_start, edges1_step, edges1_end, edges2_start, edges2_step, edges2_end};




