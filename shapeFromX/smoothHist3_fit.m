function [F, LL] = smoothHist3_fit(RGB, bin_low, bin_high, lambda_smooth, robust_cost)

addpath(genpath('./minFunc_2012'));
addpath(genpath('./matlabPyrTools/'));

X_splat = splat3(RGB, bin_low, bin_high);
X_splat.N = X_splat.N ./ sum(X_splat.N(:));

% pyrFilt = sqrt(2)*namedFilter('binom5');
pyrFilt = sqrt(2)*namedFilter('binom5');

W = zeros(size(X_splat.N));
Wpyr = buildGpyr3_simple(W, 6, pyrFilt);

lambda_pos = 0;%1000000;

% checkgrad(Wpyr, 'lossfun_hist3', 10^-5, X_splat, lambda_smooth, lambda_pos, pyrFilt);

N_ITERS = 5000;
OPTIONS = struct('Method', 'lbfgs', 'MaxFunEvals', 10*N_ITERS, 'MaxIter', N_ITERS, 'Corr', 50, 'F_STREAK', 5, 'F_PERCENT', 0.001);
for ii = 1:2
  Wpyr = minFunc('lossfun_hist3', Wpyr, OPTIONS, X_splat, lambda_smooth, lambda_pos, robust_cost, pyrFilt);
end

W = reconLpyr3_simple(Wpyr, pyrFilt);

% W = convn(padarray(W, [1,1,1], 'replicate'), ones(2,2,2)/8);
% W = W(2:end-1, 2:end-1, 2:end-1);

P = exp(-W);
Z = X_splat.bin_area*sum(P(:));
P = P ./ Z;
F = -log(P);
LL = -F;

F = F - min(F(:));


% for i = 1:size(W,1)
%   surf(squeeze(W(i,:,:)), 'EdgeColor', 'none'); axis vis3d off; drawnow;
% end
% 
% for i = 1:size(W,2)
%   surf(squeeze(W(:,i,:)), 'EdgeColor', 'none'); axis vis3d off; drawnow;
% end
% 
% for i = 1:size(W,3)
%   surf(squeeze(W(:,:,i)), 'EdgeColor', 'none'); axis vis3d off; drawnow;
% end


% function X_splat = smoothHist3_fit(RGB, sigma, bin_low, bin_high, lambda_smooth)
% 
% addpath(genpath('./minFunc_2012'));
% addpath(genpath('./matlabPyrTools/'));
% 
% X_splat = splat3(RGB, sigma, bin_low, bin_high);
% 
% pyrFilt = sqrt(2)*namedFilter('binom5');
% 
% W = zeros(size(X_splat.N));
% Wpyr = buildGpyr3_simple(W, 6, pyrFilt);
% 
% lambda_pos = 0;%1000000;
% 
% N_ITERS = 5000;
% OPTIONS = struct('Method', 'lbfgs', 'MaxFunEvals', 10*N_ITERS, 'MaxIter', N_ITERS, 'Corr', 50, 'F_STREAK', 5, 'F_PERCENT', 0.001);
% for ii = 1:4
%   Wpyr = minFunc('lossfun_hist3', Wpyr, OPTIONS, X_splat, lambda_smooth, lambda_pos, pyrFilt);
% end
% 
% W = reconLpyr3_simple(Wpyr, pyrFilt);
% 
% P = exp(-W);
% Z = X_splat.bin_area*sum(P(:));
% P = P ./ Z;
% X_splat.W = -log(P);
% 
% 
% X_splat.F = X_splat.W - min(X_splat.W(:));
% 
% % for i = 1:size(W,1)
% %   surf(squeeze(W(i,:,:)), 'EdgeColor', 'none'); axis vis3d off; drawnow;
% % end
% % 
% % for i = 1:size(W,2)
% %   surf(squeeze(W(:,i,:)), 'EdgeColor', 'none'); axis vis3d off; drawnow;
% % end
% % 
% % for i = 1:size(W,3)
% %   surf(squeeze(W(:,:,i)), 'EdgeColor', 'none'); axis vis3d off; drawnow;
% % end
