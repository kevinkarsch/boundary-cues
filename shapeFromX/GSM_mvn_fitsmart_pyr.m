function model = GSM_mvn_fitsmart(X, K, FIT_MU)

addpath(genpath('./matlabPyrTools/'));

if nargin < 3
  FIT_MU = 0;
end

vars = 1000.^((-(K-1)/2 : (K-1)/2) / ((K-1)/2));

pis = ones(size(vars)) / numel(vars);

data = [];
data.X = X;
data.vars = vars;

params = [];
params.USE_NUMERICAL = 1;
params.F_STREAK = 5;
params.F_PERCENT = 0.01;

state = [];
params.PYR_FILT = 2*namedFilter('binom5');
state.log_pis_pyr = buildGpyr_simple(zeros(size(pis)), floor(log2(K)), params.PYR_FILT, 'zero');

if FIT_MU
  state.mu = mean(X,1);
  state.Sigma_vec = PSD2vec(cov(X));
else
  state.Sigma_vec = PSD2vec((X'*X) / size(X,1));
end

addpath(genpath('./minFunc_2012'));
N_ITERS = 5000;
OPTIONS = struct('Method', 'lbfgs', 'MaxIter', N_ITERS, 'MaxFunEvals', N_ITERS*3, 'Corr', 500, 'F_STREAK', 5, 'F_PERCENT', 0.0001, 'numDiff', 1, 'useComplex', 0);
state = minFunc('GSM_mvn_gradfit_pyr_loss', state, OPTIONS, data, params);

state.log_pis = reconLpyr_simple(state.log_pis_pyr, params.PYR_FILT, 'zero');

% state = minimize_conjugate(state, 'GSM_mvn_gradfit_loss', 100, data, params);

if FIT_MU
  model.mu = state.mu;
end
model.pis = exp(state.log_pis*20);
model.pis = model.pis ./ sum(model.pis);
model.Sigma = vec2PSD(state.Sigma_vec);
model.vars = data.vars;

% LL_GSM = sum(GSM_mvn_pdf(model, X));
% LL_gaussian = sum(lmvnpdf(X, [0,0,0], (X'*X) / size(X,1)));

model.Sigma_inv = inv(model.Sigma);
model.Sigma_R = cholcov(model.Sigma,0);
model.logmults = [];

log_len = log((2*pi).^(size(X,2)/2));
for k = 1:K
  R = model.Sigma_R * sqrt(model.vars(k));
  model.logmults(k) = log(model.pis(k)) - log_len - sum(log(diag(R)));
end


n_bins = 100000;

X_icov = X * model.Sigma_inv;
mahal_dist = 0.5 * sum(X_icov .* X,2);

bin_range = [0, max(mahal_dist)];
bin_width = (bin_range(2) - bin_range(1))/(n_bins-1);
bins = [bin_range(1):bin_width:bin_range(2)]';

P = 0;
p_div = 0;

for k = 1:K
  p = exp(model.logmults(k) + bins./-model.vars(k));
  P = P + p;
  p_div = p_div + p./-model.vars(k);
end

LL_bin = log(P);

model.lut.bin_range = bin_range;
model.lut.bin_width = bin_width;
model.lut.n_bins = n_bins;
model.lut.F = LL_bin;

model.LL_zero = GSM_mvn_pdf(model, zeros(1, size(data.X,2)), 1);

bins = model.lut.bin_range(1) : model.lut.bin_width : model.lut.bin_range(2);


% sqrt_mahal_dist = sqrt(mahal_dist);
% n_hist_bins = size(X,1)/10;
% be = sqrt(bins(end));
% hist_bins = 0:(be/n_hist_bins):be;
% N = hist(sqrt_mahal_dist, hist_bins);
% N = N ./ sum(N);
% 
% logN = nan(size(N));
% logN(N > 0) = log(N(N>0))

figure(1002); clf;
plot(sqrt(bins), LL_bin); hold on;
% plot(hist_bins, logN, 'r-');
xlabel('Mahalanobis Distance');
ylabel('Log-Likelihood');
axis square

