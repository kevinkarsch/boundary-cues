% function [] = interp1_fixed_construct(range, n_bins, fname, dfname)

range = [-11.521, 10.12];
n_bins = 50;
% fname = @sin;
% dfname = @cos;

bin_width = (range(2) - range(1))/(n_bins);
low_val = range(1);

X = [range(1):bin_width:range(2)]';

% LUT_F = feval(fname, X);
% LUT_dF = feval(dfname, X);

% LUT_F = (X-3).^2;
% LUT_dF = 2*(X-3);

LUT_F = sin(X);
LUT_dF = cos(X);

X0 = [range(1)-5:bin_width/5:range(2)+5]';
% X0 = [range(1):bin_width/5:range(2)]';
[F0, dF0] = interp1_fixed(X0, low_val, bin_width, LUT_F, LUT_dF);

% plot(X0, (X0-3).^2, 'k--'); hold on;
% plot(X0, 2*(X0-3), 'k--'); hold on;
plot(X, LUT_F, 'b-'); hold on;
scatter(X0, F0, 'bo'); hold on;

plot(X, LUT_dF, 'r-'); hold on;
scatter(X0, dF0, 'ro'); hold on;
