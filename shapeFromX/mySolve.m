function [state] = mySolve(data, params, Z_init, L_init)

addpath(genpath('./minFunc_2012'));

sz = [size(data.true.im,1), size(data.true.im,2)];

S = min(params.MAX_PYR_DEPTH, 1 + maxPyrHt(sz, size(namedFilter(PYR_FILTER_NAME),1)));

if nargin < 3
  Z_init = zeros(sz);
end
Zpyr_init = buildLpyr_simple(Z_init, S, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
state.Zpyr = Zpyr_init;


if nargin < 4
  L_init = reshape(data.L_gaussian.mu, [9,3]);
end
% if ~params.SOLVE_LIGHT  
%   L_init = data.true.light;    
% end

if params.WHITEN_LIGHT
  state.L_white = ((reshape(L_init, [], 1)' - data.L_whiten_params.mean) * data.L_whiten_params.map);
else
  state.L_white = reshape(L_init, [], 1)';
end

data.L_precond = params.L_MULT{1};
state.L_white = data.L_precond*state.L_white;


start_time = clock;

params.FIX_LIGHT = ~params.SOLVE_LIGHT;
params.FIX_SHAPE = 0;

for ii = 1:params.N_SUB_ITERS 
  OPTIONS = struct('Method', 'lbfgs', ...
                   'MaxIter', params.N_ITERS_OPTIMIZE, ...
                   'Corr', params.LBFGS_NCORR, ...
                   'F_STREAK', params.F_STREAK, ...
                   'F_PERCENT', params.F_PERCENT, ...
                   'progTol', params.PROG_TOL, ...
                   'optTol', params.OPT_TOL, ...
                   'DISPLAY', params.MINFUNC_DISPLAY);
  state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
end

state_bak = state;
[loss, ~, state] = lossfun_sirfs(state_bak, data, params);

state.final_loss = loss;

fprintf(['Solved in ', time2human(etime(clock, start_time)), '\n']);

