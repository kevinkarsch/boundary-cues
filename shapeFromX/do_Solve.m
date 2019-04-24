function [state] = do_Solve(data, params, Z_init, L_init)

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
if ~params.SOLVE_LIGHT  
  L_init = data.true.light;    
end

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


  
if params.DEBUG_GRADIENT
  
  [vec, template] = cell2vector(state.Zpyr);
  vec = randn(size(vec));
  state.Zpyr = vector2cell(vec, template);
  
  state.L_white = randn(size(state.L_white));
  
%   if params.SOLVE_LIGHT
%     state.L_white(4,:) = state.L_white(1,:) + randn(size(state.L_white(1,:)))/10;
%     state.L_white(5,:) = state.L_white(2,:) + randn(size(state.L_white(1,:)))/100;
%   end
  
  %     data.temporary.Z = zeros(size(data.true.height));
  %     state = rmfield(state, 'Zdown');
  
  %     state.Lwhite = randn(size(state.Lwhite));
  %     state.Zdown = 10*randn(size(state.Zdown));
  %     data.L_active = data.L_meta.depth' <= params.LIGHT_DEPTH_START;
  [dy,dh] = checkgrad(state, params.LOSSFUN, 10^-5, data, params);
  %     keyboard
end

for ii = 1:params.N_SUB_ITERS
%   state = minimize_conjugate(state, params.LOSSFUN, 100, data, params);
%   state = minimize_simple(state, params.LOSSFUN, 400, data, params);

%   OPTIONS = struct('Method', 'lbfgs', 'MaxIter', 100, 'Corr', 1, 'F_STREAK', params.F_STREAK, 'F_PERCENT', params.F_PERCENT, 'progTol', params.PROG_TOL, 'optTol', params.OPT_TOL);
%   state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
  
  OPTIONS = struct('Method', 'lbfgs', 'MaxIter', params.N_ITERS_OPTIMIZE, 'Corr', params.LBFGS_NCORR, 'F_STREAK', params.F_STREAK, 'F_PERCENT', params.F_PERCENT, 'progTol', params.PROG_TOL, 'optTol', params.OPT_TOL);
  state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
  
%   state_bak = state;
%   
%   state = state_bak;
%   Z0 = reconLpyr_simple(state.Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
%   
%   Z = reconLpyr_simple(state.Zpyr, namedFilter(PYR_FILTER_NAME), 'zero');
%   state.Zpyr = buildLpyr_simple(Z, S, namedFilter(PYR_FILTER_NAME), 'zero');
%   Z1 = reconLpyr_simple(state.Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% %   imagesc([Z0, Z1]); figure; showLpyr_simple(state.Zpyr)  
%   
%   state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
  
end

state_bak = state;
[loss, junk, state] = lossfun_sirfs(state_bak, data, params);

state.final_loss = loss;

fprintf(['Solved in ', time2human(etime(clock, start_time)), '\n']);












% function [state] = do_Solve(data, params, Z_init, L_init)
% 
% addpath(genpath('./minFunc_2012'));
% 
% sz = [size(data.true.im,1), size(data.true.im,2)];
% 
% S = min(params.MAX_PYR_DEPTH, 1 + maxPyrHt(sz, size(namedFilter(PYR_FILTER_NAME),1)));
% 
% if nargin < 3
%   Z_init = zeros(sz);
% end
% 
% if nargin < 4
%   L_init = reshape(data.L_gaussian.mu, [9,3]);
% end
% 
% % Z_init = simpleShape_poisson(data.valid);
% 
% Zpyr_init = buildLpyr_simple(Z_init, S, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% % Zpyr_init = buildOGpyr(Z_init, S, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% state.Zpyr = Zpyr_init;
% 
% data.L_precond = params.L_MULT{1};
% 
% N_LIGHTS = params.N_LIGHTS{1};
% 
% 
% if ~params.SOLVE_LIGHT
%   
%   params.N_LIGHTS = {1};
%   
%   L = data.true.light;  
%   state.Ls_white = ((reshape(L, [], 1)' - data.L_whiten_params.mean) * data.L_whiten_params.map);
%   
% else
%   
%   if params.N_LIGHTS{1} == 1
%     
%     state.Ls_white = {};
%     L = L_init;
%     state.Ls_white{1} = ((reshape(L, [], 1)' - data.L_whiten_params.mean) * data.L_whiten_params.map);
%     state.Ls_white = cat(1,state.Ls_white{:});
%     
%   else
%     
%     %   li_map = [3,2,4,5,6:100];
%     li_map = [5,2,3,6:100];
%     
%     state.Ls_white = {};
%     for li = 1:N_LIGHTS
%       L = zeros(9,3);
%       L(1,:) = mean(L_init(1,:));
%       oi = li_map(li);
%       %   oi = li+1;
%       %   oi(oi>=4) = oi(oi>=4) + 1;
%       %   if oi >= 4
%       %     oi = oi + 1;
%       %   end
%       dim = floor((oi - 2)/2) + 2;
%       si = (mod(oi, 2)*2 - 1);
%       
%       L(dim,:) = si;
%       
%       L(3,:) = mean(L_init(3,:));
%       L(7,:) = mean(L_init(7,:));
%       %   figure; visSH_color(L)
%       
%       
%       
%       %     figure(li);
%       %     visSH_color(L);
%       %     drawnow;
%       
%       state.Ls_white{li} = ((reshape(L, [], 1)' - data.L_whiten_params.mean) * data.L_whiten_params.map);
%     end
%     state.Ls_white = cat(1,state.Ls_white{:});
%   end
% 
% end
%   
% state.Ls_white = data.L_precond*state.Ls_white;
% 
% start_time = clock;
% 
% params.FIX_LIGHT = ~params.SOLVE_LIGHT;
% params.FIX_SHAPE = 0;
% 
% 
%   
% if params.DEBUG_GRADIENT
%   
%   [vec, template] = cell2vector(state.Zpyr);
%   vec = randn(size(vec));
%   state.Zpyr = vector2cell(vec, template);
%   
%   state.Ls_white = randn(size(state.Ls_white));
%   
% %   if params.SOLVE_LIGHT
% %     state.Ls_white(4,:) = state.Ls_white(1,:) + randn(size(state.Ls_white(1,:)))/10;
% %     state.Ls_white(5,:) = state.Ls_white(2,:) + randn(size(state.Ls_white(1,:)))/100;
% %   end
%   
%   %     data.temporary.Z = zeros(size(data.true.height));
%   %     state = rmfield(state, 'Zdown');
%   
%   %     state.Lwhite = randn(size(state.Lwhite));
%   %     state.Zdown = 10*randn(size(state.Zdown));
%   %     data.L_active = data.L_meta.depth' <= params.LIGHT_DEPTH_START;
%   [dy,dh] = checkgrad(state, params.LOSSFUN, 10^-5, data, params);
%   %     keyboard
% end
% 
% N_SUB_ITERS = 1;
% for ii = 1:N_SUB_ITERS
% %   state = minimize_conjugate(state, params.LOSSFUN, 100, data, params);
% %   state = minimize_simple(state, params.LOSSFUN, 400, data, params);
% 
% %   OPTIONS = struct('Method', 'lbfgs', 'MaxIter', 100, 'Corr', 1, 'F_STREAK', params.F_STREAK, 'F_PERCENT', params.F_PERCENT, 'progTol', params.PROG_TOL, 'optTol', params.OPT_TOL);
% %   state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
%   
%   OPTIONS = struct('Method', 'lbfgs', 'MaxIter', params.N_ITERS_OPTIMIZE, 'Corr', params.LBFGS_NCORR, 'F_STREAK', params.F_STREAK, 'F_PERCENT', params.F_PERCENT, 'progTol', params.PROG_TOL, 'optTol', params.OPT_TOL);
%   state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% end
% 
% % solve_light = params.SOLVE_LIGHT;
% % params.SOLVE_LIGHT = 0;
% % state = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% % params.SOLVE_LIGHT = solve_light;
% 
% % state.height = reconLpyr_simple(state.Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% % state = rmfield(state, 'Zpyr');
% 
% 
% state_bak = state;
% [loss, junk, Z, L, Ls] = lossfun_saifs(state_bak, data, params);
% 
% % mult_bak = params.ZPYR_MULT{1};
% % params.ZPYR_MULT{1} = 0.25;%1/sqrt(2);
% % Zpyr = buildLpyr_simple(Z, S, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% % Zrecon = reconLpyr_simple(Zpyr, params.ZPYR_MULT{1}*namedFilter(PYR_FILTER_NAME), 'zero');
% % 
% % state2 = state;
% % state2.Zpyr = Zpyr;
% % state2 = minFunc(params.LOSSFUN, state2, OPTIONS, data, params);
% % [loss, junk, Z2, L] = lossfun_srfs_color_light_pyr_multi(state2, data, params);
% 
% state = [];
% state.height = Z;
% state.light = L;
% state.lights = Ls;
% 
% % Lw = state.Lwhite / data.L_precond;
% % state.light = reshape(Lw * data.L_whiten_params.inverse + data.L_whiten_params.mean, 9, []);
% % if size(state.light,2) == 1
% %   state.light = repmat(state.light, [1,3]);
% % end
% % state = rmfield(state, 'Lwhite');
% 
% 
% state.final_loss = loss;
% 
% fprintf(['Solved in ', time2human(etime(clock, start_time)), '\n']);




% function [state] = do_Solve(data, params)
% 
% do_break = false;
% 
% sz = [size(data.true.im,1), size(data.true.im,2)];
% 
% S = min(params.MAX_PYR_DEPTH, 1 + maxPyrHt(sz, size(namedFilter(PYR_FILTER_NAME),1)));
% 
% if params.DEBUG_GRADIENT
%   S = 3;
% end
% 
% [state.Zdown, data.pind] = downsample(zeros(sz), S);
% 
% pyr_blurs = S : -1 : 1;
% 
% params.FIX_LIGHT = 1;
% params.FIX_SHAPE = 0;
% 
% start_time = clock;
% 
% if params.SOLVE_LIGHT
% 
%   if isinf(params.LIGHT_DEPTH_START)
%     data.L_active = true(1,size(data.true.light,3));
%   else
%     data.L_active = data.L_meta.depth' <= params.LIGHT_DEPTH_START;
%   end
%   
% end
% 
% for pyr_blur = pyr_blurs
% 
%   fprintf('Scale '); fprintf('%f ', pyr_blur); fprintf('\n');  
% 
% %   if params.SOLVE_LIGHT && ~params.USE_GREEDY
% %     if params.DO_DISPLAY
% %       X = data.L_meta.X;
% %       Y = data.L_meta.Y;
% %       children = data.L_meta.children;
% %       V = data.L_meta.vis;
% %       keep = data.L_active;
% %       
% %       children = children(keep, keep);
% %       X = X(keep);
% %       Y = Y(keep);
% %       V = V(keep);
% %       
% %       [i,j] = find(children);
% %       
% %       figure(4);
% %       line(X([i,j])', Y([i,j])', 'color', 'k')
% %       axis equal off
% %       
% %       hold on;
% %       for vi = 1:length(V)
% %         imagesc(X(vi) - size(V{vi},1)/2, Y(vi) - size(V{vi},1)/2, V{vi});
% %       end
% %       hold off;
% %       
% %     end
% %   end
%     
%   N_SUB_ITERS = 3;
%   OPTIONS = struct('Method', 'lbfgs', 'MaxIter', 1000, 'Corr', params.LBFGS_NCORR, 'TolX', 0, 'TolFun', 0, 'F_STREAK', 5, 'F_PERCENT', params.F_PERCENT, 'F_MINITERS', 10);
%   
%   if params.DEBUG_GRADIENT
%     state.Zdown = 1*randn(size(state.Zdown));
%     checkgrad(state, params.LOSSFUN, 10^-5, data, params);
%     state = minimize_conjugate(state, params.LOSSFUN, 3, data, params);
%   end
%   
%   for ii = 1:N_SUB_ITERS
%     [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
%     if (ii < N_SUB_ITERS) && ((length(output.trace.fval) <= 5) || (exitflag == 2))
% %       state = minimize_simple(state, params.LOSSFUN, 10, data, params);
%       [state, ftrace] = minimize_conjugate(state, params.LOSSFUN, 3, data, params);
%     else
%       break;
%     end
% %     elseif (exitflag == 2)
% %       [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% %     end
%   end
%   
%   if pyr_blur ~= 1
% %     state.Zdown = imresize(state.Zdown, data.pind(end-1,:), 'bicubic');
%     state.Zdown = upsample(state.Zdown, data.pind(end-1:end,:));
%     data.pind = data.pind(1:end-1,:);
%   end
%   
%   if params.SOLVE_LIGHT
%     
%     global global_P global_loss_best
%     
%     if pyr_blur <= params.GREEDY_STOP_DEPTH
%       LL = global_P;
%       data.L_active = find(LL >= (max(LL) - params.GREEDY_LL_MIN));
%       params.USE_GREEDY = 0;
%     end
%     
%     
% %     if ~params.USE_GREEDY
% %       for prune_i = 1:10
% %         
% %         P = global_P;
% %         
% %         p_uniform = 1/sum(data.L_active);
% %         
% %         children = data.L_meta.children;
% %         
% %         do_split = P > (params.LIGHT_RELPROB_MAX * p_uniform); %params.LIGHT_P_MAX{1};
% %         
% %         do_prune = data.L_active & (P < params.LIGHT_PROB_MIN);%(P < params.LIGHT_P_MIN{1});
% %         do_prune(any(children(data.L_active,:),1)) = false;
% %         
% %         do_add = false(size(data.L_active));
% %         do_add(any(data.L_meta.children(:,do_split),2)) = true;
% %         do_add(data.L_active) = false;
% %         do_add(do_prune) = false;
% %         
% %         if any(do_prune)
% %           fprintf('Pruning ')
% %           fprintf('%d ', find(do_prune));
% %           fprintf('\n');
% %           data.L_active(do_prune) = false;
% %         end
% %         
% %         if any(do_add)
% %           fprintf('Adding ')
% %           fprintf('%d ', find(do_add));
% %           fprintf('\n');
% %           data.L_active(do_add) = true;
% %         end
% %         
% %         if any(do_add) || any(do_prune)
% %           global_loss_best = inf;
% %           feval(params.LOSSFUN, state, data, params);
% %         else
% %           fprintf('Nothing to add or prune\n');
% %           break;
% %         end
% %         
% %         
% %       end
% %       
% %       %     global global_P;
% %       %     data.L_active(global_P < params.LIGHT_P_CUTOFF) = false;
% %       
% %       fprintf('%d/%d Lights Active\n', sum(data.L_active), length(data.L_active))
% %       fprintf('P(lights): ');
% %       fprintf('%0.3f ', sort(global_P(data.L_active), 'descend'));
% %       fprintf('\n');
% %       
% %       %     if (pyr_blur == pyr_blurs(1)) & ~isinf(params.LIGHT_K_START)
% %       %       global global_loss_best
% %       %       global_loss_best = inf;
% %       %       data.L_active = true(size(data.L_active));
% %       %       feval(params.LOSSFUN, state, data, params);
% %       %     end
% %       %
% %       %     global global_P;
% %       %
% %       %     data.L_active(global_P < params.LIGHT_P_CUTOFF) = false;
% %       %     fprintf('%d/%d Lights Active\n', sum(data.L_active), length(data.L_active))
% %       %     fprintf('P(lights): ');
% %       %     fprintf('%0.3f ', sort(global_P(data.L_active), 'descend'));
% %       %     fprintf('\n');
% %       
% %     end
% 
%   end
%   
% end
% 
% if params.SOLVE_LIGHT
%   [junk, junk, L_expected, L_prob] =  lossfun_srfs_color_light(state, data, params);
%   state.light = L_expected;
%   state.light_prob = L_prob;
% end
% 
% state.height = state.Zdown;
% state = rmfield(state, 'Zdown');
% 
% state.final_loss = output.trace.fval(end);
% 
% 
% fprintf(['Solved in ', time2human(etime(clock, start_time)), '\n']);














% function [state] = do_Solve(data, params)
% 
% do_break = false;
% 
% sz = [size(data.true.im,1), size(data.true.im,2)];
% 
% S = min(params.MAX_PYR_DEPTH, 1 + maxPyrHt(sz, 5));
% 
% if params.DEBUG_GRADIENT
%   S = 3;
% end
% 
% [state.Zdown, data.pind] = downsample(zeros(sz), S);
% 
% pyr_blurs = S : -1 : 1;
% 
% params.FIX_LIGHT = 1;
% params.FIX_SHAPE = 0;
% 
% start_time = clock;
% 
% if params.SOLVE_LIGHT
%   
% %   if isinf(params.LIGHT_K_START)
%   data.L_active = true(1,size(data.true.light,3));
% %   else
% %     data.L_active = false(1,size(data.true.light,3));
% %     data.L_active(data.L_meta.clusters{params.LIGHT_K_START}) = true;
% %   end
% end
% 
% for pyr_blur = pyr_blurs
% 
%   fprintf('Pyr Blur: '); fprintf('%f ', pyr_blur); fprintf('\n');  
%   
%   N_SUB_ITERS = 3;
%   OPTIONS = struct('Method', 'lbfgs', 'MaxIter', params.N_ITERS_OPTIMIZE, 'Corr', params.LBFGS_NCORR, 'TolX', 0, 'TolFun', 0, 'F_STREAK', 5, 'F_PERCENT', params.F_PERCENT, 'F_MINITERS', 10);
%   
%   if params.DEBUG_GRADIENT
%     state.Zdown = 1*randn(size(state.Zdown));
%     checkgrad(state, params.LOSSFUN, 10^-5, data, params);
%     state = minimize_conjugate(state, params.LOSSFUN, 3, data, params);
%   end
%   
%   for ii = 1:N_SUB_ITERS
%     [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
%     if (ii < N_SUB_ITERS) && ((length(output.trace.fval) <= 5) || (exitflag == 2))
%       [state, ftrace] = minimize_conjugate(state, params.LOSSFUN, 3, data, params);
%     else
%       break;
%     end
% %     elseif (exitflag == 2)
% %       [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% %     end
%   end
%   
%   if pyr_blur ~= 1
%     state.Zdown = upsample(state.Zdown, data.pind(end-1:end,:));
%     data.pind = data.pind(1:end-1,:);
%   end
%   
%   if params.SOLVE_LIGHT
%     
% %     global global_P
% %     
% %     
% %     for prune_i = 1:10
% %       
% %       P = global_P;
% %       
% %       Dsq = data.L_meta.Dsq;
% %       params.LIGHT_DISTANCE_GAMMA = 0.01;
% %       params.LIGHT_P_MIN = 0.01;
% %       
% %       A = exp(-params.LIGHT_DISTANCE_GAMMA*Dsq);
% %       A = A ./ repmat(sum(A,1), size(A,1), 1);
% %       
% %       P = (A*P')';
% % 
% %       data.L_active = P > params.LIGHT_P_MIN;
% %       feval(params.LOSSFUN, state, data, params);
% % 
% %     end
% %       
% % %       data.true.light;
% % %       
% % %       figure; visLights(Ls, P)
% % %       figure; visLights(Ls, P2)
% % %       
% % %       p_uniform = 1/sum(data.L_active);
% % %       
% % %       do_split = P > (params.LIGHT_RELPROB_MAX * p_uniform); %params.LIGHT_P_MAX{1};
% % %       do_prune = data.L_active & (P < params.LIGHT_RELPROB_MIN * p_uniform);%(P < params.LIGHT_P_MIN{1});
% % %       
% % %       if any(do_prune)
% % %         fprintf('Pruning ')
% % %         fprintf('%d ', find(do_prune));
% % %         fprintf('\n');
% % %         
% % %         data.L_active(do_prune) = false;
% % %       
% % %       end
% % %       
% % %       do_add = false(size(data.L_active));
% % %       do_add(data.L_meta.neighbors(1:params.LIGHT_K_SPLIT,do_split)) = true;
% % %       do_add(data.L_active) = false;
% % %       data.L_active(do_add) = true;
% % % 
% % %       if any(do_add)
% % %         fprintf('Adding ')
% % %         fprintf('%d ', find(do_add));
% % %         fprintf('\n');
% % % 
% % %         feval(params.LOSSFUN, state, data, params);
% % % 
% % %       else
% % %         fprintf('Nothing to add\n');
% % %         break;
% % %       end
% % %       
% % %     end
% %     
% %     
% % %     for prune_i = 1:10
% % %       
% % %       L_active_bak = data.L_active;
% % %       
% % %       alpha = params.ALPHA_LIGHT{1};
% % %       LL = -alpha * losses_A;
% % %       LL = LL - max(LL);
% % %       P = exp(LL);
% % %       P = P ./ sum(P(~isnan(P)));
% % %       P(isnan(P)) = 0;
% % %       
% % %       p_uniform = 1/sum(data.L_active);
% % %       
% % %       do_split = P > (params.LIGHT_RELPROB_MAX * p_uniform); %params.LIGHT_P_MAX{1};
% % %       do_prune = data.L_active & (params.LIGHT_RELPROB_MIN * p_uniform);%(P < params.LIGHT_P_MIN{1});
% % %       
% % %       if any(do_prune)
% % %         fprintf('Pruning ')
% % %         fprintf('%d ', find(do_prune));
% % %         fprintf('\n');
% % %         losses_A(do_prune) = nan;
% % %         data.L_active(do_prune) = false;
% % %       end
% % %       
% % %       
% % %       do_add = false(size(data.L_active));
% % %       do_add(data.L_meta.neighbors(1:params.LIGHT_K_SPLIT,do_split)) = true;
% % %       do_add(L_active_bak) = false;
% % %       
% % %       if any(do_add)
% % %         fprintf('Adding ')
% % %         fprintf('%d ', find(do_add));
% % %         fprintf('\n');
% % %         L_active_bak2 = data.L_active;
% % %         data.L_active = do_add;
% % %         feval(params.LOSSFUN, state, data, params);
% % %         losses_A(do_add) = global_losses_A(do_add);
% % %         data.L_active = L_active_bak2;
% % %         data.L_active(do_add) = true;
% % %       else
% % %         fprintf('Nothing to add!\n');
% % %         global_P = P;
% % %         break
% % %       end
% % %       
% % %       
% % %       
% % %     end
%     
%     
%     fprintf('%d/%d Lights Active\n', sum(data.L_active), length(data.L_active))
%     
%     [P_sorted, P_sortidx] = sort(global_P(data.L_active), 'descend');
%     fprintf('Lights: ');
%     fprintf('%d\t', P_sortidx);
%     fprintf('\n');
%     fprintf('P:      ');
%     fprintf('%0.3f\t', P_sorted);
%     fprintf('\n');
% 
%   end
%   
% end
% 
% if params.SOLVE_LIGHT
%   [junk, junk, L_expected, L_prob] =  lossfun_srfs_color_light(state, data, params);
%   state.light = L_expected;
%   state.light_prob = L_prob;
% end
% 
% state.height = state.Zdown;
% state = rmfield(state, 'Zdown');
% 
% state.final_loss = output.trace.fval(end);
% 
% 
% fprintf(['Solved in ', time2human(etime(clock, start_time)), '\n']);
% 
% 
% 
% % function [state] = do_Solve(state, data, params)
% % 
% % do_break = false;
% % 
% % S = 'auto';
% % if params.CHEAT
% %   S = params.CHEAT;
% % end
% % 
% % keyboard
% % [pyr, data.pind] = buildGpyr(state.height, S);
% % 
% % S = size(data.pind,1);
% % state = rmfield(state, 'height');
% % state.Zdown = pyrBand(pyr, data.pind, S);
% % 
% % pyr_blurs = fliplr([1 : params.BLUR_STEP : S]);
% % 
% % if pyr_blurs(end) ~= 1
% %   pyr_blurs = [pyr_blurs, 1];
% % end
% % 
% % pyr_blurs
% % 
% % TolX = 0;%10^-9;
% % TolFun = 0;%10^-9;
% % 
% % params.FIX_LIGHT = 1;
% % params.FIX_SHAPE = 0;
% % 
% % start_time = clock;
% % 
% % 
% % for pyr_blur = pyr_blurs
% % 
% %   fprintf('Pyr Blur: '); fprintf('%f ', pyr_blur); fprintf('\n');
% %   data.pyramid_blur_level = pyr_blur;
% %   
% %   if(~all(data.pind(floor(pyr_blur),:) == size(state.Zdown)))
% %     
% %     Zpind = data.pind;
% %     
% %     Zpind_top = Zpind(end-1:end,:);
% %     Zpyr = zeros(sum(prod(Zpind_top,2),1),1);
% %     Zpyr(pyrBandIndices(Zpind_top, 2)) = state.Zdown(:);
% %     state.Zdown = reconLpyr(Zpyr, Zpind_top, 'all', 'binom5', 'zero');
% %     
% %     data.pind = Zpind(1:floor(pyr_blur),:);
% %     
% %   end
% %   
% % %   state.Zdown = randn(size(state.Zdown));
% % %   [dy,dh] = checkgrad(state, 'lossfun_srfs', 1e-5, data, params);
% % %   keyboard
% %   
% %   N_SUB_ITERS = 1;
% %   OPTIONS = struct('Method', 'lbfgs', 'MaxIter', 1000, 'Corr', params.LBFGS_NCORR, 'F_CONSTRAIN', params.F_CONSTRAIN, 'TolX', TolX, 'TolFun', TolFun, 'F_STREAK', 5, 'F_PERCENT', params.F_PERCENT, 'F_MINITERS', 10);
% %   
% %   for ii = 1:N_SUB_ITERS
% %     [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% %     if (length(output.trace.fval) <= 5)
% %       state = minimize_conjugate(state, params.LOSSFUN, 3, data, params);
% %       [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% %     elseif (exitflag == 2)
% %       [state, junk1, exitflag, output] = minFunc(params.LOSSFUN, state, OPTIONS, data, params);
% %     end
% %   end
% %   
% % end
% % 
% % if params.SOLVE_LIGHT
% %   [junk, junk, L_expected, L_prob] =  lossfun_srfs_light_fast(state, data, params);
% %   state.light = L_expected;
% %   state.light_prob = L_prob;
% % end
% % 
% % state.height = state.Zdown;
% % state = rmfield(state, 'Zdown');
% % 
% % state.final_loss = output.trace.fval(end);
% % 
% % 
% % fprintf(['Solved in ', time2human(etime(clock, start_time)), '\n']);
