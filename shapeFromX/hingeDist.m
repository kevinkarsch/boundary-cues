function [loss, d_loss_X, d_loss_Y] = hingeDist(X, Y, delta)
% loss = \sum_i,j max(0, Y_i - X_j - delta)

X = X + delta; % We can account for delta by just shifting X (or Y) and forgetting about it

% Sort X and Y, and integrate X
[Xs, X_sortidx] = sort(X, 'ascend');
[Ys, Y_sortidx] = sort(Y, 'ascend');
Xc = cumsum(Xs, 1);

% Concatenate and sort the sorted X and Y. This tells us how many entries
% of Xs are below each entry in Ys, and vice versa.
XY = [Xs; Ys];
[junk, XY_sortidx] = sort(XY, 'ascend');
tf = XY_sortidx <= length(X);

c = cumsum(tf);
js = c(~tf); % The count of how many Xs are below each Y

c = cumsum(~tf);
ks = c(tf); % The count of how many Ys are below each X

loss = sum(Ys .* js) - sum(Xc(js(js>0)));

d_loss_X = zeros(size(X));
d_loss_X(X_sortidx) = ks - length(Y);

d_loss_Y = zeros(size(Y));
d_loss_Y(Y_sortidx) = js;


n = 1/(length(X) * length(Y));

loss = loss * n;
d_loss_X = d_loss_X * n;
d_loss_Y = d_loss_Y * n;



% function [loss, d_loss_X, d_loss_Y] = hingeDist(X, Y, delta)
% D = bsxfun(@minus, Y, X');
% hinge = max(0, D-delta);
% 
% loss = sum(hinge(:)) /numel(hinge);
%   
% d_loss_X = -1/numel(hinge) * sum(hinge > 0,1)';
% d_loss_Y =  1/numel(hinge) * sum(hinge > 0,2);




% function [loss, d_loss] = hingeDist(XY, delta)
% X = XY{1};
% Y = XY{2};
% % function [loss, d_loss_X, d_loss_Y] = hingeDist(X, Y, delta)
% 
% X = X + delta;
% 
% % D = bsxfun(@minus, Y, X');
% % hinge = max(0, D);
% % % hinge = max(0, D-delta);
% % loss0 = sum(hinge(:));
% % 
% % d_loss_X = -sum(hinge > 0,1)';
% % d_loss_Y =  sum(hinge > 0,2);
% % d_loss = {d_loss_X, d_loss_Y};
% 
% % return
% 
% % [Xs, X_sortidx] = sort(X, 'ascend');
% % [Ys, Y_sortidx] = sort(Y, 'ascend');
% % Xc = cumsum(Xs, 1);
% % % Xc = flipud(cumsum(flipud(Xs)));
% % 
% % js = [];
% % loss1 = 0;
% % for i = 1:length(Ys)
% % %   loss = loss + sum(max(0, Ys(i) - Xs));
% %   
% % %   sum(max(0, Y(i) - X))
% % %   
% % %   sum(Y(i) - Xs(1:j))
% % %   j*Y(i) - sum(Xs(1:j))
% % 
% %   j = find(Ys(i) >= Xs, 1, 'last');
% %   loss1 = loss1 + j*Ys(i) - Xc(j);
% %   
% %   js(i) = j;
% %   
% % end
% % 
% % XY = [Xs; Ys];
% % [XYs, XY_sortidx] = sort(XY, 'ascend');
% % tf = XY_sortidx <= length(X);
% % % tf = ismember(XYs, Xs);
% % c = cumsum(tf);
% % js = c(~tf);
% % 
% % loss = sum(Ys .* js) - sum(Xc(js));
% 
% 
% [Xs, X_sortidx] = sort(X, 'ascend');
% [Ys, Y_sortidx] = sort(Y, 'ascend');
% Xc = cumsum(Xs, 1);
% 
% XY = [Xs; Ys];
% [XYs, XY_sortidx] = sort(XY, 'ascend');
% tf = XY_sortidx <= length(X);
% c = cumsum(tf);
% js = c(~tf);
% 
% YX = [Ys; Xs];
% [YXs, YX_sortidx] = sort(YX, 'ascend');
% tf = YX_sortidx <= length(Y);
% c = cumsum(tf);
% ks = c(~tf);
% 
% loss = sum(Ys .* js) - sum(Xc(js(js>0)));
% 
% % [loss0, loss]
% % loss - loss0
% 
% d_loss_X = zeros(size(X));
% d_loss_X(X_sortidx) = ks - length(Y);
% 
% d_loss_Y = zeros(size(Y));
% d_loss_Y(Y_sortidx) = js;
% 
% d_loss = {d_loss_X, d_loss_Y};



% D = bsxfun(@minus, Y, X');
% hinge = max(0, D-delta);
% 
% loss = sum(hinge(:)) /numel(hinge);
%   
% d_loss_X = -1/numel(hinge) * sum(hinge > 0,1)';
% d_loss_Y =  1/numel(hinge) * sum(hinge > 0,2);


% X = X + delta;
% 
% D = bsxfun(@minus, Y, X');
% hinge = max(0, D);
% 
% loss = sum(hinge(:)) /numel(hinge);
%   
% d_loss_X = -1/numel(hinge) * sum(hinge > 0,1)';
% d_loss_Y =  1/numel(hinge) * sum(hinge > 0,2);



% X = X + delta;
% 
% [Xs, X_sortidx] = sort(X, 'ascend');
% [Ys, Y_sortidx] = sort(Y, 'ascend');
% 
% keyboard
% 
% Xc = cumsum(Xs);
% Yc = cumsum(Ys);
% 
% XYs = [Xs;Ys];
% XYc = [Xc;Yc];
% [junk, XYs_sortidx] = sort(XYs, 'ascend');
% 
% 
% D = bsxfun(@minus, Ys, Xs');
% hinge = max(0, D);
% sum(hinge(:))
% 
% loss = sum(hinge(:)) /numel(hinge);
%   
% d_loss_X = -1/numel(hinge) * sum(hinge > 0,1)';
% d_loss_Y =  1/numel(hinge) * sum(hinge > 0,2);
