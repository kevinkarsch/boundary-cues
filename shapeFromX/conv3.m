function Xc = conv3(X, f)

% Xp = [X(1,1), X(1,:), X(1,end); 
%       X(:,1), X,      X(:,end);
%       X(end,1), X(end,:), X(end,end)];
% 
% Xc = conv2(Xp, f, 'valid');


Xc = conv2(X, f, 'same');
Xc(end,:) = Xc(end,:) + conv2(X(end,:), f(1,:), 'same');
Xc(1,:)   = Xc(1,:)   + conv2(X(1,:),   f(end,:), 'same');
Xc(:,end) = Xc(:,end) + conv2(X(:,end), f(:,1), 'same');
Xc(:,1)   = Xc(:,1)   + conv2(X(:,1),   f(:,end), 'same');

Xc(1,1)     = Xc(1,1)     + X(1,1)     .* f(end,end);
Xc(1,end)   = Xc(1,end)   + X(1,end)   .* f(end,1);
Xc(end,1)   = Xc(end,1)   + X(end,1)   .* f(1,end);
Xc(end,end) = Xc(end,end) + X(end,end) .* f(1,1);


% fx = [1,2,1]/4;
% fy = [-1;0;1]/2;
% 
% tic; for i = 1:1000; Xc = conv2(Xp, f, 'valid'); end; toc
% tic; for i = 1:1000; Xc2 = conv2(fy, fx, Xp, 'valid'); end; toc

