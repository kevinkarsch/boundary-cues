function Xc = conv3d(X, f)

Xc_full = conv2(X, f, 'full');

% Xc = Xc_full(2:end-1,2:end-1);
% 
% Xc(1,:) = Xc(1,:) + Xc_full(1,2:end-1);
% Xc(:,1) = Xc(:,1) + Xc_full(2:end-1,1);
% Xc(end,:) = Xc(end,:) + Xc_full(end,2:end-1);
% Xc(:,end) = Xc(:,end) + Xc_full(2:end-1,end);
% 
% Xc(1,1) = Xc(1,1) + Xc_full(1,1);
% Xc(1,end) = Xc(1,end) + Xc_full(1,end);
% Xc(end,1) = Xc(end,1) + Xc_full(end,1);
% Xc(end,end) = Xc(end,end) + Xc_full(end,end);

Xc = unpad1(Xc_full);
