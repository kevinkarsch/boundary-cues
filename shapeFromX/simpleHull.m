function [A,b] = simpleHull(X_raw, N_CONSTRAINTS)

[k, vol_before] = convhulln(X_raw);
X = X_raw(unique(k),:);

if N_CONSTRAINTS < size(X,1)
  [idx, X_thin] = kmeans_ann(X, N_CONSTRAINTS, 'N_RESTARTS', 50);
else
  X_thin = X;
end

X_thin = X_thin(unique(convhulln(X_thin)),:);

[A,b] = vert2con(X_thin);

mag = sqrt(sum(A.^2,2));
A = A ./ repmat(mag, 1, size(A,2));

if size(X_raw,2) == 2
  return;
end

D = A*X';
b = max(D, [], 2);

error = max(max(bsxfun(@minus, A*X', b)))

% assert(error < 10^-8);
% 
% X_recon = con2vert(A, b);
% [k, vol_after] = convhulln(X_recon);
% 
% vol_before
% vol_after
% 
% % assert(vol_after >= vol_before)
% 
% if size(X_recon,2) == 3
%   trisurf(k, X_recon(:,1), X_recon(:,2), X_recon(:,3), 'Facecolor','cyan'); axis equal;
%   drawnow;
% end





% D = bsxfun(@minus, A*X', b);
% max(D(:))
% 
% vol = nan(size(b));
% for remove_i = 1:size(A,1)
%   
%   keep = [1:(remove_i-1),(remove_i+1):length(b)];
%   A2 = A(keep,:);
%   b2 = b(keep);
%   X_recon = con2vert(A2,b2);
%   [junk v]=convhulln(unique(X_recon, 'rows'));
%   vol(remove_i) = v;
%   
% end
% 
% [junk, sortidx] = sort(vol, 'descend');
% 
% A = A(sortidx(1:min(length(sortidx),K)),:);
% b = b(sortidx(1:min(length(sortidx),K)));
% 
% mag = sqrt(sum(A.^2,2));
% A = A ./ repmat(mag, 1, size(A,2));
% b = b ./ mag;

% center = (max(P) + min(P))/2;
% A2 = bsxfun(@minus, P, center);
% b2 = -sum(A2.* P,2);
% mag = sqrt(sum(A2.^2,2));
% A2 = A2 ./ repmat(mag, 1, size(A2,2));
% b2 = b2 ./ mag;
% 
% mult = -sign(bsxfun(@minus, A2 * center', b2));
% % A2 = A2 .* repmat(mult, [1,size(A2,2)]);
% b2 = b2 .* mult;
% 
% keyboard
% 
% x = X_raw(1:1000:end,:);
% 
% max(max(bsxfun(@minus, A*x', b), 0), [], 1)
% 
% max(max(bsxfun(@minus, A2*x', b2), 0), [], 1)
% 
% A = [A; A2];
% b = [b; b2];

% [i,j,k] = ndgrid(-6:0.1:1, -6:0.1:1, -6:0.1:1);
% ijk = [i(:), j(:), k(:)];
% D = reshape(max(max(0, bsxfun(@minus, A*ijk', b)), [], 1), size(i));
% D1 = 10*D(:,:,35);
% D2 = bwdist(D1 == 0);

% trisurf(convhulln(X_recon), X_recon(:,1), X_recon(:,2), X_recon(:,3), 'Facecolor','cyan'); axis equal;
% % hold on;
% % trisurf(convhulln(P), P(:,1), P(:,2), P(:,3), 'Facecolor','magenta'); axis equal;
% % hold off;
% drawnow;
