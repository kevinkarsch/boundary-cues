function [D, min_idx]=distMat(P1, P2, thresh, interval, conserve_memory)
% squared Euclidian distances between vectors

% P1 = double(P1);
% P2 = double(P2);
% 
% X1=repmat(sum(P1.^2,2),[1 size(P2,1)]);
% X2=repmat(sum(P2.^2,2),[1 size(P1,1)]);
% R=P1*P2';
% D=X1+X2'-2*R;
% return

if nargin >= 3 && ~isempty(thresh)
  do_thresh = 1;
else
  do_thresh = 0;
end

if do_thresh
  D = false(size(P1,1), size(P2,1));
else
  D = zeros(size(P1,1), size(P2,1), class(P1));
end

if nargin < 4
  interval = 500;
end

if nargin < 5
  conserve_memory = 0;
end

persistent X1 X2 R;
if ~conserve_memory
  R=2*P1*P2';
end

last_len = -1;
t_start = cputime;
t_last = t_start;
for i = 1 : interval : size(P1,1)
    rows = i-1 + [1:interval];
    rows = rows(rows <= size(P1,1));

    X1=repmat(sum(P1(rows,:).^2,2),[1 size(P2,1)]);
    if size(X1,1) ~= last_len
        X2=repmat(sum(P2.^2,2),[1 size(X1,1)])';
    end
    last_len = size(X1,1);
    
    if ~conserve_memory
      if do_thresh
        D(rows,:) = (X1+X2-R(rows,:)) <= thresh;
      else
        D(rows,:) = (X1+X2-R(rows,:));
      end
    else
      if do_thresh
        D(rows,:) = (X1+X2-2*P1(rows,:)*P2') <= thresh;
      else
        D(rows,:) = (X1+X2-2*P1(rows,:)*P2');
      end
    end
    
    t = cputime;
    
    if i > 1 && (t - t_last) > 20
      t_last = t;
      fprintf('DISTMAT: %0.2f%% done, %0.1f minutes left\r', i/size(P1,1)*100, (1/60)*(t - t_start)/i*(size(P1,1)-i));
    end
    
end

if nargout >= 2
  D_min = min(D,[],1);
  [i,j] = find(D == repmat(D_min, size(D,1), 1));
  [uu, ui] = unique(j, 'first');
  min_idx = i(ui);
end


end


