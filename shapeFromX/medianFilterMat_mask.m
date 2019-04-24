function A = medianFilterMat(invalid, half_width)

sz = size(invalid);

width = 2*half_width + 1;
fs = {};
for i = -half_width:half_width
  for j = -half_width:half_width%-half_width:half_width
    if ((i == 0) && (j == 0))
      continue
    end
    f = zeros(width,width);
    f(i+half_width+1, j+half_width+1) = -1;
    f(half_width+1, half_width+1) = 1;

%     f = f / 2;
    
    f = f(find(any(f,2), 1, 'first') : find(any(f,2), 1, 'last'), find(any(f,1), 1, 'first') : find(any(f,1), 1, 'last'));
    fs{end+1} = f;
    
  end
end

do_remove = false(length(fs),1);
for i = 1:length(fs)
  for j = (i+1):length(fs)
    C = conv2(abs(fs{i}), abs(fs{j}), 'valid');
    if any(C>=2)
      do_remove(j) = 1;
    end
  end
end

fs = fs(~do_remove);

A = {};
for i = 1:length(fs)
  A{i} = conv2mat(sz, fs{i});
  f = double(fs{i} ~= 0);
  R = reshape(conv2(double(invalid), f, 'valid'), [], 1);
  keep = R == 0;
  A{i} = A{i}(keep,:);
  
%   [x1, y1] = find(f > 0);
%   [x2, y2] = find(f < 0);
%   [x,y] = bresenham([x1, x2], [y1, y2]);
%   fl = full(sparse(y,x, ones(size(x)), size(f,1), size(f,2)));
%   R = reshape(conv2(double(invalid), fl, 'valid'), [], 1);
%   
%   W = (R==0);
%   
%   A{i} = sparse(1:length(W), 1:length(W), W, length(W), length(W)) * A{i};
%   A{i} = A{i}(any(A{i},2),:);
  
end

A = cat(1,A{:});
