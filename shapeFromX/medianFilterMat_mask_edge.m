function A = medianFilterMat_mask_edge(valid, edges, half_width)
% valid is a binary mask of which pixels are valid
% edges is a binary mask of which pixels should break pairwise connections. Set edges = [] to disable that features.
% half_width is the manhattan distance from which we should take neighbors. (half_width == 2 give all neighbors within a 5x5 patches, for example)

invalid = ~valid;
sz = size(invalid);

width = 2*half_width + 1;
fs = {};
gs = {};
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
    
    [i1,j1] = find(f == 1);
    [i2,j2] = find(f == -1);
    vec = [i2 - i1, j2 - j1];
    vec = vec ./ sqrt(sum(vec.^2,2));
    
    [ii,jj] = ndgrid(1:size(f,1), 1:size(f,2));
    ii = ii(:) - i1;
    jj = jj(:) - j1;
    d = reshape(sqrt(sum(((([ii,jj] * vec') * vec) - [ii,jj]).^2,2)), size(f));
    g = d <= .5;

    gs{end+1} = g;

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
gs = gs(~do_remove);

A = {};
for i = 1:length(fs)
  A{i} = conv2mat(sz, fs{i});
  f = double(fs{i} ~= 0);
  g = double(gs{i});
  keep = reshape(conv2(double(invalid), f, 'valid'), [], 1) == 0;
  if ~isempty(edges)
    keep = keep & (reshape(conv2(double(edges), g, 'valid'), [], 1) == 0);
  end
  
  A{i} = A{i}(keep,:);
  
end

A = cat(1,A{:});
