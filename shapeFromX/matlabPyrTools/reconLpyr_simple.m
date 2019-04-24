function res = reconLpyr(pyr, filt2, edges)

res = pyr{end};
for d = (length(pyr)-1):-1:1
  sz_up = size(pyr{d+1});
  sz_down = size(pyr{d});
  
%   res = upConv(res, filt2, edges, [2 1], [1 1], [sz_down(1), sz_up(2)]);
%   res = upConv(res, filt2', edges, [1 2], [1 1], sz_down);
  
  if (sz_up(1) == 1)
    res = upConv(res, filt2', edges, [1 2], [1 1], [sz_up(1), sz_down(2)]);
  elseif (sz_up(2) == 1)
    res = upConv(res, filt2, edges, [2 1], [1 1], [sz_down(1), sz_up(2)]);
  else
    res = upConv(res, filt2, edges, [2 1], [1 1], [sz_down(1), sz_up(2)]);
    res = upConv(res, filt2', edges, [1 2], [1 1], sz_down);
  end
  
  res = res + pyr{d};  
end

% res_sz = ind(1,:);
% 
% if any(levs > 1)
% 
%   int_sz = [ind(1,1), ind(2,2)];
%   
%   nres = reconLpyr( pyr(prod(res_sz)+1:size(pyr,1)), ...
%       ind(2:size(ind,1),:), levs-1, filt2, edges);
%   
%   if (res_sz(1) == 1)
%     res = upConv(nres, filt2', edges, [1 2], [1 1], res_sz);
%   elseif (res_sz(2) == 1)
%     res = upConv(nres, filt2, edges, [2 1], [1 1], res_sz);
%   else
%     hi = upConv(nres, filt2, edges, [2 1], [1 1], int_sz);
%     res = upConv(hi, filt2', edges, [1 2], [1 1], res_sz);
%   end
% 
% else
%   
%   res = zeros(res_sz);
% 
% end
% 
% if any(levs == 1)
%   res = res + pyrBand(pyr,ind,1);
% end
