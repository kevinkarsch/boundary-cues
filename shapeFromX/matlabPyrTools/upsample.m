function im = upsample(im, pind, filt, edges)

if nargin < 3
  filt = namedFilter(PYR_FILTER_NAME)*2;
%   filt = namedFilter('binom5');
%   filt = [1;4;6;4;1]/8;
end

if nargin < 4
  edges = 'zero';
end

assert(all(size(im) == pind(end,:)))

for d = (size(pind,1)-1):-1:1
  
  sz_down = size(im);  
  sz_up = pind(d,:);

  
  if (sz_down(1) == 1)
    im = upConv(im, filt', edges, [1 2], [1 1], [sz_down(1), sz_up(2)]);
  elseif (sz_down(2) == 1)
    im = upConv(im, filt, edges, [2 1], [1 1], [sz_up(1), sz_down(2)]);
  else
    im = upConv(im, filt, edges, [2 1], [1 1], [sz_up(1), sz_down(2)]);
    im = upConv(im, filt', edges, [1 2], [1 1], sz_up);
  end

  
%   im = upConv(im, filt, edges, [2 1], [1 1], [sz_up(1), sz_down(2)]);
%   im = upConv(im, filt', edges, [1 2], [1 1], sz_up);
    
end
