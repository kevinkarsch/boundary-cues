function [im, pind] = downsample(im, D, filt, edges)

if nargin < 3
  filt = namedFilter(PYR_FILTER_NAME);
%   filt = namedFilter('binom5')/2;
%   filt = [-0.0761; 0.3536; 0.8593; 0.3536; -0.0761]/16;
%   filt = [1;4;6;4;1]/16;
end

if nargin < 4
  edges = 'zero';
end

pind = zeros(D,2);
for d = 1:D
  
  im_sz = size(im);
  
  pind(d,:) = im_sz;
    
  if d == D
    break;
  end
  
  if (im_sz(2) == 1)
    im = corrDn(im, filt, edges, [2 1], [1 1]);
  elseif (im_sz(1) == 1)
    im = corrDn(im, filt', edges, [1 2], [1 1]);
  else
    im = corrDn(im, filt', edges, [1 2], [1 1]);
    im = corrDn(im, filt, edges, [2 1], [1 1]);
  end

  
%   im = corrDn(im, filt', edges, [1 2], [1 1]);
%   im = corrDn(im, filt, edges, [2 1], [1 1]);
  
end
