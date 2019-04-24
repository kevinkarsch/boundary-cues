function [B,N] = getBorderNormals(V, d)

if nargin < 2
  d = 8;
end

B = V & bmorph(~V, [0,1,0;1,1,1;0,1,0]);

[Bi,Bj] = find(B);
[x,y] = meshgrid(-d:d, -d:d);
gaussian = exp(-4/d.^2*(x.^2 + y.^2));

% V_pad = padarray(V, [d,d], false);
V_pad = zeros(size(V)+2*d);
V_pad((d+1):end-d,(d+1):end-d) = V;

N = nan(length(Bi), 2);
for i = 1:length(Bi)
  
  %     patch = V(Bi(i) + [-d:d], Bj(i) + [-d:d])
  patch = V_pad(Bi(i) + [0:2*d], Bj(i) + [0:2*d]);
  patch = patch .* gaussian;
  [pi,pj, v] = find(patch);
  pi = pi - (d+1);
  pj = pj - (d+1);
  N(i,:) = -[mean(pi.*v), mean(pj.*v)];
  
end

N = N ./ repmat(sqrt(sum(N.^2,2)), [1,2]);

