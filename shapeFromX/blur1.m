function Zb = blur1(Z)

Zb = {};
for c = 1:size(Z,3)
  Zb{c} = grid2point(point2grid(Z(:,:,c)));
end
Zb = cat(3, Zb{:});



