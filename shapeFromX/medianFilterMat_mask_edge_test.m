

clear all;


[i,j] = ndgrid(-10:10, -10:10);
valid = (i.^2 + j.^2) <= 80;

edges = (j == 2) | (i == -1);

half_width = 2;


M = medianFilterMat_mask_edge(valid, edges, half_width);

[idx,junk] = find(M');
ii = idx(1:2:end);
jj = idx(2:2:end);

[i1, j1] = ind2sub(size(valid), ii);
[i2, j2] = ind2sub(size(valid), jj);

imagesc(double(valid) + double(valid .* edges)); colormap('gray');
hold on;
colors = 'bgrcmy';
                        
for i = 1:length(i1)
  line([j1(i),j2(i)] + randn(1,2)/20, [i1(i),i2(i)] + randn(1,2)/20, 'color', colors(mod(i-1, length(colors))+1))
end
axis image;
