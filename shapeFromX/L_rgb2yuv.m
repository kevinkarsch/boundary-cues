function Ls_yuv = L_rgb2yuv(Ls)

Y = mean(Ls,2);
U = Ls(:,3,:) - Ls(:,2,:);
V = Ls(:,1,:) - Ls(:,2,:);

Ls_yuv = [Y, U, V];
