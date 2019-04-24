function Ls = L_yuv2rgb(Ls_yuv)

G = Ls_yuv(:,1,:) - (Ls_yuv(:,2,:) + Ls_yuv(:,3,:))/3;
R = Ls_yuv(:,3,:) + G;
B = Ls_yuv(:,2,:) + G;

Ls = [R,G,B];
