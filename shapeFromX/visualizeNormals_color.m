function out = visualizeNormals(N)

N = N(:,:,[3,1,2]);
N(:,:,2:3) = N(:,:,2:3) / 1.25;
V = max(0, min(1, yuv2rgb_simple(N)));
out = V;
