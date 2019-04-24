function yuv = rgb2yuv(rgb)
%     Wr = 0.299;
%     Wb = 0.114;
%     Wg = 0.587;
%     Umax = 0.436;
%     Vmax = 0.615;
% 
%     rgb = double(rgb);
%     y = Wr * rgb(:,:,1) + Wg * rgb(:,:,2) + Wb * rgb(:,:,3);
%     u = (Umax/(1 - Wb)) * (rgb(:,:,3) - y);
%     v = (Vmax/(1 - Wr)) * (rgb(:,:,1) - y);
%     
%     yuv = cat(3, y, u, v);
    
    A = [ 0.299,    0.587,    0.114;
         -0.14713, -0.28886,  0.436;
          0.615,   -0.51499, -0.10001];
    yuv = reshape(reshape(rgb, [], 3) * A', size(rgb));
    
end