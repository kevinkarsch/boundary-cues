function [S, dS_Z, N, dN_Z] = renderSH(Z, L)

assert(size(L,1) == 9)
assert((size(L,2) == 1) || (size(L,2) == 3))

[N, dN_Z] = getNormals_conv(Z);
N_vec = reshape(N, [], 3);

S = {};
dS_Z = {};
for c = 1:size(L,2)
  [s, ds_n] = renderSH_helper(N_vec, L(:,c));
  S{c} = reshape(s, size(Z));
  dS_N = reshape(ds_n, [size(Z), 3]);
  
  dS_Z{c}.F1 = dN_Z.F1_1 .* dS_N(:,:,1) + dN_Z.F1_2 .* dS_N(:,:,2) + dN_Z.F1_3 .* dS_N(:,:,3);
  dS_Z{c}.F2 = dN_Z.F2_1 .* dS_N(:,:,1) + dN_Z.F2_2 .* dS_N(:,:,2) + dN_Z.F2_3 .* dS_N(:,:,3);
end
S = cat(3, S{:});



% function [S, dS_Z, dS_L] = renderSH(Z, L, grid_neighbors)
% 
% if nargin < 3
%   grid_neighbors = imageNeighbors(size(Z)+1);
% end
% 
% RETURN_GRADIENT = (nargout >= 2); 
% 
% Z_grid = point2grid(Z);
% 
% if RETURN_GRADIENT
%   
%   [N1, N2, dN1_Z, dN2_Z] = getNormals_safe(Z_grid, grid_neighbors);
%   
%   [S1, dS1_L, dS1_N] = renderSH_helper(N1, L);
%   [S2, dS2_L, dS2_N] = renderSH_helper(N2, L);
% 
%   S1_valid = S1 > 0;
%   S2_valid = S2 > 0;
%     
%   dS_Z.Z00 = 0.5 * (S1_valid .* sum([dN1_Z.X00, dN1_Z.Y00, dN1_Z.Z00] .* dS1_N,2));
%   dS_Z.Z01 = 0.5 * ((S1_valid .* sum([dN1_Z.X01, dN1_Z.Y01, dN1_Z.Z01] .* dS1_N,2)) + (S2_valid .* sum([dN2_Z.X01, dN2_Z.Y01, dN2_Z.Z01] .* dS2_N,2)));
%   dS_Z.Z10 = 0.5 * ((S1_valid .* sum([dN1_Z.X10, dN1_Z.Y10, dN1_Z.Z10] .* dS1_N,2)) + (S2_valid .* sum([dN2_Z.X10, dN2_Z.Y10, dN2_Z.Z10] .* dS2_N,2)));
%   dS_Z.Z11 = 0.5 * (S2_valid .* sum([dN2_Z.X11, dN2_Z.Y11, dN2_Z.Z11] .* dS2_N,2));
% 
%   dS_L = (dS1_L + dS2_L)/2;
% 
% else
%   
%   [N1, N2] = getNormals_safe(Z_grid, grid_neighbors);
%   
%   [S1] = renderSH_helper(N1, L);
%   [S2] = renderSH_helper(N2, L);
% 
% end
% 
% S = reshape(0.5*(S1 + S2), size(Z));
% 
% 
% 
