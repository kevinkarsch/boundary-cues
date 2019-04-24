function [E, dS_dN, dS_dL] = renderSH_helper(N, L, NO_SHADOW)

RETURN_DE = nargout >= 2;
RETURN_DL = nargout >= 3;

ORDER = 2;

if ORDER == 1  
  L(5:end) = 0;
end


c1 = 0.429043;
c2 = 0.511664;
c3 = 0.743125;
c4 = 0.886227;
c5 = 0.247708;

M = [ c1 * L(9),  c1 * L(5),  c1 * L(8),  c2 * L(4);
      c1 * L(5), -c1 * L(9),  c1 * L(6),  c2 * L(2);
      c1 * L(8),  c1 * L(6),  c3 * L(7),  c2 * L(3);
      c2 * L(4),  c2 * L(2),  c2 * L(3),  c4 * L(1) - c5 * L(7)];

assert(size(N,2) == 3)

I = ones(size(N,1), 1);
Na = [N, I];
NM = Na * M;
E = sum(NM .* Na,2);

% if (nargin < 3) || NO_SHADOW
%   E = max(0, E);
% end
% 
% iszero = E == 0;

if RETURN_DL
  dS_dL = [c4 * I, (2*c2)*N(:,2), (2*c2)*N(:,3), (2*c2)*N(:,1), (2*c1)*(N(:,1).*N(:,2)), (2*c1)*(N(:,2).*N(:,3)), c3*(N(:,3).^2) - c5, (2*c1)*(N(:,1).*N(:,3)), c1*(N(:,1).^2 - N(:,2).^2)];
%   dS_dL(iszero,:) = 0;
end

if RETURN_DE
  dS_dN = 2*NM(:,1:3);
%   dS_dN(iszero,:) = 0;
end




% function [E, dS_dL, dS_dN] = renderSH_helper(N, L)
% 
% RETURN_DL = nargout >= 2;
% RETURN_DE = nargout >= 3;
% 
% c1 = 0.429043;
% c2 = 0.511664;
% c3 = 0.743125;
% c4 = 0.886227;
% c5 = 0.247708;
% 
% M = [ c1 * L(9),  c1 * L(5),  c1 * L(8),  c2 * L(4);
%       c1 * L(5), -c1 * L(9),  c1 * L(6),  c2 * L(2);
%       c1 * L(8),  c1 * L(6),  c3 * L(7),  c2 * L(3);
%       c2 * L(4),  c2 * L(2),  c2 * L(3),  c4 * L(1) - c5 * L(7)];
% 
% 
% assert(size(N,2) == 3)
% 
% I = ones(size(N,1), 1);
% Na = [N, I];
% NM = Na * M;
% E = max(0, sum(NM .* Na,2));
% 
% % Mrep = repmat(M(end,:), size(N,1),1);
% % NM = N * M(1:3,:) + Mrep;
% % E = max(0, sum(NM(:,1:3) .* N,2) + NM(:,4));
% 
% 
% % E = sum(NM .* Na,2);
% iszero = E == 0;
% 
% if RETURN_DL
%   dS_dL = [c4 * I, (2*c2)*N(:,2), (2*c2)*N(:,3), (2*c2)*N(:,1), (2*c1)*(N(:,1).*N(:,2)), (2*c1)*(N(:,2).*N(:,3)), c3*(N(:,3).^2) - c5, (2*c1)*(N(:,1).*N(:,3)), c1*(N(:,1).^2 - N(:,2).^2)];
%   %     E = max(0, dS_dL * L);
%   dS_dL(iszero,:) = 0;
% end
% 
% if RETURN_DE
%   dS_dN = 2*NM(:,1:3);
% %   dS_dN = N * M(:,1:3);
%   dS_dN(iszero,:) = 0;
% end
% 
% 
% % NORMALIZE = false;
% %
% % if NORMALIZE
% %   % Divide by the max
% %
% %   E_unsigned = sum((N * M) .* N,2);
% %   E = max(0, E_unsigned);
% %   nonzero = E > 0;
% %
% %   Z_idx = find(E == max(E), 1, 'first');
% %   Z = max(eps, E(Z_idx));
% %
% %   %   Z2 = (4.*(c2.*L(3)+c1.*(L(6).^2+L(8).^2).^(1/2)).^2+(2.*c2.*(L(2).^2+L(4).^2).^(1/2)+c1.*(L(5).^2+L(9).^2).^(1/2)).^2).^(1/2)+(1/2).*(2.*c4.*L(1)+(c3+(-2).*c5).*L(7)+abs(c3.*L(7)));
% %   %   [Z2/Z]
% %
% %   if RETURN_DL
% %     dL = [c4 * (N(:,4).^2), (2*c2)*(N(:,2).*N(:,4)), (2*c2)*(N(:,3).*N(:,4)), (2*c2)*(N(:,1).*N(:,4)), (2*c1)*(N(:,1).*N(:,2)), (2*c1)*(N(:,2).*N(:,3)), c3*(N(:,3).^2) - c5*(N(:,4).^2), (2*c1)*(N(:,1).*N(:,3)), c1.*(N(:,1).^2+(-1).*N(:,2).^2)];
% %     dL = dL/Z - E * (dL(Z_idx,:)/(Z^2));
% %     dL = repmat(nonzero, 1, size(dL,2)) .* dL;
% %   end
% %
% %   if RETURN_DE
% %     dE = (2/Z)*(N * M(:, 1:3));
% %     dE = repmat(nonzero, 1, size(dE,2)) .* dE;
% %   end
% %
% %   E = E/Z;
% %
