function [E, dS_dN, dS_dL] = renderSH_helper(N, L)

try
  
  if nargout >= 3
    [E, dS_dN, dS_dL] = renderSH_helper_mex(N, L);
  elseif nargout >= 2
    [E, dS_dN] = renderSH_helper_mex(N, L);
  else
    [E] = renderSH_helper_mex(N, L);
  end
  
catch ME
  fprintf('ERR in renderSH_helper: %s\n', ME.message);
  
  if nargout >= 3
    [E, dS_dN, dS_dL] = renderSH_helper_mat(N, L);
  elseif nargout >= 2
    [E, dS_dN] = renderSH_helper_mat(N, L);
  else
    [E] = renderSH_helper_mat(N, L);
  end
  
end
