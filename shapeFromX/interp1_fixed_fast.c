#include "math.h"
#include "mex.h"

void mexFunction(int nargout, mxArray *varargout[], int nargin, const mxArray *varargin[])
{
  
  mxArray *X_mat, *LUT_F_mat;
  double *X, *LUT_F;
  double low_val, bin_width, v, v1, v2;
  int n, LUT_n, i, below_bin, above_bin;
  double *F, *dF;
  double above_F, above_dF, below_F, below_dF;
  bool return_dF;
  
  return_dF = nargout >= 2;
  
  X_mat = varargin[0];
  X = mxGetPr(X_mat);
  
  low_val = *mxGetPr(varargin[1]);
  bin_width = *mxGetPr(varargin[2]);
  
  LUT_F_mat = varargin[3];
  LUT_F = mxGetPr(LUT_F_mat);
  
  n = mxGetDimensions(X_mat)[0];
  LUT_n = mxGetDimensions(LUT_F_mat)[0];
  
  above_F = LUT_F[LUT_n-1];
  above_dF = (LUT_F[LUT_n-1] - LUT_F[LUT_n-2])/bin_width;
  
  below_F = LUT_F[0];
  below_dF = (LUT_F[1] - LUT_F[0])/bin_width;

  varargout[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  F = mxGetPr(varargout[0]);

  if(return_dF){
    varargout[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    dF = mxGetPr(varargout[1]);
  }

  for( i = 0; i < n; i = i + 1){

    v = (X[i]-low_val) / bin_width;
    
    if (v >= (LUT_n-1)){ /* Above */
      
      v = v - (LUT_n - 1);
      F[i] = LUT_F[LUT_n-1] + bin_width*v*above_dF;
      
      if(return_dF){
        dF[i] = above_dF;
      }
      
    }
    else if (v <= 0){ /* Below */
      
      F[i] = LUT_F[0] + bin_width*v*below_dF;
      
      if(return_dF){
        dF[i] = below_dF;
      }
      
    }
    else{ /* in between */
      
      below_bin = floor(v);
      above_bin = below_bin+1;

      v1 = v - below_bin;
      v2 = 1.0 - v1;
      
      F[i] = LUT_F[below_bin] * v2 + LUT_F[above_bin] * v1;
      
      if(return_dF){
        dF[i] = (LUT_F[above_bin] - LUT_F[below_bin])/(bin_width);
      }
      
    } 
  }  
}
