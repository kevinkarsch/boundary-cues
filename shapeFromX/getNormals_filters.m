
% f1 = [0, 0, 0;
%       0, 1, 1;
%       0,-1,-1]/4;
% 
% f2 = [0, 0, 0;
%       0, 1,-1;
%       0, 1,-1]/4;


f1 = [1, 2, 1;
      0, 0, 0;
     -1, -2, -1]/8;

f2 = [1, 0, -1;
      2, 0, -2;
      1, 0, -1]/8;

    
% f1 = [0, 1, 0;
%       0, 0, 0;
%       0,-1, 0]/2;
% 
% f2 = [0, 0, 0;
%       1, 0, -1;
%       0, 0, 0]/2;


f1m  = reshape(f1(end:-1:1), [3,3]);
f2m  = reshape(f2(end:-1:1), [3,3]);
