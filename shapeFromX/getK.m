function [K, dKZ, Zmag] = getK(Z, tap)

getK_filters

if nargin < 2
  tap = GETK_TAP;
end

% Zp = [Z(1,1), Z(1,:), Z(1,end); 
%       Z(:,1), Z,      Z(:,end);
%       Z(end,1), Z(end,:), Z(end,end)];
% 

assert( (tap == 3) || (tap == 4) || (tap == 5) )
if tap == 4
  assert(1==0)
%  Z = conv2(Z, ones(2,2)/4, 'same');
elseif tap == 5
%   [1;4;6;4;1]/(16*sqrt(2))
  Z = conv3(Z, ([1;2;1] * [1,2,1])/12);
%   Z = conv2(Z, ([1;2;1] * [1,2,1])/12, 'same');
end

% Zx = conv2(Z, f1, 'same');
% Zy = conv2(Z, f2, 'same');
% Zyy = conv2(Z, f22, 'same');
% Zxx = conv2(Z, f11, 'same');
% Zxy = conv2(Z, f12, 'same');

% Zx = conv3(Z, f1);
% Zy = conv3(Z, f2);
% Zyy = conv3(Z, f22);
% Zxx = conv3(Z, f11);
% Zxy = conv3(Z, f12);

Zp = pad1(Z);

Zx = conv2(Zp, f1, 'valid');
Zy = conv2(Zp, f2, 'valid');
Zyy = conv2(Zp, f22, 'valid');
Zxx = conv2(Zp, f11, 'valid');
Zxy = conv2(Zp, f12, 'valid');


% 
% f1 = [1, 1;-1,-1]/2;
% f2 = [1,-1; 1,-1]/2;
%     
% f1m  = reshape(f1(end:-1:1), size(f1));
% f2m  = reshape(f2(end:-1:1), size(f2));
% 
% fb = [1,1;1,1]/4;



ZxSq = Zx.^2;
ZySq = Zy.^2;
ZxSq_p = 1 + ZxSq;
ZySq_p = 1 + ZySq;
ZxZy = -2*Zx.*Zy;
ZmagSq = ZxSq_p + ZySq;
Zmag = sqrt(ZmagSq);

denom = max(eps, 2*ZmagSq.*Zmag);
numer = ZxSq_p.*Zyy + ZxZy.*Zxy + ZySq_p.*Zxx;

K = numer ./ denom;

% if nargout >= 3
%   K1 = zeros(size(Z));
%   K2 = zeros(size(Z));
%   Kd = zeros(size(Z));
%   K1(2:end-1,2:end-1) = (ZxSq_p.*Zyy) ./ denom;
%   K2(2:end-1,2:end-1) = (ZySq_p.*Zxx) ./ denom;
%   Kd(2:end-1,2:end-1) = (ZxZy.*Zxy) ./ denom;
% end

if nargout >= 2

  dKZ.denom = denom;
  
  B = 3*(numer./ZmagSq);
  dKZ.f1 = 2*(Zx.*Zyy - Zy.*Zxy) - (Zx.*B);
  dKZ.f2 = 2*(Zy.*Zxx - Zx.*Zxy) - (Zy.*B);
  
  dKZ.f11 = ZySq_p;
  dKZ.f22 = ZxSq_p;
  dKZ.f12 = ZxZy;

end



% function [K, dKZ, K1, K2, Kd] = getK(Z, tap)
% 
% getK_filters
% 
% if nargin < 2
%   tap = GETK_TAP;
% end
% 
% assert( (tap == 3) || (tap == 4) || (tap == 5) )
% if tap == 4
%   Z = conv2(Z, ones(2,2)/4, 'same');
% elseif tap == 5
%   Z = conv2([1,2,1]/4, [1,2,1]'/4, Z, 'same');
% %   Z = conv2(Z, ([1;2;1] * [1,2,1])/16, 'same');
% end
% 
% 
% Zx = conv2(Z, f1, 'same');
% Zy = conv2(Z, f2, 'same');
% Zyy = conv2(Z, f22, 'same');
% Zxx = conv2(Z, f11, 'same');
% Zxy = conv2(Z, f12, 'same');
% 
% ZxSq = Zx.^2;
% ZySq = Zy.^2;
% ZxSq_p = 1 + ZxSq;
% ZySq_p = 1 + ZySq;
% ZxZy = -2*Zx.*Zy;
% ZmagSq = ZxSq_p + ZySq;
% Zmag = sqrt(ZmagSq);
% 
% denom = max(eps, 2*ZmagSq.*Zmag);
% numer = ZxSq_p.*Zyy + ZxZy.*Zxy + ZySq_p.*Zxx;
% 
% K = numer ./ denom;
% 
% if nargout >= 3
%   K1 = zeros(size(Z));
%   K2 = zeros(size(Z));
%   Kd = zeros(size(Z));
%   K1(2:end-1,2:end-1) = (ZxSq_p.*Zyy) ./ denom;
%   K2(2:end-1,2:end-1) = (ZySq_p.*Zxx) ./ denom;
%   Kd(2:end-1,2:end-1) = (ZxZy.*Zxy) ./ denom;
% end
% 
% if nargout >= 2
% 
%   dKZ.denom = denom;
%   
%   B = 3*(numer./ZmagSq);
%   dKZ.f1 = 2*(Zx.*Zyy - Zy.*Zxy) - (Zx.*B);
%   dKZ.f2 = 2*(Zy.*Zxx - Zx.*Zxy) - (Zy.*B);
%   
%   keyboard
%   
%   dKZ.f11 = ZySq_p;
%   dKZ.f22 = ZxSq_p;
%   dKZ.f12 = ZxZy;
% 
% end



% function [K, dKZ, K1, K2, Kd] = getK(Z, tap)
% 
% getK_filters
% 
% if nargin < 2
%   tap = GETK_TAP;
% end
% 
% assert( (tap == 3) || (tap == 4) || (tap == 5) )
% if tap == 4
%   Z = conv2(Z, ones(2,2)/4, 'same');
% elseif tap == 5
%   Z = conv2([1,2,1]/4, [1,2,1]'/4, Z, 'same');
% end
% 
% 
% Zx = conv2(Z, f1, 'valid');
% Zy = conv2(Z, f2, 'valid');
% Zyy = conv2(Z, f22, 'valid');
% Zxx = conv2(Z, f11, 'valid');
% Zxy = conv2(Z, f12, 'valid');
% 
% ZxSq = Zx.^2;
% ZySq = Zy.^2;
% Zxp = (1 + ZxSq);
% Zyp = (1 + ZySq);
% ZxZy = -2*Zx.*Zy;
% ZmagSq = Zxp + ZySq;
% Zmag = sqrt(ZmagSq);
% 
% denom = max(eps, 2*ZmagSq.*Zmag);
% numer = Zxp.*Zyy + ZxZy.*Zxy + Zyp.*Zxx;
% 
% K = zeros(size(Z));
% K(2:end-1, 2:end-1) = numer ./ denom;
% 
% if nargout >= 3
%   K1 = zeros(size(Z));
%   K2 = zeros(size(Z));
%   Kd = zeros(size(Z));
%   K1(2:end-1,2:end-1) = (Zxp.*Zyy) ./ denom;
%   K2(2:end-1,2:end-1) = (Zyp.*Zxx) ./ denom;
%   Kd(2:end-1,2:end-1) = (ZxZy.*Zxy) ./ denom;
% end
% 
% if nargout >= 2
% 
%   dKZ.denom = denom;
%   
%   B = 3*(numer./ZmagSq);
%   dKZ.f1 = 2*(Zx.*Zyy - Zxy.*Zy) - (Zx.*B);
%   dKZ.f2 = 2*(Zxx.*Zy - Zx.*Zxy) - (Zy.*B);
%   
%   dKZ.f11 = Zyp;
%   dKZ.f22 = Zxp;
%   dKZ.f12 = ZxZy;
% 
% end


% function [K, dKZ, K1, K2, Kd] = getK(Z)
% 
% getK_filters
% 
% Zx = conv2(Z, f1, 'valid');
% Zy = conv2(Z, f2, 'valid');
% Zyy = conv2(Z, f22, 'valid');
% Zxy = conv2(Z, f12, 'valid');
% Zxx = conv2(Z, f11, 'valid');
% 
% Zmag = sqrt(1 + Zx.^2 + Zy.^2);
% 
% denom = max(eps, 2*(Zmag.^3));
% numer = (1 + Zx.^2).*Zyy + -2*Zx.*Zy.*Zxy + (1 + Zy.^2).*Zxx;
% 
% K = zeros(size(Z));
% K(2:end-1, 2:end-1) = numer ./ denom;
% 
% if nargout >= 2
% 
%   dKZ.denom = denom;
%   
%   dKZ.f1 = 2*(Zx.*Zyy - Zxy.*Zy) - (Zx.*3.*numer)./(Zmag.^2);
%   dKZ.f2 = 2*(Zxx.*Zy - Zx.*Zxy) - (Zy.*3.*numer)./(Zmag.^2);
%   
%   dKZ.f11 = (1 + Zy.^2);
%   dKZ.f22 = (1 + Zx.^2);
%   dKZ.f12 = (-2*Zx.*Zy);
% 
% end


% function [K, dKZ, K1, K2, Kd] = getK(Z)
% 
% getK_filters
% 
% Zx = conv2(Z, f1, 'valid');
% Zy = conv2(Z, f2, 'valid');
% Zyy = conv2(Z, f22, 'valid');
% Zxy = conv2(Z, f12, 'valid');
% Zxx = conv2(Z, f11, 'valid');
% 
% Zmag = sqrt(1 + (Zx.^2) + (Zy.^2));
% 
% numer = ((1 + (Zx.^2)).*Zyy - 2*(Zx.*Zy).*Zxy + (1 + (Zy.^2)).*Zxx);
% Kv = numer ./ (2.*(Zmag.^3));
% K = zeros(size(Z));
% K(2:end-1, 2:end-1) = Kv;
% 
% 
% if nargout >= 2
%   
%   B = (3/2)*(numer./(Zmag.^5));
%   
%   dKZ.f1 = (Zx.*Zyy - Zxy.*Zy)./(Zmag.^3) - (Zx.*B);
%   
%   dKZ.f2 = (Zxx.*Zy - Zx.*Zxy)./(Zmag.^3) - (Zy.*B);
%   
%   dKZ.f11 = (1 + (Zy.^2))./(2.*(Zmag.^3));
%   dKZ.f22 = (1 + (Zx.^2))./(2.*(Zmag.^3));
%   dKZ.f12 = -((Zx.*Zy)./(Zmag.^3));
% 
% end


