function [mask,normals] = getSelfOcclusionNormals(occ, h, w)

%Rasterize occluding contours to mask, calculate normals
mask = false(h,w);
nx = zeros(h,w);
ny = zeros(h,w);
for i=1:numel(occ)
    for j=1:numel(occ(i).x)-1
        n_ij = [occ(i).y(j)-occ(i).y(j+1), occ(i).x(j)-occ(i).x(j+1)];
        numPts = 2.*norm(n_ij);
        x = linspace(occ(i).x(j), occ(i).x(j+1), numPts);
        y = linspace(occ(i).y(j), occ(i).y(j+1), numPts);
        occ_idx = sub2ind([h,w],round(y),round(x));
        mask(occ_idx) = true;
        n_ij = n_ij./norm(n_ij);
        nx(occ_idx) = n_ij(2);
        ny(occ_idx) = -n_ij(1);
    end
end
normals = [nx(mask), ny(mask)];

% %Use Barron's code to compute normals at occluding contours
% B = mask;
% [Bi,Bj] = find(B);
% d=8;
% [x,y] = meshgrid(-d:d, -d:d);
% gaussian = exp(-4/d.^2*(x.^2 + y.^2));
% 
% % V_pad = padarray(V, [d,d], false);
% V_pad = zeros(size(B)+2*d);
% V_pad((d+1):end-d,(d+1):end-d) = B;
% 
% N = nan(length(Bi), 2);
% for i = 1:length(Bi)
%   
%   %     patch = V(Bi(i) + [-d:d], Bj(i) + [-d:d])
%   patch = V_pad(Bi(i) + [0:2*d], Bj(i) + [0:2*d]);
%   patch = patch .* gaussian;
%   [pi,pj, v] = find(patch);
%   pi = pi - (d+1);
%   pj = pj - (d+1);
%   N(i,:) = -[mean(pi.*v), mean(pj.*v)];
%   
% end
% N = N ./ repmat(sqrt(sum(N.^2,2)), [1,2]);
% 
% %Make sure normals are oriented properly
% TODO