function [mask,normals] = getSelfOcclusionNormals(occ, h, w)

%Rasterize occluding contours to mask, calculate normals
mask = false(h,w);
nx = zeros(h,w);
ny = zeros(h,w);
for i=1:numel(occ)
    for j=1:size(occ{i},1)-1
        n_ij = [occ{i}(j,2)-occ{i}(j+1,2), occ{i}(j,1)-occ{i}(j+1,1)];
        numPts = 2*norm(n_ij);
        x = linspace(occ{i}(j,1), occ{i}(j+1,1), numPts);
        y = linspace(occ{i}(j,2), occ{i}(j+1,2), numPts);
        occ_idx = sub2ind([h,w],round(y),round(x));
        mask(occ_idx) = true;
        n_ij = n_ij./norm(n_ij);
        nx(occ_idx) = n_ij(2);
        ny(occ_idx) = -n_ij(1);
    end
end
normals = [nx(mask), ny(mask)];
