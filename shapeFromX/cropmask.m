function [cm, crop_idx1, crop_idx2] = cropmask( mask, img )
%CROPMASK Returns the cropped version of the mask (cm), and an axis-aligned
%bounding box of the crop [xmin xmax ymin ymax] (aabb)
    binmask = imdilate(mask(:,:,1)>=1, strel('disk',1));
    [X,Y] = meshgrid(1:size(binmask,2), 1:size(binmask,1));
    X = X.*binmask;
    Y = Y.*binmask;
    aabb = [min(X(X(:)>0)), max(X(:)), min(Y(Y(:)>0)), max(Y(:))];
    crop_idx1 = aabb(3):aabb(4);
    crop_idx2 = aabb(1):aabb(2);
    cm = mask(crop_idx1, crop_idx2, :);
    if(nargin>1)
        cm = img(crop_idx1, crop_idx2, :); %If input is an image, return the cropped version
    end
end
