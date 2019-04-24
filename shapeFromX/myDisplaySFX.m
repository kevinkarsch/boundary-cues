function v = myDisplaySFX(sfx, outfile)
%     minheight = min(sfx.height(sfx.height(:)>0));
%     maxheight = max(sfx.height(:));
%     visHeight = (sfx.height-minheight)./(maxheight-minheight);
    visHeight = visualizeDEM(sfx.height);
    invalid = repmat(all(visHeight == 1,3), [1,1,3]);
    visHeight(invalid) = 0;
    
%     visNormals = (sfx.normals+1)./2;
    visNormals = visualizeNormals_color(sfx.normals);
    visNormals(invalid) = 0;
    
    h = figure;
    subplot(2,3,1), imshow(sfx.im), title('Image');
    subplot(2,3,2), imshow(visHeight), title('Heightfield');
    subplot(2,3,3), imshow(sfx.albedo), title('Albedo');
    subplot(2,3,4), imshow(sfx.mask), title('Mask');
    subplot(2,3,5), imshow(visNormals), title('Normals');
    subplot(2,3,6), imshow(sfx.shading), title('Shading');
    
    if(nargout>0)
        v = [sfx.im, visHeight, sfx.albedo; repmat(sfx.mask,[1,1,3]), visNormals, sfx.shading];
    end
    if(nargin>1)
        [~,~,ext] = fileparts(outfile);
        print(h, outfile, ['-d' ext(2:end)]);
    end
end
