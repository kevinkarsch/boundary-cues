function [new_mask, new_normals] = removeSharpBorderNormals(border, extbd)
    [h,w] = size(border.mask);
    smooth_mask = false(h,w);
    for i=1:numel(extbd)
        for j=1:size(extbd{i},1)-1
            %Only remove boundary if both boundary vertices are labeled 'smooth'
            if( extbd{i}(j,3)~=1 || extbd{i}(j+1,3)~=1 )
                continue;
            end
            %Remove smooth boundary normals
            numPts = 2*norm( extbd{i}(j,:) - extbd{i}(j+1,:) );
            x = linspace(extbd{i}(j,1), extbd{i}(j+1,1), numPts);
            y = linspace(extbd{i}(j,2), extbd{i}(j+1,2), numPts);
            smooth_idx = sub2ind([h,w],clip(round(y),1,h),clip(round(x),1,w));       
            for k=-1:1;
                blk3x3_smoothidx = bsxfun(@plus, smooth_idx, k*h + (-1:1)');
                blk3x3_smoothidx(blk3x3_smoothidx<1) = [];
                blk3x3_smoothidx(blk3x3_smoothidx>w*h) = [];
                smooth_mask( blk3x3_smoothidx(:) ) = true;
            end
        end
    end
    new_mask = border.mask & ~smooth_mask;
    normals = zeros([h,w,2]);
    normals(repmat(border.mask,[1,1,2])) = border.normal;
    new_normals = reshape(normals(repmat(new_mask,[1,1,2])),[],2);
end

function y = clip(x,lb,ub)
    y=x;
    y(x<lb) = lb;
    y(x>ub) = ub;
end