function fold = getFoldPixels(fold_convex, fold_concave)
    fold = {};
    for i=1:numel(fold_convex)
        fold = [fold; getFoldInfo(fold_convex{i}, true)]; %#ok<*AGROW>
    end
    for i=1:numel(fold_concave)
        fold = [fold; getFoldInfo(fold_concave{i}, false)];
    end
end

function info = getFoldInfo(fold_i, convex)
    info = cell(size(fold_i,1)-1, 1);
    for j=1:size(fold_i,1)-1
        e_p = [fold_i(j,1)-fold_i(j+1,1), fold_i(j,2)-fold_i(j+1,2)]; %Vector parallel to fold
        if(e_p(1) < 0)
            e_p = -e_p;
        elseif(e_p(1) == 0)
            e_p(2) = abs(e_p(2));
        end
        numPts = 2*norm(e_p);
        x = linspace(fold_i(j,1), fold_i(j+1,1), numPts);
        y = linspace(fold_i(j,2), fold_i(j+1,2), numPts);
        t = e_p ./ sqrt(sum(e_p.*e_p));
        
        info{j}.position = [y' x'];
        info{j}.tangent = repmat(fliplr(t), [numel(x), 1]);
        info{j}.normal = [-info{j}.tangent(:,2), info{j}.tangent(:,1)]; % 90deg rotation of the tangent
        if(~convex)
            info{j}.normal = -info{j}.normal;
        end
    end
end
