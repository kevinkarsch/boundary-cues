function hinge = getHingeConstraints(contact_pts, horizon, h, w)
    hinge = {};
    
    if( isempty(contact_pts) || isempty(horizon) )
        return;
    end
    
    h_t = horizon(1:2)-horizon(3:4); h_t = h_t./norm(h_t);
    h_n = [h_t(2), -h_t(1)];
    h_d = -dot(h_n, horizon(1:2));
    
    for i=1:size(contact_pts,1)
        for j=i+1:size(contact_pts,1)
            cp_ray = contact_pts(i,:)-contact_pts(j,:);
            cp_ray = cp_ray./norm(cp_ray);
            denom = dot(h_n, cp_ray);
            if(abs(denom)>1e-6)
                ti = -(dot(h_n, contact_pts(i,:)) + h_d) / denom;
                tj = -(dot(h_n, contact_pts(j,:)) + h_d) / denom;
                dist_ij = 10.* (ti/tj - 1);
            else %Points parallel to horizon
                dist_ij = 0;
            end
            hinge{end+1} = []; %#ok<AGROW>
            hinge{end}.mask1 = getHingeMask(contact_pts(i,:), h, w);
            hinge{end}.mask2 = getHingeMask(contact_pts(j,:), h, w);
            hinge{end}.delta = dist_ij;
        end
    end
end

function m = getHingeMask(pt, h, w)
    m = false(h,w);
    m(max(min(round(pt(2))+(-1:1),h),1),max(min(round(pt(1))+(-1:1),w),1)) = true;
end
