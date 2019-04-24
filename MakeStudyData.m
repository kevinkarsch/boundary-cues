function MakeStudyData( labelmat )
    
    fprintf('\n\n PROCESSING: %s \n\n', labelmat);

    load(labelmat);
    [path, matname, ~] = fileparts(labelmat);
    
%if(exist(fullfile(path,[matname '-result.mat']),'file')); load(fullfile(path,[matname '-result.mat'])); end
    
    cd shapeFromX/;
    
    %%%%%%%%%%%%%%% SFC all cues %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfc_all = 'Shape from contour + all cues';
    fprintf('Solving %s\n', heightname.sfc_all); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 1;', ...
                    'params.USE_SHARP_BDRY = 1;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 1.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };']; %Not using contact points
    [sfx info] = myShapeFromX(param_string, labeldata, false);
    height.sfc_all = sfx.height;
    
    %%%%%%%%%%%%%%% Baron+Malik SFS %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfs = 'Shape from shading (Baron+Malik)';
    fprintf('Solving %s\n', heightname.sfs); 
    param_string = ['params.SOLVE_ALBEDO = 1;', ...
                    'params.SOLVE_LIGHT = 1;', ...
                    'params.USE_SELF_OCCLUSION = 0;', ...
                    'params.USE_SHARP_BDRY = 0;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 0.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfs = sfx.height;
    
    %%%%%%%%%%%%%%% SFS all cues %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfs_all = 'Shape from shading + all cues';
    fprintf('Solving %s\n', heightname.sfs_all); 
    param_string = ['params.SOLVE_ALBEDO = 1;', ...
                    'params.SOLVE_LIGHT = 1;', ...
                    'params.USE_SELF_OCCLUSION = 1;', ...
                    'params.USE_SHARP_BDRY = 1;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 1.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };']; %Not using contact points
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfs_all = sfx.height;
    
    %%%%%%%%%%%%%%% SFC %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfc = 'Shape from contour';
    fprintf('Solving %s\n', heightname.sfc); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 0;', ...
                    'params.USE_SHARP_BDRY = 0;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 0.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfc = sfx.height;
    
    %%%%%%%%%%%%%%% SFC + smooth/sharp %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfc_ss = 'Shape from contour + sharp bdry';
    fprintf('Solving %s\n', heightname.sfc_ss); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 0;', ...
                    'params.USE_SHARP_BDRY = 1;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 0.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfc_ss = sfx.height;
    
    %%%%%%%%%%%%%%% SFC + int bdry %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfc_selfocc = 'Shape from contour + self occlusions';
    fprintf('Solving %s\n', heightname.sfc_selfocc); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 1;', ...
                    'params.USE_SHARP_BDRY = 0;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 0.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfc_selfocc = sfx.height;
    
    %%%%%%%%%%%%%%% SFC + self occlusion + smooth/sharp %%%%%%%%%%%%%%%%%
    heightname.sfc_selfocc_ss = 'Shape from contour + self occlusions + sharp bdry';
    fprintf('Solving %s\n', heightname.sfc_selfocc_ss); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 1;', ...
                    'params.USE_SHARP_BDRY = 1;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 0.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfc_selfocc_ss = sfx.height;
    
    %%%%%%%%%%%%%%% SFC + folds %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfc_folds = 'Shape from contour + folds';
    fprintf('Solving %s\n', heightname.sfc_folds); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 0;', ...
                    'params.USE_SHARP_BDRY = 1;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 1.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 0.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfc_folds = sfx.height;
    
    %%%%%%%%%%%%%%% SFC + contact pts %%%%%%%%%%%%%%%%%%%%%%
    heightname.sfc_contactpts = 'Shape from contour + contact pts';
    fprintf('Solving %s\n', heightname.sfc_contactpts); 
    param_string = ['params.SOLVE_ALBEDO = 0;', ...
                    'params.SOLVE_LIGHT = 0;', ...
                    'params.USE_SELF_OCCLUSION = 0;', ...
                    'params.USE_SHARP_BDRY = 1;', ...
                    'params.MULT_OPTS.saifs.height.fold = { 0.0 };', ...
                    'params.MULT_OPTS.saifs.height.hinge = { 1.0 };'];
    sfx = myShapeFromX(param_string, labeldata, false);
    height.sfc_contactpts = sfx.height;
        
    %Save visualization
    ha = tight_subplot(3,3,[.03 .03],[.01 .1],[.01 .01]);
    axes(ha(1)); imshow(visualizeDEM(height.sfc)); title('sfc');
    axes(ha(2)); imshow(visualizeDEM(height.sfc_selfocc)); title('sfc selfocc');
    axes(ha(3)); imshow(visualizeDEM(height.sfc_ss)); title('sfc ss');
    axes(ha(4)); imshow(visualizeDEM(height.sfc_selfocc_ss)); title('sfc selfocc ss');
    axes(ha(5)); imshow(visualizeDEM(height.sfc_folds)); title('sfc folds');
    axes(ha(6)); imshow(visualizeDEM(height.sfc_contactpts)); title('sfc contactpts');
    axes(ha(7)); imshow(visualizeDEM(height.sfc_all)); title('sfc all');
    axes(ha(8)); imshow(visualizeDEM(height.sfs)); title('sfs');
    axes(ha(9)); imshow(visualizeDEM(height.sfs_all)); title('sfs all');
    
    cd ..;
    
    print('-dpdf', fullfile(path,[matname '-result.pdf']));
    close;
    %%%%%%%%%%%%%%% Create/save result struct %%%%%%%%%%%%%%%%%
    im = labeldata.im;
    mask = sfx.mask;
    annotations.self_occ = false(size(mask));
    annotations.self_occ(info.data.crop_idx1,info.data.crop_idx2) = info.data.self_occ_mask;
    
    acont_pts = false(size(info.data.valid));
    for i=1:numel(info.data.hingeDist)
        acont_pts = acont_pts | info.data.hingeDist{i}.mask1;
        acont_pts = acont_pts | info.data.hingeDist{i}.mask2;
    end
    annotations.contact_pts = false(size(mask));
    annotations.contact_pts(info.data.crop_idx1,info.data.crop_idx2) = acont_pts;
        
    afolds = false(size(info.data.valid));
    for i=1:numel(info.data.fold)
        p = info.data.fold{i}.position;
        afolds(sub2ind(size(afolds), round(p(:,1)), round(p(:,2)))) = true;
    end
    annotations.folds = false(size(mask));
    annotations.folds(info.data.crop_idx1,info.data.crop_idx2) = afolds;
    
    %save(fullfile(path,[matname '-result.mat']), '-append', 'im', 'mask', 'height', 'heightname', 'labeldata', 'annotations'); 
    save(fullfile(path,[matname '-result.mat']), 'im', 'mask', 'height', 'heightname', 'labeldata', 'annotations'); 
end

