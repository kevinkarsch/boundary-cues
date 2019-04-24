function [results, state, data, params] = myGo(eval_string, labeldata)

%Create valid mask from labeldata
im = im2double(labeldata.im);
[h,w,~] = size(im);
if(islogical(labeldata.extbd))
    mask = labeldata.extbd;
else
    if( ~iscell(labeldata.extbd) ); labeldata.extbd = {labeldata.extbd}; end
    %create initial mask from first contour polygon
    mask = poly2mask(labeldata.extbd{1}(:,1), labeldata.extbd{1}(:,2), h, w);
    for i=2:numel(labeldata.extbd) %subtract all other contours (these are "holes" in the object)
        mask = mask & ~poly2mask(labeldata.extbd{i}(:,1), labeldata.extbd{i}(:,2), h, w);
    end
end

addpath(genpath('./IPPConvsToolbox'));
addpath(genpath('./matlabPyrTools/'));
addpath(genpath('./bmorph/'));

rand('twister',5489);
randn('state',0);

params.EVAL_STRING = eval_string;
CONSTANTS;
eval(params.EVAL_STRING);
myCONSTANTS;
eval(params.EVAL_STRING);

[mask_crop, crop_idx1, crop_idx2] = cropmask(mask);
data.crop_idx1 = crop_idx1; data.crop_idx2 = crop_idx2;
im_crop = im(crop_idx1, crop_idx2, :);
im_crop(repmat(~mask_crop, [1,1,3])) = nan;
data.true.im = im_crop;
data.true.log_im = log(im_crop+0.001);
data.valid = mask_crop;
%Crop labeldata as well
labelfields = fieldnames(labeldata);
crop_offset = [crop_idx2(1)-1, crop_idx1(1)-1];
for i=1:numel(labelfields)
    fldname = labelfields{i};
    if( strcmpi(fldname, 'im') ); continue; end
    if( strcmpi(fldname, 'horizonline') && ~isempty(labeldata.(fldname)))
        labeldata.horizonline = labeldata.horizonline - [crop_offset, crop_offset];
    elseif( ~isempty(labeldata.(fldname)) && ~iscell(labeldata.(fldname)) ) 
        labeldata.(fldname)(:,1:2) = bsxfun(@minus, labeldata.(fldname)(:,1:2), crop_offset);
    else
        for j=1:numel(labeldata.(fldname))
            labeldata.(fldname){j}(:,1:2) = bsxfun(@minus, labeldata.(fldname){j}(:,1:2), crop_offset);
        end
    end
end

%%%

[data.border.mask, data.border.normal] = getBorderNormals(data.valid);

%%%%%%%%%%%%%%%%%% Modifications to handle new constraints - KK 9/13/2012
[h,w] = size(data.valid);

%Remove border normals from mask if exterior boundary is sharp
if(params.USE_SHARP_BDRY)
    [data.border.mask, data.border.normal] = removeSharpBorderNormals(data.border, labeldata.extbd);
end
    
%Append self occlusion info to above mask/normals
if(params.USE_SELF_OCCLUSION)
    [self_occ_mask, self_occ_normals] = getSelfOcclusionNormals(labeldata.intbd, h, w);
    normal_im = zeros([size(data.border.mask),2]);
    normal_im(repmat(data.border.mask,[1,1,2])) = data.border.normal;
    normal_im(repmat(self_occ_mask,[1,1,2])) = self_occ_normals;
    data.border.mask = data.border.mask | self_occ_mask;
    data.border.normal = reshape(normal_im(repmat(data.border.mask,[1,1,2])), [], 2);
else
    self_occ_mask = false(size(data.valid));
end
data.self_occ_mask = self_occ_mask;

%Get fold pixels
data.fold = getFoldPixels(labeldata.convexcnr, labeldata.concavecnr);

%Get hinge distances (contact point constraints)
data.hingeDist = getHingeConstraints(labeldata.contact_pts, labeldata.horizonline, h, w);
%%%%%%%%%%%%%%%%%%%%%%%%%%


%data.Z_median_filter_mat = medianFilterMat_mask_edge(data.valid, self_occ_mask, params.Z_MEDIAN_HALFWIDTH);
data.Z_median_filter_mat = medianFilterMat_mask(~data.valid, params.Z_MEDIAN_HALFWIDTH);
data.A_median_filter_mat = medianFilterMat_mask(~data.valid, params.A_MEDIAN_HALFWIDTH);

data.Z_median_filter_mat_T = data.Z_median_filter_mat';
data.A_median_filter_mat_T = data.A_median_filter_mat';

%%%

load('prior.mat');
data.prior = prior;
for v = params.GLOBAL_VARS
    eval(['global ', v{1}, ';']);
    eval([v{1}, ' = [];']);
end
load('Ldata.mat');
  
if params.DISABLE_LIGHT_PRIOR
    data.L_whiten_params.mean = zeros(1,27);
    data.L_whiten_params.map = eye(27,27);
    data.L_whiten_params.inverse = eye(27,27);
    data.L_gaussian.mu = zeros(1,27);
else  
    data.L_whiten_params = Ldata.whiten.(params.L_WHITEN_PARAMS);
    data.L_gaussian = Ldata.gaussian.(params.L_WHITEN_PARAMS);
end
    
%%%%%%%%%
if params.MAKE_VIDEO
    mkdir video_frames;
    data.video_folder = ['video_frames/', num2str(params.EXPERIMENT), '_', num2str(params.SOLVE_LIGHT), '_', name, '/'];
    mkdir(data.video_folder)
    system(['rm -rf ', data.video_folder, '*']);
    %data.video_log = fopen([data.video_folder, 'log.txt'], 'w');
    global GLOBAL_VIDEO_COUNT; %#ok<*TLEV>
    GLOBAL_VIDEO_COUNT = 1;
end
%%%%%%%%%

start_time = clock;
params.LOSSFUN = 'my_lossfun_sirfs';
[state] = mySolve(data, params);
state.solve_time = etime(clock, start_time);
  
losses = [];  
if params.DO_DISPLAY;
    [~, ~, losses.us.height] = priorZ(state.height, data, params);
    [~, ~, losses.us.albedo] = priorA(log(state.albedo), data, params);
    [~, ~, losses.us.light] = priorL(state.light, data, params);
    dispStruct('losses.us.height', losses.us.height)
    dispStruct('losses.us.albedo', losses.us.albedo)
    dispStruct('losses.us.light', losses.us.light)
    fprintf('Total us loss: %f\n', sum(struct2vector(losses.us)))
end
    
state_init = [];
state_init.light = reshape(data.L_gaussian.mu, [9,3]);
state_init.height = zeros(size(state.height));
[S, ~, N] = renderSH(state_init.height, state_init.light);
A = data.true.log_im - S;       
state_init.shading = exp(S)-0.001;
state_init.albedo = exp(A)-0.001;
state_init.normal = N;
  
results.mult = params.MULT_OPTS;
results.losses = losses;
results.solve_time = state.solve_time;

state_pad = state;
state_pad.height = nan(size(mask));
state_pad.albedo = nan([size(mask),3]);
state_pad.shading = nan([size(mask),3]);

state_pad.height(crop_idx1, crop_idx2) = state.height;
state_pad.albedo(crop_idx1, crop_idx2,:) = state.albedo;
state_pad.shading(crop_idx1, crop_idx2,:) = state.shading;
state_pad.mask = mask;
  
invalid = true(size(mask));
invalid(crop_idx1, crop_idx2) = ~data.valid;
state_pad.height(invalid) = nan;
state_pad.albedo(repmat(invalid,[1,1,3])) = nan;
state_pad.shading(repmat(invalid,[1,1,3])) = nan;
  
state = state_pad;


