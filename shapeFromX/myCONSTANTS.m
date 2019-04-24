%Shape from contour always enabled
params.SOLVE_SHAPE = 1;

params.MULT_OPTS.saifs.height.contour =                     { 0.7071 }; %1
params.MULT_OPTS.saifs.height.slant =                       { 1.0905 }; %1
params.MULT_OPTS.saifs.height.smooth =                      { 3.668 }; %1
params.DISABLE_LIGHT_PRIOR = 1;

params.USE_SELF_OCCLUSION = 1;
params.USE_SHARP_BDRY = 1;
params.MULT_OPTS.saifs.height.fold =                        { 1.0 };
params.MULT_OPTS.saifs.height.hinge =                       { 1.0 };

%Shape from shading
if(params.SOLVE_ALBEDO)
    params.DISABLE_LIGHT_PRIOR = 0;
    params.ALBEDO_MODEL = 'ours';
    params.MULT_OPTS.saifs.albedo.entropy_sigma =           { 0.0625 };
    params.MULT_OPTS.saifs.albedo.entropy =                 { 1.1892 };
    params.MULT_OPTS.saifs.albedo.hist =                    { 9.5137 };
    params.MULT_OPTS.saifs.albedo.smooth =                  { 6.1688 };
    params.ALBEDO_SMOOTH_MODEL = 0;
else
    params.MULT_OPTS.saifs.albedo.entropy =                 { 0 };
    params.MULT_OPTS.saifs.albedo.hist =                    { 0 };
    params.MULT_OPTS.saifs.albedo.smooth =                  { 0 };
end

%Solve for illumination and shape
if(params.SOLVE_LIGHT)
    params.DISABLE_LIGHT_PRIOR = 0;
    params.FORCE_MONOCHROME_LIGHT = 0;
    params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 11.3137 };
    params.MULT_OPTS.saifs.light.lab_color.GSM =            { 0 };
    params.MULT_OPTS.saifs.light.lab_color.MOG =            { 0 }; %{ 12.3377 };
    params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 2.181 };
    params.MULT_OPTS.saifs.light.natural_color.GSM =        { 0 };
    params.MULT_OPTS.saifs.light.natural_color.MOG =        { 0 };
    params.WHITEN_LIGHT = 1;
    params.L_WHITEN_PARAMS = 'lab_color'; %'natural_color' %'natural_gray', 'lab_gray', 'lab_color'
else
    params.MULT_OPTS.saifs.light.lab_color.gaussian =       { 0 };
    params.MULT_OPTS.saifs.light.lab_color.GSM =            { 0 };
    params.MULT_OPTS.saifs.light.lab_color.MOG =            { 0 };
    params.MULT_OPTS.saifs.light.natural_color.gaussian =   { 0 };
    params.MULT_OPTS.saifs.light.natural_color.GSM =        { 0 };
    params.MULT_OPTS.saifs.light.natural_color.MOG =        { 0 };
end

params.MINFUNC_DISPLAY = ''; %'OFF', 'FINAL', or leave blank for per-iteration details
