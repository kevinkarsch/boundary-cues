function [sfx, info] = myShapeIllumReflFromShading(labeldata, display)
    param_string = ['params.SOLVE_ALBEDO = 1;', ...
                    'params.SOLVE_LIGHT = 1;'];
    if(~exist('display','var'))
        display = false;
    end
    [sfx, info] = myShapeFromX(param_string, labeldata, display);
end
