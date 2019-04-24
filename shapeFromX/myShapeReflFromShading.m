function [sfx, info] = myShapeReflFromShading(labeldata, display)
    param_string = ['params.SOLVE_ALBEDO = 1;', ...
                    'params.SOLVE_LIGHT = 0;'];
    if(~exist('display','var'))
        display = false;
    end
    [sfx, info] = myShapeFromX(param_string, labeldata, display);
end
