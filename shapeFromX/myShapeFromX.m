function [sfx, info] = myShapeFromX(param_string, labeldata, display)
    if(display)
        param_string = [param_string, 'params.DO_DISPLAY = 1;'];
    end
    [results, state, data, params] = myGo(param_string, labeldata);
    sfx.im = im2double(labeldata.im);
    sfx.height = state.height;
    sfx.light = state.light;
    [logshading, ~, sfx.normals] = renderSH(sfx.height, sfx.light);
    sfx.shading = exp(logshading);
    sfx.albedo = exp(log(sfx.im)-logshading);
    sfx.mask = state.mask;
    sfx.renderFunc = @(height,light) renderSH(height, light);
    info = struct('results', results, 'state', state, 'data', data, 'params', params);
    if(display)
        myDisplaySFX(sfx);
    end
end
