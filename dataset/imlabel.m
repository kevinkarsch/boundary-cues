% 1. Choose the type of labels from the radiobutton group
%
% 2. Click to add a point, type 'u' to undo, 's' to save a FINISHED polygon
% (interior boudary, convex corner, concave corner) and start a new one, 'q'
% to save and exit the current editing type. 'u', 's', 'q' are not used in 
% 'Horizon' mode (automatically exit after two clicks), 's' is not used for
% 'Exterior boundary' and 'Contact Pts'.
%
% 3. Delete a saved label object (polygon, line, contact point, etc) by: (a) 
% exit from current label mode, (b) click the 'Delete' btn, (c) choose the
% target object by clicking around any of its points (NN search), (d) type 
% 'y' to confirm deleting or other keys to cancel (then the procedure goto 
% (c) again). Type 'q' to exit 'Delete' mode.
%
% 4. When finished, exit from labeling or delete mode first, then click 'Done'

function labeldata = imlabel(im, labeldata)%imgid, VOCopts)
    hf=figure(1);
    set(hf, 'Position', [400 400 600 500]);

    h1 = uibuttongroup('parent', hf, 'Position',[0 0.9 1 0.1]);
    u1 = uicontrol('Style','Radio','String','Ext Bd','Position',[10 25 90 20],'parent',h1);
    u2 = uicontrol('Style','Radio','String','Interior Bd','pos',[100 25 90 20],'parent',h1);
    u3 = uicontrol('Style','Radio','String','Contact Pts','pos',[190 25 90 20],'parent',h1);
    u4 = uicontrol('Style','Radio','String','Horizon Ln','pos',[280 25 90 20],'parent',h1);
    u5 = uicontrol('Style','Radio','String','Parallel Ln','pos',[370 25 90 20],'parent',h1);
    u1 = uicontrol('Style','Radio','String','Ext Bd2','Position',[10 5 90 20],'parent',h1);
    u7 = uicontrol('Style','Radio','String','Convex Cnr','pos',[100 5 90 20],'parent',h1);
    u8 = uicontrol('Style','Radio','String','Concave Cnr','pos',[190 5 100 20],'parent',h1);
    
    hdel = uicontrol('Style','togglebutton','String','Delete','pos',[460 10 70 25],'parent',h1);
    hexit = uicontrol('Style','togglebutton','String','Done','pos',[540 10 50 25], 'parent', h1);   
    set(h1,'SelectionChangeFcn',@selcbk);
    set(h1,'SelectedObject',[]);  %select the first opt by default
    set(h1,'Visible','on');
    
    convexcnrclr = [0.33 0.1 0.55];
    concavecnrclr = [1 0.5 0.1];
    sharpextbdclr = [0 1 1];
    paralineclr = [0.58 0 0.83];
    markersize = 6;
    lwidth = 3;
    %initialize labeldata
    if ~exist('labeldata', 'var')
        labeldata.im = im;
        labeldata.extbd = cell(0);    %{[x y issharp]}
        labeldata.intbd = cell(0);    %{[x y]}
        labeldata.contact_pts = [];     %[x y]
        labeldata.horizonline = [];     %[x1 y1 x2 y2]
        labeldata.paraline = cell(0);        %{[x y]}
        labeldata.convexcnr = cell(0);
        labeldata.concavecnr = cell(0);
        labeldataexists = 0;
    else
        labeldataexists = 1;
        labeldata.im = im;
    end
    % backwards compatability with the previously stupid version that only
    % stored one external boundary
    if ~iscell(labeldata.extbd)
        labeldata.extbd = {labeldata.extbd};
    end
    %paralineclr = [rand 0.8 0.8]; %%close to cyan with randomness
    bd = [];
    plotscene(labeldata);
    
    %%check VOC groundtruth exterior boundary
    %%skip this step

    if(labeldataexists)
        axis off;
        return;
    end
    finished = 0;   %%value flips when pressing the Exit button.
    while(finished == 0)
        pause(1);
    end
    close(hf);
    
    function selcbk(source,eventdata)
        %plotscene(labeldata);
        opt = get(get(source,'SelectedObject'),'String');
        disp(opt);
        if strcmp(opt, 'Ext Bd') == 1 || strcmp(opt, 'Ext Bd2') == 1
            while(1)
                [x y k] = ginput(1);
                if y < -25
                    disp('exiting Ext Bd');
                    plotscene(labeldata);
                    break;
                elseif k == 'q'
                    if ~isempty(bd) 
                        labeldata.extbd{numel(labeldata.extbd)+1} = bd; 
                        disp('added an exterior boundary');
                        disp(bd); bd = [];
                    end
                    disp('Exiting Ext Bd');
                    break;
                elseif k == 's'
                    if ~isempty(bd) 
                        labeldata.extbd{numel(labeldata.extbd)+1} = bd; 
                        disp('added an exterior boundary');
                        disp(bd); bd = [];
                    end
                elseif k == 'u'
                    if size(bd,1) > 0 
                    %if size(labeldata.extbd,1) > 0
                        disp('delete the last contact point');
                        bd(end,:) = [];
                        plotscene(labeldata);
                    end
                else
                    if strcmp(opt, 'Ext Bd') == 1 && k ~= 2
                        bd = [bd; x y 0];
                        disp('added a smooth ext boundary point');
                    else
                        bd = [bd; x y 1];
                        disp('added a sharp ext boundary point');
                    end
                    disp([x y]);
                    drawextbd(bd,2);
                end
            end
            

        elseif strcmp(opt, 'Interior Bd') == 1
            pts = [];
            while(1)
                [x y k] = ginput(1);
                if y < -25
                    disp('exiting Interior Bd');
                    plotscene(labeldata);
                    break;
                elseif k == 'q'
                    if ~isempty(pts) 
                        labeldata.intbd{numel(labeldata.intbd)+1} = pts; 
                        disp('added an interior boundary');
                        disp(pts);
                    end
                    disp('exiting Interior Bd');
                    break;
                elseif k == 's'
                    if ~isempty(pts) 
                        labeldata.intbd{numel(labeldata.intbd)+1} = pts; 
                        disp('added an interior boundary');
                        disp(pts);
                    end
                    pts = [];
                elseif k == 'u' 
                    if size(pts, 1) > 0
                        disp('deleting the last interior boundary point');
                        pts(end,:) = [];
                        plotscene(labeldata);
                        if(~isempty(pts)) drawintbd(pts); end
                    end
                else
                    pts = [pts; [x y]];
                    drawintbd(pts);
                end
            end

        elseif strcmp(opt, 'Contact Pts') == 1
            while(1)
                [x y k] = ginput(1);
                if k== 'q' || y < -25
                    disp('exiting Contact Pts');
                    break;
                elseif k == 'u'
                    if size(labeldata.contact_pts,1) > 0
                        disp('delete the last contact point');
                        labeldata.contact_pts(end,:) = [];
                        plotscene(labeldata);
                    end
                else
                    labeldata.contact_pts = [labeldata.contact_pts; x y];
                    disp('added a contact point');
                    disp([x y]);
                    drawcontactpts([x y]);
                end
            end
            
        elseif strcmp(opt, 'Horizon Ln') == 1
            [x y] = ginput(1);
            if y > 0
                labeldata.horizonline = [x y];
                drawhorizonline(labeldata.horizonline);
                [x y] = ginput(1);
                
                labeldata.horizonline = [labeldata.horizonline [x y]];
                disp('add a line');
                disp(labeldata.horizonline);
            else
                disp('did nothing');
            end
            disp('exiting Horizon Line');
            plotscene(labeldata);
            
        elseif strcmp(opt, 'Parallel Ln') == 1
            pts = [];
            while(1)
                [x y k] = ginput(1);
                if y < -25
                    disp('exiting Parallel Line');
                    plotscene(labeldata);
                    break;
                elseif k == 'q'
                    if ~isempty(pts) 
                        labeldata.paraline{numel(labeldata.paraline)+1} = reshape(pts', 4, [])'; 
                        disp('added a set of parallel lines and exiting');
                        disp(reshape(pts', 4, [])');
                    end
                    disp('exiting Parallel Line');
                    break;
                elseif k == 's'
                    if ~isempty(pts) 
                        labeldata.paraline{numel(labeldata.paraline)+1} = reshape(pts', 4, [])'; 
                        disp('added a set of parallel lines');
                        disp(reshape(pts', 4, [])');
                    end
                    pts = [];
                    %paralineclr = [rand 0.8 0.8];
                    continue;
                elseif k == 'u'
                    if size(pts,1) > 0
                        pts(end,:) = [];
                        plotscene(labeldata);
                        drawparallellines(pts);
                    end
                else
                    pts = [pts; [x y]];
                    drawparallellines(pts);
                end
                
                [x y k] = ginput(1);
                if k == 'u'
                    if size(pts,1) > 0
                        pts(end,:) = [];
                        plotscene(labeldata);
                        drawparallellines(pts);
                    end
                else
                    pts = [pts; [x y]];
                end
                plotscene(labeldata);
                drawparallellines(pts);
            end
            plotscene(labeldata);

        elseif strcmp(opt, 'Convex Cnr')
            pts = [];
            while(1)
                [x y k] = ginput(1);
                if y < -25      %%25 = 500 (window height) * 0.05 (gap between two sub windows)
                    disp('exiting Convex Cnr');
                    plotscene(labeldata);
                    break;
                elseif k == 'q'
                    if ~isempty(pts) 
                        labeldata.convexcnr{numel(labeldata.convexcnr)+1} = pts; 
                        disp('added a convex corner');
                        disp(pts);
                    end
                    disp('exiting Convex Cnr');
                    break;
                elseif k == 's'
                    if ~isempty(pts) 
                        labeldata.convexcnr{numel(labeldata.convexcnr)+1} = pts; 
                        disp('added an convex boundary');
                        disp(pts);
                    end
                    pts = [];
                elseif k == 'u'
                    if size(pts, 1) > 0
                        disp('deleting the last convex corner point');
                        pts(end,:) = [];
                        plotscene(labeldata);
                        if(~isempty(pts)) drawconvexcnr(pts); end
                    end
                else
                    pts = [pts; [x y]];
                    drawconvexcnr(pts);
                end
            end
            
        elseif strcmp(opt, 'Concave Cnr')
             pts = [];
            while(1)
                [x y k] = ginput(1);
                if y < -25      %%25 = 500 (window height) * 0.05 (gap between two sub windows)
                    disp('exiting Convex Cnr');
                    plotscene(labeldata);
                    break;
                elseif k == 'q'
                    if ~isempty(pts) 
                        labeldata.concavecnr{numel(labeldata.concavecnr)+1} = pts; 
                        disp('added a concave corner');
                        disp(pts);
                    end
                    disp('exiting Concave Cnr');
                    break;
                elseif k == 's'
                    if ~isempty(pts) 
                        labeldata.concavecnr{numel(labeldata.concavecnr)+1} = pts; 
                        disp('added an concave corner');
                        disp(pts);
                    end
                    pts = [];
                elseif k == 'u'
                    if size(pts, 1) > 0
                        disp('deleting the last concave corner point');
                        pts(end,:) = [];
                        plotscene(labeldata);
                        if(~isempty(pts)) drawconcavecnr(pts); end
                    end
                else
                    pts = [pts; [x y]];
                    drawconcavecnr(pts);
                end
            end
            
        elseif strcmp(opt, 'Delete') == 1   %Delete elements
            while(1)
                [x y k ] = ginput(1);
                if k == 'q'|| y < -25
                    disp('exiting Delete');
                    break;
                else
                    [field ind] = deleteNN([x y]);
                    if ~isempty(field)
                        if strcmp(field, 'intbd') || strcmp(field, 'convexcnr')...
                                  || strcmp(field, 'concavecnr') || strcmp(field, 'paraline')
                            selected = labeldata.(field){ind};
                        elseif strcmp(field, 'extbd')
                            selected = labeldata.(field){ind}(:,1:2);
                        elseif strcmp(field, 'contact_pts')
                            selected = labeldata.(field)(ind, :);
                        elseif strcmp(field, 'horizonline')
                            selected = labeldata.(field);
                        end
                        selected = reshape(selected', 2, [])';
                        for i = 1:size(selected,1)
                            plot(selected(i,1), selected(i,2), 'wx', 'MarkerSize', markersize*1.5, 'LineWidth', lwidth*1.5);
                        end
                        disp('Press [y] to delete, other keys to skip.');
                        [x y k] = ginput(1);
                        if k == 'y'
                            if strcmp(field, 'extbd') || strcmp(field, 'intbd') || strcmp(field, 'convexcnr')...
                                  || strcmp(field, 'concavecnr') || strcmp(field, 'paraline')
                                labeldata.(field)(ind) = [];
                            elseif strcmp(field, 'contact_pts')
                                labeldata.(field)(ind, :) = [];
                            elseif strcmp(field, 'horizonline')
                                labeldata.(field) = [];
                            end
                        end
                        plotscene(labeldata);
                    else
                        disp('nothing selected.');
                    end
                end
            end

        elseif strcmp(opt, 'Done') == 1
            finished = 1;
        end
    end

%     function drawextbd(pts)
%         plot(pts(1,1), pts(1,2), 'y*', 'MarkerSize', markersize, 'LineWidth', lwidth);
%         plot(pts(:,1), pts(:,2), 'r', 'LineWidth', lwidth);
%     end
    function drawextbd(pts, lw)
        if ~exist('lw', 'var')
            lw = lwidth;
        end
        if ~isempty(pts)
            %plot(pts(1,1), pts(1,2), 'y*', 'MarkerSize', markersize, 'LineWidth', lw);
            for i = 1:size(pts,1)-1
                x = [pts(i,1) pts(i+1,1)];
                y = [pts(i,2) pts(i+1,2)];
                if pts(i+1,3) > 0.5
                    plot(x, y, 'color', sharpextbdclr, 'LineWidth', lw);
                else
                    plot(x, y, 'r', 'LineWidth', lw);
                end
            end
        end
    end
    function drawintbd(pts)
        %plot(pts(1,1), pts(1,2), 'y*', 'MarkerSize', markersize, 'LineWidth', lwidth);
        plot(pts(:,1), pts(:,2), 'g', 'LineWidth', lwidth);
    end
    function drawconvexcnr(pts)
        %plot(pts(1,1), pts(1,2), 'y*', 'MarkerSize', markersize, 'LineWidth', lwidth);
        plot(pts(:,1), pts(:,2), 'color', convexcnrclr, 'LineWidth', lwidth);
    end
    function drawconcavecnr(pts)
        %plot(pts(1,1), pts(1,2), 'y*', 'MarkerSize', markersize, 'LineWidth', lwidth);
        plot(pts(:,1), pts(:,2), 'color', concavecnrclr, 'LineWidth', lwidth);
    end
    function drawcontactpts(pt)
        plot(pt(1), pt(2), 'bs', 'MarkerSize', markersize, 'LineWidth', lwidth);
    end
    function drawhorizonline(pts)
        if numel(pts) < 4
            plot(pts(1), pts(2), 'y*', 'MarkerSize', markersize, 'LineWidth', lwidth); 
        else
            xs = [pts(1) pts(3)];
            ys = [pts(2) pts(4)];
            plot(xs, ys, 'm', 'LineWidth', lwidth+1);
        end
    end
    function drawparallellines(pts) 
        if mod(size(pts,1),2) == 1
            plot(pts(end,1), pts(end,2), 'y*', 'MarkerSize', markersize, 'LineWidth', lwidth);
        end
        for i = 1:2:size(pts,1)
            if i+1 <= size(pts,1)
                xs = [pts(i,1) pts(i+1,1)];
                ys = [pts(i,2) pts(i+1,2)];
                plot(xs, ys, 'color', paralineclr, 'LineWidth', lwidth);
            end
        end
    end
    function plotscene(labeldata)
        h2 = axes('Parent',hf,'Position',[0.05, 0.05, 0.91, 0.8]);
        imagesc(labeldata.im);
        hold on;
        axis image;
        
        %horizon line
        if ~isempty(labeldata.horizonline)
            drawhorizonline(labeldata.horizonline); 
        end
        
        %parallel lines
        for i = 1:numel(labeldata.paraline)
            lines = labeldata.paraline{i};
            for j = 1:size(lines,1)
                xs = [lines(j,1) lines(j,3)];
                ys = [lines(j,2) lines(j,4)];
                plot(xs, ys, 'color', paralineclr, 'LineWidth', lwidth);
            end
        end
        
        %convex corner
        for i = 1:numel(labeldata.convexcnr)
            drawconvexcnr(labeldata.convexcnr{i});
        end
        
        %concave corner
        for i = 1:numel(labeldata.concavecnr)
            drawconcavecnr(labeldata.concavecnr{i});
        end
       
        for i = 1:numel(labeldata.extbd)
            drawextbd(labeldata.extbd{i});
        end
        drawextbd(bd, 2);
%         x = [labeldata.extbd(end,1) labeldata.extbd(1,1)];
%         y = [labeldata.extbd(end,2) labeldata.extbd(1,2)];
%         if labeldata.extbd(end,3) > 0.5 && labeldata.extbd(1,3) > 0.5
%             plot(x, y, 'color', sharpextbdclr, 'LineWidth', lwidth);
%         else
%             plot(x, y, 'r', 'LineWidth', lwidth);
%         end
        
        %int boundary
        for i = 1:numel(labeldata.intbd)
            drawintbd(labeldata.intbd{i});
        end
        
        %contact points
        for i = 1:size(labeldata.contact_pts,1)
            drawcontactpts(labeldata.contact_pts(i,:));
        end
        
        
    end
    function [field, ind] = deleteNN(pt) 
        field = ''; ind = NaN;
        mindist = 1000;
        
        for i = 1:numel(labeldata.extbd)
            pts = labeldata.extbd{i}(:,1:2);
            for j = 1:size(pts,1)
                if norm(pt - pts(j,:)) < mindist && norm(pt - pts(j,:)) < 10
                    mindist = norm(pt - pts(j,:));
                    field = 'extbd';
                    ind = i;
                end
            end
        end
        
        for i = 1:numel(labeldata.intbd)
            pts = labeldata.intbd{i};
            for j = 1:size(pts,1)
                if norm(pt - pts(j,:)) < mindist && norm(pt - pts(j,:)) < 10
                    mindist = norm(pt - pts(j,:));
                    field = 'intbd';
                    ind = i;
                end
            end
        end
        
        for i = 1:numel(labeldata.convexcnr)
            pts = labeldata.convexcnr{i};
            for j = 1:size(pts,1)
                if norm(pt - pts(j,:)) < mindist && norm(pt - pts(j,:)) < 10
                    mindist = norm(pt - pts(j,:));
                    field = 'convexcnr';
                    ind = i;
                end
            end
        end
        
        for i = 1:numel(labeldata.concavecnr)
            pts = labeldata.concavecnr{i};
            for j = 1:size(pts,1)
                if norm(pt - pts(j,:)) < mindist && norm(pt - pts(j,:)) < 10
                    mindist = norm(pt - pts(j,:));
                    field = 'concavecnr';
                    ind = i;
                end
            end
        end
        
        pts = labeldata.contact_pts;
        if ~isempty(pts)
            for j = 1:size(pts,1)
                if norm(pt - pts(j,:)) < mindist && norm(pt - pts(j,:)) < 10
                    mindist = norm(pt - pts(j,:));
                    field = 'contact_pts';
                    ind = j;
                end
            end
        end
        
        pts = labeldata.horizonline;
        if ~isempty(pts)
            if norm(pt - pts(1:2)) < mindist && norm(pt - pts(1:2)) < 10
                mindist = norm(pt - pts(1:2));
                field = 'horizonline';
                ind = NaN;
            elseif norm(pt - pts(3:4)) < mindist && norm(pt - pts(3:4)) < 10
                mindist = norm(pt - pts(3:4));
                field = 'horizonline';
                ind = NaN;
            end
        end
        
        for i = 1:numel(labeldata.paraline)
            pts = labeldata.paraline{i};
            for j = 1:size(pts, 1)
                if norm(pt - pts(j,1:2)) < mindist && norm(pt - pts(j, 1:2)) < 10
                    mindist = norm(pt - pts(j,1:2));
                    field = 'paraline';
                    ind = i;
                elseif norm(pt - pts(j, 3:4)) < mindist && norm(pt - pts(j, 3:4)) < 10
                    mindist = norm(pt - pts(j, 3:4));
                    field = 'paraline';
                    ind = i;
                end
            end
        end
    end   
end