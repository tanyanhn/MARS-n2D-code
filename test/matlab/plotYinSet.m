function plotYinSet(sf, lineColor, lineWidth, offset)
    if numel(sf) == 0
        return;
    end
    
    if nargin < 4
        offset = [0 0];
    end
    if nargin < 3
        lineWidth = 1;
    end
    if nargin < 2
        lineColor = '';
    end
    filled = numel(lineColor) == 0;
    
    fillcolor = cell(numel(sf), 1);
    linestyle = cell(numel(sf), 1);
    
    % bouding box and orientation
    amax = 0;
    for c=1:numel(sf)
        a = signedArea(linearizeSpline(sf{c})');
        if abs(a) > abs(amax), amax = a; end
        if a > 0
            fillcolor{c} = 'y';
            linestyle{c} = '-';
        else
            fillcolor{c} = 'w';
            linestyle{c} = ':';
        end
    end
    
    % fill the background
    hold on;
    if filled && amax < 0
        set(gca, 'Color', 'y');
    end
    
    % plot the boundary or fill the interior
    for c=1:numel(sf)
        verts = linearizeSpline(sf{c})';
        verts = verts + offset;
        if filled
            fill(verts(:,1),verts(:,2), fillcolor{c}, ...
                'EdgeColor', 'none');
        else
            plot(verts(:,1),verts(:,2), ...
                [lineColor,linestyle{c}], ...
                'LineWidth', lineWidth, ...
                'MarkerSize', 10);
        end
    end
    
    %setAxisProperties();
    hold off;
end

% function setAxisProperties()
%     axis equal tight;
%     ax = gca;
%     ax.XColor = ax.Color;
%     ax.YColor = ax.Color;
%     ax.XTick = [];
%     ax.YTick = [];
% end