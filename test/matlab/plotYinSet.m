function plotYinSet(sf, lineColor, fillColor, lineWidth, offset)
    if numel(sf) == 0
        return;
    end
    minesVolume = -1e-8;
    
    if nargin < 5
        offset = [0 0];
    end
    if nargin < 4
        lineWidth = 1;
    end
    if nargin < 2
        fillColor = 'y';
    end
    if nargin < 2
        lineColor = '';
    end
    filled = numel(lineColor) == 0;
    
    fillcolor = cell(numel(sf), 1);
    linestyle = cell(numel(sf), 1);
    fillAlpha = cell(numel(sf), 1);
    
    % bouding box and orientation
    amax = 0;
    for c=1:numel(sf)
        a(c) = signedArea(linearizeSpline(sf{c})');
        if abs(a(c)) > abs(amax), amax = a(c); end
        if a(c) > 0
            fillcolor{c} = fillColor;
            linestyle{c} = '-';
            fillAlpha{c} = 0.7;
        elseif a(c) < minesVolume
            fillcolor{c} = 'w';
            linestyle{c} = ':';
            fillAlpha{c} = 1;
        end
    end
    
    % fill the background
    hold on;
    if filled && amax < minesVolume
        set(gca, 'Color', fillColor);
    end
    
    % plot the boundary or fill the interior
    for c=1:numel(sf)
        verts = linearizeSpline(sf{c}, 2)';
        verts = verts + offset;
        if filled && abs(a(c)) > -minesVolume
            fill(verts(:,1),verts(:,2), fillcolor{c}, ...
                'EdgeColor', 'k', 'FaceAlpha', fillAlpha{c});
        else
            plot(verts(:,1),verts(:,2), ...
                [lineColor,linestyle{c}], ...
                'LineWidth', lineWidth, ...
                'MarkerSize', 10);
        end
    end
    
    setAxisProperties();
    hold off;
end

function setAxisProperties()
    axis equal tight;
    ax = gca;
    ax.FontSize = 14;
    % ax.XColor = ax.Color;
    % ax.YColor = ax.Color;
    % ax.XTick = [];
    % ax.YTick = [];
end