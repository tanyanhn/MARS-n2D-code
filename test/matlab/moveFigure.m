function panAxis(direction, stepRatio)
    % panAxis 平移当前坐标轴的可视范围
    % 输入参数:
    %   direction : 移动方向 ('up','down','left','right')
    %   stepRatio : 移动步长比例（可选，默认1.0）
    %
    % 示例:
    %   panAxis('right')       % 向右移动一个视图宽度
    %   panAxis('up', 0.5)     % 向上移动半个视图高度
    
    if nargin < 2
        stepRatio = 1.0; % 默认移动整个视图范围
    end
    subRatio = stepRatio * 0.23;

    ax = gca;
    
    % 获取当前坐标范围
    xl = xlim(ax);
    yl = ylim(ax);
    newMinX = xl(1);
    newMaxX = xl(2);
    newMinY = yl(1);
    newMaxY = yl(2);
    
    % 计算视图尺寸
    viewWidth = diff(xl);
    viewHeight = diff(yl);
    
    
    % 根据方向调整范围
    switch lower(direction)
        case 'up'
            newMinY = yl(1) + viewHeight*stepRatio;
            newMaxY = yl(2) + viewHeight*stepRatio;
            newMinX = xl(1) + viewWidth*subRatio;
            newMaxX = xl(2) + viewWidth*subRatio;
            % newMaxY = min(newMaxY, dataYL(2)); % 限制不超过数据上限
            
        case 'down'
            newMinY = yl(1) - viewHeight*stepRatio;
            newMaxY = yl(2) - viewHeight*stepRatio;
            newMinX = xl(1) - viewWidth*subRatio;
            newMaxX = xl(2) - viewWidth*subRatio;
            % newMinY = max(newMinY, dataYL(1)); % 限制不低于数据下限
            
        case 'right'
            newMinX = xl(1) + viewWidth*stepRatio;
            newMaxX = xl(2) + viewWidth*stepRatio;
            % newMaxX = min(newMaxX, dataXL(2)); % 限制不超过数据上限
            
        case 'left'
            newMinX = xl(1) - viewWidth*stepRatio;
            newMaxX = xl(2) - viewWidth*stepRatio;
            % newMinX = max(newMinX, dataXL(1)); % 限制不低于数据下限
            
        otherwise
            error('无效方向，可用选项: up, down, left, right');
    end

    % 应用新范围
    % if contains(lower(direction), {'up','down'})
        ylim(ax, [newMinY newMaxY])
    % else
        xlim(ax, [newMinX newMaxX])
    % end

end
