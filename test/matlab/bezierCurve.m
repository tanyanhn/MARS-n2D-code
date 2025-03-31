function [points t_list] = bezierCurve(control_points, max_step, t_list)
% 自适应生成贝塞尔曲线点列
% 输入：
%   control_points - 4x2矩阵，每行代表控制点的[x, y]坐标
%   max_step       - 相邻点最大允许间距
% 输出：
%   curvePoints    - Nx2矩阵，满足间距要求的曲线点列

% 参数校验
if size(control_points,1) ~= 4
    error('控制点必须为4行');
end

% 初始化包含起点和终点的列表
if nargin < 3
    t_list = [0 1];
end
points = bezierPoint(control_points, t_list);

% 循环细分直到满足间距要求
while true
    % 计算所有相邻点间距
    distances = vecnorm(diff(points), 2, 2);
    
    % 找到需要细分的区间（从后往前处理）
    exceed_idx = find(distances > max_step);
    if isempty(exceed_idx), break; end
    
    % 对每个需要细分的区间插入中点
    new_t_list = [];
    new_points = [];
    for i = 1:length(t_list) - 1
        if (distances(i) >max_step)
            num = ceil(distances(i) / max_step) + 1;
            t_mid = linspace(t_list(i), t_list(i + 1), num + 1);
            t_mid = t_mid(1:end-1)';
            mid_point = bezierPoint(control_points, t_mid);
        else
            t_mid = t_list(i);
            mid_point = points(i, :);
        end
        
        % 插入新点到列表中
        new_t_list = [new_t_list; t_mid];
        new_points = [new_points; mid_point];
    end
    new_t_list = [new_t_list; t_list(end)];
    new_points = [new_points; points(end, :)];
    t_list = new_t_list;
    points = new_points;
end

points = points';

end
