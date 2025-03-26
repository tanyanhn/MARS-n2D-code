function curvePoints = bezierCurve(control_points, max_step)
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
t_list = [0; 1];
points = [bezierPoint(control_points, 0); 
          bezierPoint(control_points, 1)];

% 循环细分直到满足间距要求
while true
    % 计算所有相邻点间距
    distances = vecnorm(diff(points), 2, 2);
    
    % 找到需要细分的区间（从后往前处理）
    exceed_idx = find(distances > max_step);
    if isempty(exceed_idx), break; end
    
    % 对每个需要细分的区间插入中点
    for i = flip(exceed_idx)'
        t_mid = (t_list(i) + t_list(i+1)) / 2;
        mid_point = bezierPoint(control_points, t_mid);
        
        % 插入新点到列表中
        t_list = [t_list(1:i); t_mid; t_list(i+1:end)];
        points = [points(1:i,:); mid_point; points(i+1:end,:)];
    end
end

% 确保包含终点
if ~isequal(points(end,:), control_points(4,:))
    points(end+1,:) = control_points(4,:);
end

curvePoints = points;
end

% 三次贝塞尔曲线单点计算函数
function pt = bezierPoint(control_points, t)
    c0 = (1 - t)^3;
    c1 = 3*(1 - t)^2*t;
    c2 = 3*(1 - t)*t^2;
    c3 = t^3;
    pt = c0*control_points(1,:) + ...
         c1*control_points(2,:) + ...
         c2*control_points(3,:) + ...
         c3*control_points(4,:);
end