function t_opt = inverseBezier(control_points, target_point, t_init)
if nargin < 3
    % 采样参数t以找到初始候选点
    num_samples = 100;
    t_samples = linspace(0, 1, num_samples);

    % 计算每个样本点的距离平方
    dist_sq = zeros(1, num_samples);
    for i = 1:num_samples
        pt = bezierPoint(control_points, t_samples(i));
        dist_sq(i) = sum((pt - target_point).^2);
    end

    % 找到最小距离的样本点
    [~, min_index] = min(dist_sq);
    t_init = t_samples(min_index);
end
    
    % 定义目标函数（距离平方）
    obj = @(t) sum((bezierPoint(control_points, t) - target_point').^2);
    
    % 设置优化选项
    options = optimset('TolX', 1e-40, 'Display', 'on');
    
    % 在初始点附近定义搜索窗口
    window = 0.1;
    lower = max(0, t_init - window);
    upper = min(1, t_init + window);
    
    % 执行局部优化
    t_local = fminbnd(obj, lower, upper, options);
    
    % 比较并返回最优结果
    if obj(t_local) < obj(t_init)
        t_opt = t_local;
    else
        t_opt = t_init;
    end
end
