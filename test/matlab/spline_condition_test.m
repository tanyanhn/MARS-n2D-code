%% 周期样条条件数实验：模块化版
% 目的：对比不同几何形状下，扰动距离 r_tiny 对条件数的影响

clear; clc; close all;

%% 1. 实验配置
N = 100;                        % 基础点数
r_start = 0.1;                  % r_tiny 起始值
r_end = 0.00001;                % r_tiny 结束值
rate = 0.9;                     % 衰减率

% --- [关键选项] 选择曲线形状 ---
% 选项: 'circle' 或 'sharp_bump'
curve_type = 'sharpBump'; 
% curve_type = 'circle'; 

% 尖锐参数: [基础半径, 凸起高度, 尖锐度]
% 尖锐度建议: 5(圆润) ~ 50(极尖)
bump_params = [1.0, 1.0, 10000]; 

%% 2. 初始化
r_tiny_vals = [r_start];
cond_results = [];

% 基础参数空间分布 (0 到 2*pi)
delta_theta = 2*pi / N;
theta_base = linspace(0, 2*pi - delta_theta, N);

fprintf('开始计算... N=%d, 形状: %s\n', N, curve_type);

%% 3. 循环计算
while r_tiny_vals(end) > r_end
    % 获取当前 r
    r = r_tiny_vals(end);
    
    % --- 几何参数构建 ---
    % 在 theta=0 附近插入扰动参数
    theta_tiny = delta_theta * r; 
    
    % 合并并排序参数 t
    thetas = sort([theta_base, theta_tiny]);
    
    % --- [调用模块] 生成几何点 ---
    % points 是一个 M x 2 的矩阵
%     points = get_curve_points(thetas, curve_type, bump_params);
    [points, ~, lookup] = get_adaptive_curve_points(N, r, curve_type, bump_params);
    
    % --- [调用模块] 构建矩阵 ---
    % 这里内部会自动计算物理弦长并构建矩阵
    A = build_periodic_matrix(points);
    
    % --- 计算条件数 ---
    cond_results = [cond_results, cond(full(A))]; 
    
    % 更新 r (准备下一次循环)
    r_tiny_vals = [r_tiny_vals, r * rate];
end

% 移除最后多生成的一个 r 值(循环跳出条件导致的)
r_tiny_vals = r_tiny_vals(1:end-1);

%% 4. 结果可视化

f = figure('Color', 'w', 'Position', [100, 100, 800, 500]);

% 子图1: 条件数变化曲线
% subplot(1, 2, 1);
plot(-log10(r_tiny_vals), cond_results, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
xlabel('-log_{10}(r_{tiny})', 'FontSize', 12);
ylabel('Condition Number \kappa(A)', 'FontSize', 12);
title(['Condition Number Growth (' curve_type ')'], 'FontSize', 14);
grid on;
addpath export_fig-3.54\
export_fig(f, "condNum.png");

% 子图2: 几何形状展示 (画出最后一个极近点的情况)
% subplot(1, 2, 2);
f= figure;
points = [points; points(1, :)];
plot(lookup(:,1), lookup(:,2), '-', 'LineWidth', 1); hold on;
plot(points(:,1), points(:,2), '.', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
% 标记那个极近的点 (通常是前两个点)
plot(points(1:2,1), points(1:2,2), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
axis equal; grid on;
title('Geometry Shape (Red: Tiny Interval)', 'FontSize', 14);
export_fig(f, "condNumShape.png");

% 简单的数值检查
fprintf('最小 r_tiny = %.5f 时, 条件数 = %.2e\n', min(r_tiny_vals), max(cond_results));


function [points, h_avg, lookup] = get_adaptive_curve_points(N, r_tiny_ratio, type, sharp_params)
    % GET_ADAPTIVE_CURVE_POINTS 生成基于物理弧长分布的点
    % 
    % 输入:
    %   N:            基础分段数 (生成 N 个均匀分布点 + 1 个扰动点)
    %   r_tiny_ratio: 扰动比例 (相对于平均物理步长 h_avg)
    %   type:         'circle' 或 'sharp_bump'
    %   sharp_params: 形状参数
    %
    % 输出:
    %   points:       生成的点坐标 (N+1) x 2
    %   h_avg:        计算出的平均物理步长 (用于验证)

    %% 1. 定义几何形状函数 (极坐标 r(theta))
    switch type
        case 'circle'
            geom_func = @(t) ones(size(t));
        case 'sharpBump'
            if nargin < 4, sharp_params = [1.0, 3.0, 20]; end
            R_base = sharp_params(1);
            H_bump = sharp_params(2);
            K_sharp = sharp_params(3);
            % Von Mises 分布形状
            geom_func = @(t) R_base + H_bump * exp(K_sharp * (cos(t) - 1));
        otherwise
            error('未知的曲线类型');
    end

    %% 2. 高密度采样 (计算累积弧长)
    % 为了保证插值精度，采样点数应远大于 N (例如 100倍)
    num_dense = max(10000, N * 100); 
    t_dense = linspace(0, 2*pi, num_dense + 1)'; % 闭合区间 [0, 2pi]
    t_dense(end) = []; % 去掉最后一个重复点 (0 和 2pi 重合)
    
    % 计算高密度点坐标
    r_dense = geom_func(t_dense);
    x_dense = r_dense .* cos(t_dense);
    y_dense = r_dense .* sin(t_dense);
    
    % 为了闭合计算弧长，临时加上起点到末尾
    x_calc = [x_dense; x_dense(1)];
    y_calc = [y_dense; y_dense(1)];
    
    % 计算相邻线段长度
    dists = sqrt(diff(x_calc).^2 + diff(y_calc).^2);
    
    % 累积弧长 (从 0 开始)
    cum_arc = [0; cumsum(dists)];
    total_length = cum_arc(end);
    
    % 去掉最后那个为了闭合加的点，保留 num_dense 个对应的弧长
    % 用于插值: x_dense 对应 cum_arc(1:end-1)
    s_dense = cum_arc(1:end-1); 
    
    % 处理插值的周期性边界问题：
    % 为了 interp1 能处理跨越 2pi 的情况，我们在数据末尾拼接一段
    s_lookup = [s_dense; total_length + s_dense(1:10)]; % 延长一小段
    x_lookup = [x_dense; x_dense(1:10)];
    y_lookup = [y_dense; y_dense(1:10)];

    %% 3. 确定目标采样位置 (基于物理距离)
    
    % 平均物理步长
    h_avg = total_length / N;
    
    % 生成 N 个均匀位置: 0, h, 2h, ..., (N-1)h
    s_uniform = (0:N-1)' * h_avg;
    
    % 生成扰动点位置: 在第1个点 (s=0) 之后 r_tiny * h_avg 处
    % 也可以选择在峰值最尖锐的地方(通常 t=0 处即 s=0 处)加扰动
    s_tiny = r_tiny_ratio * h_avg;
    
    % 合并所有目标弧长位置
    s_target = [s_uniform; s_tiny];
    
    % 排序 (为了画图好看，虽然矩阵构建不依赖顺序，但排序是个好习惯)
    s_target = sort(s_target);
    
    %% 4. 插值映射回坐标系
    % 使用 linear 插值即可，因为 dense 足够密
    x_out = interp1(s_lookup, x_lookup, s_target, 'linear');
    y_out = interp1(s_lookup, y_lookup, s_target, 'linear');
    lookup = [x_lookup, y_lookup];
    
    points = [x_out, y_out];
end

function A = build_periodic_matrix(points)
    % BUILD_PERIODIC_MATRIX 根据点集构建周期性三弯矩矩阵
    % 注意：这里使用弦长(物理距离)作为步长 h，这样几何形状才会影响条件数
    
    M = size(points, 1);
    
    % 1. 计算物理步长 h (欧几里得距离)
    % h(i) = dist(P_i, P_{i+1})
    diffs = diff(points);
    h_steps = sqrt(sum(diffs.^2, 2));
    
    % 处理闭合段: dist(P_M, P_1)
    d_close = sqrt(sum((points(1,:) - points(end,:)).^2));
    h = [h_steps; d_close]; 
    
    % 2. 构建稀疏矩阵
    A = sparse(M, M);
    
    for i = 1:M
        % 周期性索引处理
        idx_prev = mod(i - 2, M) + 1; 
        idx_curr = i;
        idx_next = mod(i, M) + 1;     
        
        % 获取相邻步长
        h_prev_val = h(idx_prev);
        h_curr_val = h(idx_curr); % 注意：根据三弯矩定义，这里对应的是 h_i
        
        sum_h = h_prev_val + h_curr_val;
        
        % 样条权重参数
        mu     = h_prev_val / sum_h;
        lambda = h_curr_val / sum_h;
        
        % 填充矩阵
        A(i, idx_prev) = mu;
        A(i, idx_curr) = 2;
        A(i, idx_next) = lambda;
    end
end