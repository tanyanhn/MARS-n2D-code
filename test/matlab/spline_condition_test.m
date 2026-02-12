%% 周期样条条件数实验：N 列表对比版
% 目的：输入一个 N 列表，为每种形状生成一张图，对比不同 N 下的条件数变化

clear; clc; close all;

%% 1. 实验配置
% --- [关键修改] 输入 N 的列表 ---
% 建议 N 不要太大(如 >2000)，否则 cond() 计算会非常耗时
N_list = [50, 100, 200, 500, 1000, 2000]; 

r_exponents = linspace(-1, -5, 50); % 指数范围
r_tiny_vals = 10.^r_exponents;      % r_tiny 实际值

% --- 形状列表 ---
shape_list = {'circle', 'sharpBump'}; 
type_list = {'not-a-knot', 'periodic'};
display_names = {'Unit Circle', 'Sharp Bump'};

% 尖锐参数: [基础半径, 凸起高度, 尖锐度]
bump_params = [1.0, 1.0, 10000]; 

fprintf('开始计算... N_list = [%s]\n', num2str(N_list));

%% 2. 循环计算

% --- 第一层循环：遍历形状 (每个形状一张图) ---
for s_idx = 2:length(type_list)
%     current_type = shape_list{s_idx};
    current_type = shape_list{1};
    bc_type = type_list{s_idx};
    shape_name = display_names{s_idx};
    
    % 为当前形状新建一个图形窗口
    f = figure('Color', 'w', 'Position', [100 + (s_idx-1)*50, 100, 800, 500]);
    hold on;
    
    % 获取颜色列表 (让不同的 N 显示不同颜色)
    colors = lines(length(N_list));
    
    fprintf('>>> 正在处理形状: %s \n', current_type);
    fprintf('>>> 正在处理条件: %s \n', bc_type);
    
    % --- 第二层循环：遍历 N 列表 (每条线对应一个 N) ---
    for n_idx = 1:length(N_list)
        N = N_list(n_idx);
        current_color = colors(n_idx, :);
        
        cond_results = zeros(size(r_tiny_vals));
        
        % --- 第三层循环：计算该 N 下的条件数曲线 ---
        for k = 1:length(r_tiny_vals)
            r = r_tiny_vals(k);
            
            % [调用模块] 生成自适应分布的点
            % 注意：r 是物理距离比例
            [points, ~, lookup] = get_adaptive_curve_points(N, r, current_type, bump_params);
            
            % [调用模块] 构建矩阵
            A = build_spline_matrix(points, bc_type);
            fullA = full(A);
            
            % 计算条件数
            cond_results(k) = cond(fullA); 
        end
        
        % --- 绘图 (单条曲线) ---
        plot(-log10(r_tiny_vals), cond_results, '.-', ...
            'Color', current_color, ...
            'LineWidth', 1.5, ...
            'MarkerSize', 10, ...
            'DisplayName', ['N = ', num2str(N)]);
            
        fprintf('    完成 N = %d\n', N);
    end
    
    % --- 单张图的美化与保存 ---
%     xlabel('-log_{10}(r_{tiny})', 'FontSize', 12);
%     ylabel('Condition Number \kappa(A)', 'FontSize', 12);
%     title(['Condition Number Growth: ' shape_name], 'FontSize', 14);
    grid on;
    legend('show', 'Location', 'east', 'FontSize', 14);
    
    % 自动保存图片
    file_name = sprintf('CondNum_%s_%s.png', current_type, bc_type);
    try
        addpath export_fig-3.54\ % 如果你有这个库
        ax = gca;
        ax.FontSize = 14;
        export_fig(f, file_name);
        fprintf('    图片已保存: %s\n', file_name);
    catch
        saveas(f, file_name); % 使用 MATLAB 自带保存作为备选
        fprintf('    图片已保存 (saveas): %s\n', file_name);
    end
    fprintf('\n');
end

fprintf('所有计算完成。\n');

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
    s_uniform = (0:N)' * h_avg;
    
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

function A = build_spline_matrix(points, bc_type)
    % BUILD_SPLINE_MATRIX 构建三次样条插值系数矩阵
    % 输入:
    %   points:  Nx2 点集坐标
    %   bc_type: 'periodic' (默认, 封闭圆环) 或 'not-a-knot' (开放曲线)
    %
    % 输出:
    %   A:       系数矩阵 (对于 periodic 为循环矩阵，对于 not-a-knot 为三对角矩阵)
    
    if nargin < 2
        bc_type = 'periodic';
    end

    M = size(points, 1);
    
    % --- 1. 计算物理步长 h ---
    % h_steps(i) = dist(P_i, P_{i+1})
    % 共 M-1 个段
    diffs = diff(points);
    h_steps = sqrt(sum(diffs.^2, 2));
    
    h = h_steps;
    if strcmp(bc_type, 'periodic')
        M = M - 1;
    else
    end
    
    % --- 2. 构建稀疏矩阵 ---
    A = sparse(M, M);
    
    if strcmp(bc_type, 'periodic')
        % ==============================
        % Mode A: Periodic (周期性闭合)
        % ==============================
        for i = 1:M
            idx_prev = mod(i - 2, M) + 1; 
            idx_curr = i;
            idx_next = mod(i, M) + 1;     
            
            % 周期性下，h 数组长度为 M，索引循环使用
            h_prev_val = h(idx_prev);
            h_curr_val = h(idx_curr); 
            
            sum_h = h_prev_val + h_curr_val;
            
            % 三弯矩方程系数
            mu     = h_prev_val / sum_h;
            lambda = h_curr_val / sum_h;
            
            A(i, idx_prev) = mu;
            A(i, idx_curr) = 2;
            A(i, idx_next) = lambda;
        end
        
    elseif strcmp(bc_type, 'not-a-knot')
        % ==============================
        % Mode B: Not-a-Knot (开放端点)
        % ==============================
        
        % 1. 内部节点 (2 到 M-1)：使用标准三弯矩方程
        for i = 2:M-1
            % 对于开放曲线:
            % h 索引: h(1)是P1->P2, h(i-1)是Pi-1->Pi, h(i)是Pi->Pi+1
            h_prev_val = h(i-1);
            h_curr_val = h(i);
            
            sum_h = h_prev_val + h_curr_val;
            
            mu     = h_prev_val / sum_h;
            lambda = h_curr_val / sum_h;
            
            A(i, i-1) = mu;
            A(i, i)   = 2;
            A(i, i+1) = lambda;
        end
        
        % 2. 边界条件处理 (基于三阶导数连续推导)
        % 方程: h_2 * M_1 - (h_1 + h_2) * M_2 + h_1 * M_3 = 0
        % 这保证了第一段和第二段的三次系数相同
        
        % --- 第一个点 (Row 1) ---
        h1 = h(1);
        h2 = h(2);
        A(1, 1) = h2 / (h1 + h2);
        A(1, 2) = -1;
        A(1, 3) = h1/ (h1 + h2);
%         A(1, 1) = h2;
%         A(1, 2) = -(h1 + h2);
%         A(1, 3) = h1;
        
        % --- 最后一个点 (Row M) ---
        % 方程对称: h_{end} * M_{end-2} - (...) * M_{end-1} + h_{end-1} * M_{end} = 0
        h_end_1 = h(end-1); % h_{M-2}
        h_end   = h(end);   % h_{M-1}
        
        A(M, M-2) = h_end/ (h_end_1 + h_end);
        A(M, M-1) = -1;
        A(M, M)   = h_end_1 / (h_end_1 + h_end);
        
    else
        error('未知的边界条件类型: %s', bc_type);
    end
end