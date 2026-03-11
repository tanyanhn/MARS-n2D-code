function multi_sample_spline_demo()
    % multi_sample_spline_demo 
    % 在同一张图上，对比基于不同样本点序列 (x_sequences) 的样条插值效果

    clc; clear; close all;

    % ===========================
    % 1. 参数配置区
    % ===========================
    
    % 1.1 定义目标函数 f(x)
    % 这里的例子：Runge 现象函数变体，或者周期函数
    % 参数设置：sigma 越小，0.1 处的二阶导数越大
    sigma = 0.1; 
    x0 = 0.2;

    % 构造 C4（实际上是 C_inf）连续函数
    % 第一项是高斯尖峰，第二项是平滑的背景波动（类似你提供的格式）
%     f = @(x) exp(-(x - x0).^2 / (2 * sigma^2));
%     f = @(x) sin(2 * pi * x);
    
    % 1.2 定义多组不同的样本点序列 (使用 cell 数组)
    % 可以在这里随意添加不同的采样策略
    x_sequences = {};
    error = [];
    rh = 0.01;
    n = 10;
    
    ngrid = 5;
    for ig = 1:ngrid
        t = linspace(0, 1, n+1);
        t = [t(1), t(2) * (1 - rh), t(2), t(2) * (1 + rh), t(3:end)];
        x_sequences{end+1} = t;
        n = n * 2;
    end
    

    % 1.3 选择插值模式
    % 可选: 'not-a-knot' (默认, 适用于非周期) 或 'periodic' (适用于周期)
    interp_method = 'not-a-knot'; 
%     interp_method = 'periodic'; 

    % ===========================
    % 2. 计算与绘图
    % ===========================
    
    % 准备画图
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on; grid on; box on;
    
    % 2.1 画由细网格生成的“真值”曲线 (作为背景水印)
%     x_fine = linspace(0, 1, 1000);
%     y_true = f(x_fine);
    
    % 技巧：使用 patch 代替 plot 来支持透明度 (EdgeAlpha)
    % [x_fine, NaN] 的写法是为了防止 patch 自动将首尾封闭填充
%     p_true = patch([x_fine, NaN], [y_true, NaN], 'k'); 
    
%     set(p_true, ...
%         'EdgeColor', [0.1 0.4 0.8], ... % 颜色：深灰色 (或你可以选 [0 0.4 0.8] 蓝色)
%         'FaceColor', 'none', ...        % 填充：无
%         'LineWidth', 6, ...             % 宽度：6 (非常粗，像背景马克笔)
%         'EdgeAlpha', 0.15);             % 透明度：0.15 (15% 不透明，很淡，网格线能透出来)
    
    % 获取颜色映射，为每组样本分配不同颜色
    colors = lines(length(x_sequences));
    legend_str = {'True curve f(x)'};
%     legend_str = {};
    
    % 2.2 循环处理每一组样本序列
    for i = 1:length(x_sequences)
        % 获取当前样本点
        x_sample = x_sequences{i};
        rb = (x_sample(3) - x_sample(2)) / (x_sample(2) - x_sample(1));
        x0 = x0 / 2;
        f = @(x) exp(-(x - x0).^2 / (2 * sigma^2));
        x_fine = linspace(0, 1, 1000);
        y_true = f(x_fine);
        
        % 如果是 Periodic 模式，且输入包含 0 和 1，
        % 求解器通常需要去重（因为 f(0)=f(1)），或者求解器内部处理。
        % 这里为了简单，如果用 periodic，我们确保传入 my_cubic_spline 的点包含边界
        % 具体的去重逻辑在 my_cubic_spline 内部处理或保持原样即可。
        
        y_sample = f(x_sample);
        
        % === 调用插值函数 ===
        [y_interp, cond] = my_cubic_spline(x_sample, y_sample, x_fine, interp_method);
        
        % 绘图 - 插值曲线
        plot(x_fine, y_interp, '--', 'LineWidth', 1.5, 'Color', colors(i,:));

        error(end+1) = max(abs(y_interp - y_true));
        
        % 绘图 - 样本点 (实心圆点)
        if (i == 1)
            plot(x_sample, y_sample, 'o', ...
                'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'w', ...
                'MarkerSize', 7, ...
                'LineWidth', 1);
        else
%             plot(x_sample(2), y_sample(2), 'o', ...
%                     'MarkerFaceColor', colors(i,:), ...
%                     'MarkerEdgeColor', 'w', ...
%                     'MarkerSize', 7, ...
%                     'LineWidth', 1);
%             text(x_sample(2), y_sample(2), '  Extra Point', ... % 前面加空格微调间距
%             'FontSize', 13, ...
%             'FontWeight', 'bold', ...
%             'Color', colors(i,:), ...
%             'HorizontalAlignment', 'left', ... % 水平靠左
%             'VerticalAlignment', 'middle');    % 垂直居中
        end
        
        % 更新图例文字
        legend_str{end+1} = sprintf('Spline with r_b=%.2g, cond=%.2g', rb, cond);
    end
    % ===========================
    % 3. 图表修饰
    % ===========================
%     title(['Cubic Spline Interpolation: ' interp_method], 'FontSize', 20);
%     xlabel('x'); ylabel('f(x)');
    
    % 只为“真值线”和“插值曲线”显示图例，忽略样本点的句柄
    % 这是一个常用的 MATLAB 技巧：获取所有子对象，选取特定几个显示
    children = get(gca, 'Children');
    % 由于 plot 顺序 (真值 -> (线,点) -> (线,点)... )，我们需要倒序或手动指定
    % 最简单的方法是上面构造 legend_str 时偷懒，或者如下只标注 Line
    
    % 重新简单的图例生成：
    % 这种方法只给线加图例，不给点加
%     h_legend = findobj(gca, 'Type', 'Line', 'LineStyle', '--');
%     h_true = findobj(gca, 'Type', 'Line', 'Color', [0.2 0.2 0.2]);
%     legend([p_true; flip(h_legend)], [legend_str(1:1:end)], 'FontSize', 14);
    
%     ylim([min(y_true)-0.5, max(y_true)+0.5]);
    axis tight 
    f = gca;
    f.FontSize = 14;
    export_fig(f, "Not-a-knot_illResult.png");
    % === 2. 计算收敛阶 ===
n = length(error);
orders = zeros(1, n); % 初始化收敛阶数组

% 向量化计算: log2( 上一次误差 / 当前误差 )
% 注意: 第一个误差没有“上一次误差”，所以保持为 0
orders(2:end) = log2(error(1:end-1) ./ error(2:end));

% === 3. 格式化输出 ===
fprintf('%-15s %-15s\n', '误差(Error)', '收敛阶(Order)');
fprintf('-------------------------------\n');

for i = 1:n
    % 输出格式: 科学计数法显示误差，浮点数显示收敛阶
    if i == 1
        % 第一个点没有收敛阶，显示 - 或 0
        fprintf('%-15.4e %-15s\n', error(i), '-');
    else
        fprintf('%-15.4e %-15.4f\n', error(i), orders(i));
    end
end
end

% ---------------------------------------------------------
% 核心插值函数 (包含 Not-a-knot 和 Periodic 实现)
% ---------------------------------------------------------
function [y_query, co] = my_cubic_spline(x, y, xq, type)
x = x(:); y = y(:);
n = length(x);
h = x(2:end) - x(1:end-1);

[A, rhs] = cubic_spline_matrix(h, type, y);
co = cond(A);

% 解出二阶导数 M
M = A \ rhs;

% 在 xq 处插值
y_query = zeros(size(xq));

% 扩展节点以便查找 (Wrapping)
x_ext = [x];
y_ext = [y];
M_ext = [M];
h_ext = h;

for k = 1:length(xq)
    val = xq(k);
    % 简单处理：将查询点映射回 [x(1), x(end)] 范围内(如果是纯周期函数)
    % 这里假设 xq 就在 [0, 1] 范围内

    % 查找区间
    idx = find(x_ext <= val, 1, 'last');
    if isempty(idx), idx = 1; end
    if idx >= length(x_ext), idx = length(x_ext) - 1; end

    hi = h_ext(idx);
    del = val - x_ext(idx);

    % 样条公式
    yy = y_ext(idx) ...
        + del * ( (y_ext(idx+1)-y_ext(idx))/hi - (2*M_ext(idx)+M_ext(idx+1))*hi/6 ) ...
        + del^2 * M_ext(idx)/2 ...
        + del^3 * (M_ext(idx+1)-M_ext(idx)) / (6*hi);
    y_query(k) = yy;
end

end