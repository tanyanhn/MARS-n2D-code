function [xsp, ysp, t] = construct_quintic_spline(points)
% CONSTRUCT_QUINTIC_SPLINE 构建5次样条(C4连续)
% 输入：points是2×n矩阵
% 输出：xsp, ysp为pp形式的5次样条，Breaks与输入点完全对齐

% --- 1. 基础处理 ---
[m, n] = size(points);
if m ~= 2, error('需2xn矩阵'); end
if n < 2, error('至少需2个点'); end

% 计算累积距离 t
diffs = diff(points, 1, 2);
d = sqrt(sum(diffs.^2, 1));
t = [0, cumsum(sqrt(d))];

% --- 2. 准备 5 次样条的节点 (Knots) ---
k = 6; % 5次样条对应阶数 k=6
% 使用 augknt 在每个 t 上放置节点，端点重复 k 次
% 这保证了转换后的 pp-form 的 breaks 刚好就是 t
knots = augknt(t, k); 

% --- 3. 估算边界导数 (核心步骤) ---
% 5次样条需要 n+4 个条件。我们有 n 个点，缺 4 个。
% 我们缺的是：起点的一阶、二阶导数；终点的一阶、二阶导数。
% 策略：用 'variational' 三次样条来估算这些物理导数。
temp_x = csape(t, points(1,:), 'variational');
temp_y = csape(t, points(2,:), 'variational');

% 提取起点的 1阶(v) 和 2阶(a) 导数
[vx_start, ax_start] = get_derivatives(temp_x, t(1));
[vy_start, ay_start] = get_derivatives(temp_y, t(1));

% 提取终点的 1阶(v) 和 2阶(a) 导数
[vx_end, ax_end] = get_derivatives(temp_x, t(end));
[vy_end, ay_end] = get_derivatives(temp_y, t(end));

% --- 4. 构建 Hermite 插值数据 ---
% 构造 sites 向量：端点重复3次(表示匹配值、1阶导、2阶导)，中间点1次
% 结构：[t1, t1, t1, t2, t3, ..., tn-1, tn, tn, tn]
sites = [t(1), t(1), t(1), t(2:end-1), t(end), t(end), t(end)];

% 构造 X 维度的值向量：对应 sites 的顺序
% [x1, vx1, ax1, x2, x3, ..., xn-1, xn, vxn, axn]
vals_x = [points(1,1), vx_start, ax_start, ...      % 起点 (3个条件)
          points(1,2:end-1), ...                     % 中间点 (n-2个条件)
          points(1,end), vx_end, ax_end];            % 终点 (3个条件)
          
% 构造 Y 维度的值向量
vals_y = [points(2,1), vy_start, ay_start, ...
          points(2,2:end-1), ...
          points(2,end), vy_end, ay_end];

% --- 5. 生成 5 次样条 ---
% spapi 支持重复 sites，会自动识别为导数插值
sp_x = spapi(knots, sites, vals_x);
sp_y = spapi(knots, sites, vals_y);

% --- 6. 转换为 pp 形式 ---
xsp = fn2fm(sp_x, 'pp');
ysp = fn2fm(sp_y, 'pp');

% 验证 (调试用)
% fprintf('Points数: %d, Breaks数: %d, Order: %d\n', length(t), length(xsp.breaks), xsp.order);
end

function [d1, d2] = get_derivatives(pp, site)
    % 辅助函数：获取某点的 1阶 和 2阶 导数值
    d1 = fnval(fnder(pp, 1), site);
    d2 = fnval(fnder(pp, 2), site);
end