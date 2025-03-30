function quintic_periodic_spline(points)
% 输入点列，构造五次周期样条插值
% points: N×2矩阵，每行代表一个点的x,y坐标
% 自动闭合点列

if ~isequal(points(1,:), points(end,:))
    points(end+1,:) = points(1,:);
    fprintf('提示：自动闭合点列，添加第%d个点\n', size(points,1));
end

n = size(points, 1) - 1; % 实际控制点数量
chords = cumsum([0; sqrt(sum(diff(points,1,1).^2,2))]); % 弦长参数化
t = chords/chords(end);
t(end) = []; % 移除重复点的参数

x = points(1:n, 1);
y = points(1:n, 2);

order = 6; % 五次样条阶数
ppx = periodic_spline(t, x, order);
ppy = periodic_spline(t, y, order);

s = linspace(0, 1, 1000);
xx = ppval(ppx, s);
yy = ppval(ppy, s);

figure;
plot(xx, yy, 'b-', 'LineWidth', 2);
hold on;
plot(points(:,1), points(:,2), 'ro', 'MarkerFaceColor', 'r');
title('五次周期样条插值');
legend('插值曲线', '原始点列');
axis equal;
grid on;
end

function pp = periodic_spline(t, vals, order)
% 正确构造周期样条节点
n = length(t);
knots = linspace(0, 1, n); % 基础节点

% 手动延拓节点实现周期性
ext_front = knots(end) - (1:(order-1))*(knots(2)-knots(1));
ext_end = knots(1) + (1:(order-1))*(knots(2)-knots(1));
full_knots = sort([ext_front, knots, ext_end]);
full_knots = mod(full_knots, 1); % 归一化到[0,1)

% 构造Collocation矩阵
A = spcol(full_knots, order, t(:)');

% 添加周期性导数约束
for d = 0:4
    A_deriv0 = spcol(full_knots, order, 0, 'sparse', d);
    A_deriv1 = spcol(full_knots, order, 1, 'sparse', d);
    A = [A; A_deriv0 - A_deriv1];
end
b = [vals(:); zeros(5,1)];

% 解约束方程组
coeffs = A \ b;

% 创建样条对象
pp = spmak(full_knots, coeffs.');
end
