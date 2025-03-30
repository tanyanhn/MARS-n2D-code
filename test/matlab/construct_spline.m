function [xsp, ysp, t] = construct_spline(points, k)
% 输入：points是2×n矩阵，每列为插值点的x,y坐标
% 输出：xsp, ysp为pp样条，以累积距离t为参数

if nargin < 2
    k = 4;
end

% 检查输入有效性
[m, n] = size(points);
if m ~= 2
    error('输入必须是2行n列的坐标矩阵');
end
if n < 2
    error('至少需要2个点进行插值');
end

% 计算相邻点间欧氏距离
diffs = diff(points, 1, 2);          % 列方向差分
d = sqrt(sum(diffs.^2, 1));          % 各段距离
if any(d < eps*1e4)                  % 排除重合点
    error('存在重合相邻点，请检查输入数据');
end

% 生成累积距离参数t (从0开始)
t = [0, cumsum(d)];
knots = aptknt(t,k);
spline_x = spapi(knots, t, points(1,:));
spline_y = spapi(knots, t, points(2,:));

% 转换为pp格式
xsp = fn2fm(spline_x, 'pp');
ysp = fn2fm(spline_y, 'pp');

end
