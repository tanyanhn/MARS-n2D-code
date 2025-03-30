% testConstruct_spline
clear  
close all

t = linspace(0, 2 * pi, 400);
r = 1;

x = r * cos(t);
y = r * sin(t);
pts = [x; y];
[xsp, ysp] = quinticSpline(pts);
t_total = xsp.breaks(end);
t_eval = linspace(0, t_total, 100);
x_eval = ppval(xsp, t_eval);
y_eval = ppval(ysp, t_eval);
% 可视化
figure
plot(pts(1,:), pts(2,:), 'ro', x_eval, y_eval, 'b-');
legend('原始点', '样条轨迹');