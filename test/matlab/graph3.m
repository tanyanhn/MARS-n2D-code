% 椭圆参数设置
a = 2;          % 长轴
b = 1;          % 短轴
k_values = [0, 0.5, 1, 1.5];  % 变形参数k的取值
theta = linspace(0, 2*pi, 1000);  % 角度采样

% 绘制不同k值的椭圆曲线
figure;
hold on;
colors = lines(length(k_values));  % 生成不同颜色

for i = 1:length(k_values)
    k = k_values(i);
    % 参数方程生成坐标
    y = b * sin(theta);
    x = a * cos(theta) + k * (b * sin(theta)).^2;
    % 绘制曲线
    plot(x, y, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', ['k=', num2str(k)]);
end

hold off;

% 图形标注
legend('show');
title('椭圆曲线向右变形效果 (a=2, b=1)');
xlabel('x'); ylabel('y');
grid on;
axis equal;
xlim([-3, 5]);  % 手动调整x轴范围以完整显示
