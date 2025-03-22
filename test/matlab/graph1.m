% 参数设置
theta = linspace(0, 2*pi, 1000); % 角度采样
r_base = 0.5;       % 基础半径

% 第一个凸起参数（左侧平缓，右侧陡峭）
h1 = 0.2;         % 凸起高度
theta1 = pi/2 - 0.5;  % 凸起中心位置
k1_left = 2;     % 左侧衰减系数（值越小衰减越慢）
k1_right = 8;    % 右侧衰减系数（值越大衰减越快）

% 第二个凸起参数（左侧陡峭，右侧平缓）
h2 = 0.2;
theta2 = pi/2 + 0.5;
k2_left = 7;
k2_right = 3;

% 计算周期角度差
delta1 = mod(theta - theta1 + pi, 2*pi) - pi;
delta2 = mod(theta - theta2 + pi, 2*pi) - pi;

% 非对称高斯函数
bump1 = h1 * exp( - ( (delta1 < 0)*k1_left + (delta1 >= 0)*k1_right ) .* delta1.^2 );
bump2 = h2 * exp( - ( (delta2 < 0)*k2_left + (delta2 >= 0)*k2_right ) .* delta2.^2 );

% 合成半径
r = r_base + bump1 + bump2;

% 转换为笛卡尔坐标
x = r .* cos(theta);
y = r .* sin(theta);

% 可视化
figure;
plot(x, y, 'LineWidth', 2);
axis equal;
title('Asymmetric Double Bump Curve');
xlabel('X'); ylabel('Y');
grid on;
hold on;
plot(0, 0, 'r+'); % 标记原点

% 绘制凸起中心参考线
plot(r_base*cos(theta1), r_base*sin(theta1), 'bo');
plot(r_base*cos(theta2), r_base*sin(theta2), 'bo');
