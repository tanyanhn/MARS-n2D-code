close all 
clc
num = 1e6;
% 参数设置
theta = linspace(0, 2*pi, num); % 角度采样
r_base = [0.2, 0.25];       % 基础半径
center = [0.3 0.5];

% 第一个凸起参数（左侧平缓，右侧陡峭）
h1 = 0.15;         % 凸起高度
theta1 = 0 - 0.6;  % 凸起中心位置
k1_left = 3;     % 左侧衰减系数（值越小衰减越慢）
k1_right = 20;    % 右侧衰减系数（值越大衰减越快）

% 第二个凸起参数（左侧陡峭，右侧平缓）
h2 = 0.25;
theta2 = 0 + 0.5;
k2_left = 10;
k2_right = 3;

% === 新增：尖点参数 ===
theta_spike = pi;    % 尖点角度
C = -0.4;             % 尖点强度（控制尖角锐度）
k_sharp = 10;      % 衰减系数（越大尖点越局部化）

% 计算周期角度差
delta1 = mod(theta - theta1 + pi, 2*pi) - pi;
delta2 = mod(theta - theta2 + pi, 2*pi) - pi;
delta_spike = mod(theta - theta_spike + pi, 2*pi) - pi; % 尖点角度差

% 非对称高斯函数
bump1 = h1 * exp( - ( (delta1 < 0)*k1_left + (delta1 >= 0)*k1_right ) .* delta1.^2 );
bump2 = h2 * exp( - ( (delta2 < 0)*k2_left + (delta2 >= 0)*k2_right ) .* delta2.^2 );
sharp_term = C * abs(delta_spike) .* exp(-k_sharp * delta_spike.^2);
sharp_term = 0;

% 合成坐标（添加尖点）
x = center(1) + r_base(1)*cos(theta) + (bump1 + bump2 + sharp_term) .* cos(theta);
y = center(2) + r_base(2)*sin(theta) + (bump1 + bump2 + sharp_term) .* sin(theta);

% 可视化
figure;
plot(x, y, 'LineWidth', 1);
axis equal;
title('Asymmetric Double Bump Curve');
xlabel('X'); ylabel('Y');
grid on;
hold on;
hold on;


% 曲线二
theta = linspace(0, 2*pi, num); % 角度采样
r_base = [0.03 0.2];       % 基础半径
center = [0.83 0.5];

% 第一个凸起参数（左侧平缓，右侧陡峭）
h1 = 0.24;         % 凸起高度
theta1 = 0 - 2.6;  % 凸起中心位置
k1_left = 20;     % 左侧衰减系数（值越小衰减越慢）
k1_right = 4;    % 右侧衰减系数（值越大衰减越快）

% 第二个凸起参数（左侧陡峭，右侧平缓）
h2 = 0.1352294955;
theta2 = 0 + 2.5;
k2_left = 10;
k2_right = 20;

% === 新增：尖点参数 ===
theta_spike = 1.5;    % 尖点角度
C = -0.05;             % 尖点强度（控制尖角锐度）
k_sharp_left = 5;       % 衰减系数（越大尖点越局部化）
k_sharp_right = 5;
% 计算周期角度差
delta1 = mod(theta - theta1 + pi, 2*pi) - pi;
delta2 = mod(theta - theta2 + pi, 2*pi) - pi;
delta_spike = mod(theta - theta_spike + pi, 2*pi) - pi; % 尖点角度差

% 非对称高斯函数
bump1 = h1 * exp( - ( (delta1 < 0)*k1_left + (delta1 >= 0)*k1_right ) .* delta1.^2 );
bump2 = h2 * exp( - ( (delta2 < 0)*k2_left + (delta2 >= 0)*k2_right ) .* delta2.^2 );
sharp_term = C * abs(delta_spike) .* exp(- ...
    ( (delta_spike < 0)*k_sharp_left + (delta_spike >= 0)*k_sharp_right ) .* (2 *delta_spike).^2);

% 合成坐标（添加尖点）
x = center(1) + r_base(1)*cos(theta) + (bump1 + bump2 + sharp_term) .* cos(theta);
y = center(2) + r_base(2)*sin(theta) + (bump1 + bump2 + sharp_term) .* sin(theta);

plot(x, y, 'LineWidth', 1);

axis([0 1 0 1]);
