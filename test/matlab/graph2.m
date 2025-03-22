% 闭合曲线参数设置
R = 1;              % 基础半径
theta1 = pi/3;      % 第一个凸起位置（60度）
theta2 = 2*pi/3;    % 第二个凸起位置（120度）
A1 = 0.3;           % 第一个凸起幅度
A2 = 0.3;           % 第二个凸起幅度
k_left1 = 100;      % 第一个凸起左侧陡峭度
k_right1 = 50;      % 第一个凸起右侧陡峭度
k_left2 = 50;       % 第二个凸起左侧陡峭度
k_right2 = 100;     % 第二个凸起右侧陡峭度
theta_sharp = pi;   % 尖锐点位置（180度）
C = -0.5;           % 尖锐点强度（调整后需增大）
k_sharp = 5;       % 尖锐点宽度控制（值越大影响范围越小）

% 生成角度数组
theta = linspace(0, 2*pi, 1000);

% 计算凸起项（保持不变）
calcBump = @(t, t0, A, kl, kr) A*exp(-(t < t0).*kl.*(t0 - t).^2 - (t >= t0).*kr.*(t - t0).^2);
term1 = calcBump(theta, theta1, A1, k_left1, k_right1);
term2 = calcBump(theta, theta2, A2, k_left2, k_right2);

% 计算改进的尖锐点项
theta_diff = mod(theta - theta_sharp + pi, 2*pi) - pi;  % 规范到[-π, π)
sharp_term = C * abs(theta_diff) .* exp(-k_sharp * theta_diff.^2); % 指数衰减修正

% 计算总半径
r = R + term1 + term2 + sharp_term;

% 转换为笛卡尔坐标并绘制
x = r .* cos(theta);
y = r .* sin(theta);

figure;
plot(x, y, 'LineWidth', 1.5);
axis equal; grid on;
title('改进的单尖锐点闭合曲线');
xlabel('X轴'); ylabel('Y轴');

% 标注特征点
hold on;
plot(r(theta == theta1)*cos(theta1), r(theta == theta1)*sin(theta1), 'ro', 'MarkerSize', 8);
plot(r(theta == theta2)*cos(theta2), r(theta == theta2)*sin(theta2), 'mo', 'MarkerSize', 8);
plot(r(theta == theta_sharp)*cos(theta_sharp), r(theta == theta_sharp)*sin(theta_sharp),...
     'kx', 'MarkerSize', 12, 'LineWidth', 2);
legend('曲线', '凸起1', '凸起2', '尖锐点');

% 验证对称位置（可添加检查代码）
% sym_theta = mod(theta_sharp + pi, 2*pi);
% sym_index = find(abs(theta - sym_theta) == min(abs(theta - sym_theta)), 1);
% plot(x(sym_index), y(sym_index), 'g*', 'MarkerSize', 10); % 对称位置标记