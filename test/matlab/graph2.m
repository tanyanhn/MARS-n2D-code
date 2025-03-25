close all
% 闭合曲线参数设置
R = 1;              % 基础半径

% 用户指定凸起最高点坐标 (示例坐标)
P1 = [1.125, 0.65];  % 第一个凸起最高点坐标（绝对坐标系）
P2 = [1.125, -0.65]; % 第二个凸起最高点坐标（绝对坐标系）
center = [-0.5 0];   % 曲线中心坐标

% 计算相对坐标
P1_rel = P1 - center; % 转换为相对于中心的坐标
P2_rel = P2 - center;

% 计算极坐标参数
theta1 = mod(atan2(P1_rel(2), P1_rel(1)), 2*pi);  % 角度规范到[0, 2pi)
H1 = norm(P1_rel);
theta2 = mod(atan2(P2_rel(2), P2_rel(1)), 2*pi);  % 角度规范到[0, 2pi)
H2 = norm(P2_rel);

% 陡峭度参数
k_left1 = 1000;     k_right1 = 1000; % 增大陡峭度参数
k_left2 = 1000;      k_right2 = 1000;
theta_sharp = pi;    % 尖锐点位置
C = -0.5;            k_sharp = 50;   % 增大尖锐度参数

% 计算幅度参数
A1 = H1 - R;
A2 = H2 - R;

% 生成角度数组
theta = linspace(0, 2*pi, 10000); % 增加采样点数

% 凸起计算函数（向量化计算）
calcBump = @(t, t0, A, kl, kr) A*exp(-(t < t0).*kl.*(t0 - t).^2 - (t >= t0).*kr.*(t - t0).^2);
term1 = calcBump(theta, theta1, A1, k_left1, k_right1);
term2 = calcBump(theta, theta2, A2, k_left2, k_right2);

% 尖锐点修正项
theta_diff = mod(theta - theta_sharp + pi, 2*pi) - pi;
sharp_term = C * abs(theta_diff) .* exp(-k_sharp * theta_diff.^2);

% 合成曲线
r = R + term1 + term2 + sharp_term;
x = r .* cos(theta) + center(1);
y = r .* sin(theta) + center(2);

% 计算x的数值导数
dxdtheta = gradient(x, theta);

% 寻找导数由正变负的交叉点
signs = sign(dxdtheta);
crossing_down = find(diff(signs) < 0);

% 标量计算函数
calcBump_single = @(t, t0, A, kl, kr) A * exp(-((t < t0) * kl * (t0 - t)^2 + (t >= t0) * kr * (t - t0)^2));
sharp_term_single = @(t) C * abs(mod(t - theta_sharp + pi, 2*pi) - pi) * exp(-k_sharp * (mod(t - theta_sharp + pi, 2*pi) - pi)^2);

% 寻找x极值点
maxima_theta = [];
for i = 1:length(crossing_down)
    idx = crossing_down(i);
    if idx >= length(theta)-1
        continue;
    end
    t_low = theta(idx);
    t_high = theta(idx + 1);
    
    % 定义x关于theta的函数（包含中心偏移）
    compute_x = @(t) (R + calcBump_single(t, theta1, A1, k_left1, k_right1) + ...
                      calcBump_single(t, theta2, A2, k_left2, k_right2) + ...
                      sharp_term_single(t)) * cos(t) + center(1);
    [t_opt, ~] = fminbnd(@(t) -compute_x(t), t_low, t_high, optimset('TolX', 1e-10));
    maxima_theta(end+1) = t_opt;
end

if isempty(maxima_theta)
    error('未找到极大值点');
end

% 绘制曲线
figure;
plot(x, y, 'LineWidth', 1.5);
axis equal; grid on;
title('闭合曲线及特征点');
xlabel('X轴'); ylabel('Y轴');
hold on;

% 绘制X轴极值点
for i = 1:length(maxima_theta)
    theta_max = maxima_theta(i);
    r_max = R + calcBump_single(theta_max, theta1, A1, k_left1, k_right1) + ...
             calcBump_single(theta_max, theta2, A2, k_left2, k_right2) + ...
             sharp_term_single(theta_max);
    x_max = r_max * cos(theta_max) + center(1);
    y_max = r_max * sin(theta_max) + center(2);
    plot(x_max, y_max, 'g*', 'MarkerSize', 10, 'LineWidth', 2);
end

% 精确计算特征点坐标
feature_points = [theta1, theta2, theta_sharp];
colors = ['r', 'm', 'k'];
markers = {'o', 'o', 'x'};

for i = 1:3
    theta_feature = feature_points(i);
    r_feature = R + calcBump_single(theta_feature, theta1, A1, k_left1, k_right1) + ...
                calcBump_single(theta_feature, theta2, A2, k_left2, k_right2) + ...
                sharp_term_single(theta_feature);
    x_feature = r_feature * cos(theta_feature) + center(1);
    y_feature = r_feature * sin(theta_feature) + center(2);
    plot(x_feature, y_feature, markers{i}, 'Color', colors(i), ...
         'MarkerSize', 12, 'LineWidth', 2);
end

legend('曲线', 'X轴极值点', '凸起1', '凸起2', '尖锐点', 'Location', 'best');