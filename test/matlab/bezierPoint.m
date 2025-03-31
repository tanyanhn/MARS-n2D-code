
% 三次贝塞尔曲线单点计算函数
function pt = bezierPoint(control_points, t)
    % 将t转换为列向量以确保维度正确
    t = t(:);
    
    % 计算各贝塞尔基函数的值（均为列向量）
    c0 = (1 - t).^3;
    c1 = 3*(1 - t).^2.*t;
    c2 = 3*(1 - t).*t.^2;
    c3 = t.^3;
    
    % 确保控制点是正确的维度（4×dim）
    % 通过列方向扩展进行矩阵乘法，最终得到n×dim矩阵
    pt = c0 * control_points(1,:) + ...
         c1 * control_points(2,:) + ...
         c2 * control_points(3,:) + ...
         c3 * control_points(4,:);
end
