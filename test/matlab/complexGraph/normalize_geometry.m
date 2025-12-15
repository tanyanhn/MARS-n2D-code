function [Vn, En, S] = normalize_geometry(V, E, S, rect, theta)
% NORMALIZE_GEOMETRY 将几何坐标旋转（可选）后缩放/平移到指定矩形
% 输入:
%   V    2xM 顶点矩阵
%   E    1xN Cell，单元 struct(ctrl 8xk, k, v1, v2, idx)
%   S    光滑结构 (原样返回)
%   rect [lr, lc, rr, rc] 目标范围 (x_min, y_min, x_max, y_max)
%   theta (可选) 旋转角度（弧度），绕源包围盒中心逆时针旋转
% 输出:
%   Vn, En 为正则化后的几何，S 原样透传

    if nargin < 5
        theta = 0;
    end

    x0 = rect(1); y0 = rect(2); x1 = rect(3); y1 = rect(4);
    tgtW = x1 - x0; tgtH = y1 - y0;
    if tgtW <= 0 || tgtH <= 0
        error('目标范围无效: 宽高需为正 (rect = [%g %g %g %g])', rect);
    end

    % 收集所有点
    allPts = V;
    for i = 1:numel(E)
        ctrl = E{i}.ctrl;
        a = ctrl(1:4, :);
        b = ctrl(5:8, :);
        pts = [a(:)'; b(:)']; % 2 x (4*k)
        allPts = [allPts, pts]; %#ok<AGROW>
    end

    % 原始包围盒与旋转中心
    srcMinX = min(allPts(1, :)); srcMaxX = max(allPts(1, :));
    srcMinY = min(allPts(2, :)); srcMaxY = max(allPts(2, :));
    center = [(srcMinX + srcMaxX) / 2; (srcMinY + srcMaxY) / 2];

    % 旋转矩阵
    cth = cos(theta); sth = sin(theta);
    rot = [cth, -sth; sth, cth];

    % 旋转后的全局包围盒
    rotAll = rot * (allPts - center) + center;
    rotMinX = min(rotAll(1, :)); rotMaxX = max(rotAll(1, :));
    rotMinY = min(rotAll(2, :)); rotMaxY = max(rotAll(2, :));
    rotW = rotMaxX - rotMinX; rotH = rotMaxY - rotMinY;

    if rotW < eps
        scaleX = 0; constX = x0 + tgtW / 2;
    else
        scaleX = tgtW / rotW; constX = x0;
    end
    if rotH < eps
        scaleY = 0; constY = y0 + tgtH / 2;
    else
        scaleY = tgtH / rotH; constY = y0;
    end

    % 统一的变换：先绕 center 旋转，再按旋转后包围盒缩放/平移
    function ptsT = apply_transform(pts)
        ptsR = rot * (pts - center) + center;
        x = ptsR(1, :); y = ptsR(2, :);
        if scaleX == 0
            xT = constX * ones(size(x));
        else
            xT = constX + (x - rotMinX) * scaleX;
        end
        if scaleY == 0
            yT = constY * ones(size(y));
        else
            yT = constY + (y - rotMinY) * scaleY;
        end
        ptsT = [xT; yT];
    end

    % V
    Vn = apply_transform(V);

    % E
    En = E;
    for i = 1:numel(E)
        ctrl = E{i}.ctrl;
        a = ctrl(1:4, :);
        b = ctrl(5:8, :);
        pts = [a(:)'; b(:)']; % 2 x (4*k)
        ptsT = apply_transform(pts);
        ctrl(1:4, :) = reshape(ptsT(1, :), 4, []);
        ctrl(5:8, :) = reshape(ptsT(2, :), 4, []);
        En{i} = E{i};
        En{i}.ctrl = ctrl;
    end
end
