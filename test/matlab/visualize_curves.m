function visualize_curves(curves, E, V)
% VISUALIZE_CURVES 展示光滑曲线集合
% curves: constructSmoothCurves 的输出
% E: 边列表（包含 v1/v2）
% V: 顶点矩阵

    if isempty(curves)
        warning('没有曲线数据可视化');
        return;
    end
    hold on;
    axis equal;
    grid on;
    title('Curves with Smooth Connections');
    colors = lines(numel(curves));

    for i = 1:numel(curves)
        c = curves{i};
        col = colors(i, :);

        % 绘制曲线采样
        plot(c.samples(:, 1), c.samples(:, 2), 'Color', col, 'LineWidth', 1.5);

        % 顶点序列
        vSeq = curve_vertex_sequence(c.edges, E);

        if c.isClosed
            endpoints = [];
            internal = unique(vSeq(1:end-1)); % 去掉重复闭合点
        else
            if numel(vSeq) >= 2
                endpoints = [vSeq(1), vSeq(end)];
                internal = vSeq(2:end-1);
            else
                endpoints = vSeq;
                internal = [];
            end
        end

        % 内部光滑连接点（圆点）
        if ~isempty(internal)
            coords = V(:, internal + 1)';
            scatter(coords(:, 1), coords(:, 2), 40, col, 'filled', 'o', ...
                'DisplayName', sprintf('Internal #%d', i));
        end

        % 端点（方块）
        if ~isempty(endpoints)
            coords = V(:, endpoints + 1)';
            scatter(coords(:, 1), coords(:, 2), 60, 'r', 'filled', 's', ...
                'DisplayName', sprintf('Endpoints #%d', i));
        end
    end

    legend('show');
    hold off;
end

function vSeq = curve_vertex_sequence(edgeList, E)
% 生成曲线的顶点序列（0-based）
    vSeq = [];
    for idx = 1:numel(edgeList)
        eInfo = edgeList(idx);
        e = E{eInfo.edgeId + 1};
        if eInfo.dir == 1
            startV = e.v1; endV = e.v2;
        else
            startV = e.v2; endV = e.v1;
        end
        if isempty(vSeq)
            vSeq = [startV, endV];
        else
            vSeq = [vSeq, endV]; %#ok<AGROW>
        end
    end
end
