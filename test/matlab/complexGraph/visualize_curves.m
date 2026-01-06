function visualize_curves(curves, E, V, F, fillIdx, samplesPerBreak)
% VISUALIZE_CURVES 展示光滑曲线集合
% curves: constructSmoothCurves 的输出
% E: 边列表（包含 v1/v2）
% V: 顶点矩阵
% F: 颜色分组信息（来自 parse_geometry_json，可选）
% fillIdx: 可选，F 下标数组。仅渲染这些下标对应的 boundary，按 boundary 的定向决定实心/空洞。
% samplesPerBreak: 每两个相邻 break 之间的采样点数（含端点，可选）

    if nargin < 4, F = []; end
    if nargin < 5, fillIdx = []; end
    if nargin < 6 || isempty(samplesPerBreak)
        samplesPerBreak = 10;
    end

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

        % 绘制曲线（优先使用 spline 直接取点）
        if isfield(c, 'spline') && isstruct(c.spline) && ...
                isfield(c.spline, 'xpp') && isfield(c.spline, 'ypp')
            sp = c.spline;
            if isfield(sp, 'breaks') && ~isempty(sp.breaks)
                breaks = sp.breaks;
            else
                breaks = sp.xpp.breaks;
            end
            pts = sample_spline(sp.xpp, sp.ypp, breaks, samplesPerBreak);
            plot(pts(:, 1), pts(:, 2), 'Color', col, 'LineWidth', 1.5);
        elseif isfield(c, 'xpp') && isfield(c, 'ypp')
            breaks = c.xpp.breaks;
            pts = sample_spline(c.xpp, c.ypp, breaks, samplesPerBreak);
            plot(pts(:, 1), pts(:, 2), 'Color', col, 'LineWidth', 1.5);
        else
            plot(c.samples(:, 1), c.samples(:, 2), 'Color', col, ...
                'LineWidth', 1.5);
        end

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
            scatter(coords(:, 1), coords(:, 2), 60, 'b', 'o');
                % 'DisplayName', sprintf('Internal #%d', i));
        end

        % 端点（方块）
        if ~isempty(endpoints)
            coords = V(:, endpoints + 1)';
            scatter(coords(:, 1), coords(:, 2), 30, 'r', 's');
                % 'DisplayName', sprintf('Endpoints #%d', i));
        end
    end

    legend('show');
    
    % 按 F 和 fillIdx 进行区域填充（使用 boundary 定向区分空洞）
    if nargin >= 4 && ~isempty(F) && nargin >= 5 && ~isempty(fillIdx)
        alphaVal = 0.15;
        sel = unique(fillIdx);
        sel = sel(sel >= 1 & sel <= numel(F)); % 合法索引
        for fi = sel(:).'
            faceColor = double(F(fi).color) / 255;
            xs = {}; ys = {};
            for bi = 1:numel(F(fi).boundaries)
                loop = F(fi).boundaries{bi}.vertexLoop;
                coords = V(:, loop + 1);
                xs{end+1} = coords(1, :); %#ok<AGROW>
                ys{end+1} = coords(2, :); %#ok<AGROW>
            end
            try
                pg = polyshape(xs, ys, 'Simplify', false); % 定向区分空洞/实体
                plot(pg, 'FaceColor', faceColor, 'FaceAlpha', alphaVal, 'EdgeColor', 'none');
            catch
                % 回退到逐 boundary 填充（不处理空洞）
                for bi = 1:numel(xs)
                    fill(xs{bi}, ys{bi}, faceColor, ...
                        'FaceAlpha', alphaVal, 'EdgeColor', 'none');
                end
            end
        end
    end

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

function pts = sample_spline(xpp, ypp, breaks, samplesPerBreak)
    pts = zeros(0, 2);
    if numel(breaks) < 2
        return;
    end
    samplesPerBreak = max(2, round(samplesPerBreak));
    for bi = 1:(numel(breaks) - 1)
        t0 = breaks(bi);
        t1 = breaks(bi + 1);
        t = linspace(t0, t1, samplesPerBreak);
        segPts = [ppval(xpp, t)', ppval(ypp, t)'];
        if ~isempty(pts)
            segPts = segPts(2:end, :);
        end
        pts = [pts; segPts]; %#ok<AGROW>
    end
end
