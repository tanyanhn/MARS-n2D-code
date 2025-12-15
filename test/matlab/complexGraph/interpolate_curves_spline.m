function [curvesOut, splineE] = interpolate_curves_spline(curves, E)
% INTERPOLATE_CURVES_SPLINE 对曲线采样做五次样条插值并按边切分
% 输入:
%   curves: constructSmoothCurves 的输出 (需包含 samples, sampleBreaks, edges, isClosed)
%   E:      原始边列表 (1xN cell)
% 输出:
%   curvesOut: 在 curves 基础上增加 spline 字段 (xpp, ypp, breaks, tSamples)
%   splineE:   1xN cell, 每个单元对应 E{i} 的样条子段元数据

    numEdges = numel(E);
    splineE = cell(1, numEdges);
    curvesOut = curves;

    for ci = 1:numel(curves)
        c = curves{ci};
        pts = c.samples;
        if size(pts, 1) < 2
            warning('曲线 %d 采样点不足，跳过样条插值', ci);
            continue;
        end

        % 弦长参数 (0..1)
        d = sqrt(sum(diff(pts, 1, 1).^2, 2));
        if sum(d) < eps
            tSamples = zeros(size(pts, 1), 1);
        else
            s = [0; cumsum(d)];
            tSamples = s / s(end);
        end

        % 五次样条插值（使用已有 quinticSpline）
        try
            if c.isClosed
                [xpp, ypp, breaks] = quinticSpline(pts');
            else
                [xpp, ypp, breaks] = construct_spline(pts', 6);
            end
        catch ME
            error('曲线 %d 样条插值失败: %s', ci, ME.message);
        end

        curvesOut{ci}.spline = struct( ...
            'xpp', xpp, ...
            'ypp', ypp, ...
            'breaks', breaks ...
        );

        % 按边切分到 splineE
        segStart = c.sampleBreaks(:)';   % 每段 Bezier 起始索引
        segStart(end + 1) = size(pts, 1) + 1;
        segIdx = 1;

        for ei = 1:numel(c.edges)
            eInfo = c.edges(ei);
            k = eInfo.segCount;
            if segIdx + k > numel(segStart)
                warning('曲线 %d 的 sampleBreaks 与 segCount 不一致', ci);
                break;
            end
            sStart = segStart(segIdx);
            sEnd   = segStart(segIdx + k) - 1;
            segIdx = segIdx + k;

            t0 = tSamples(sStart);
            t1 = tSamples(sEnd);
            dir = eInfo.dir;
            if dir == -1
                [t0, t1] = deal(t1, t0);
            end

            edgeId = eInfo.edgeId + 1; % MATLAB 索引
            if edgeId < 1 || edgeId > numEdges
                warning('曲线 %d 引用非法 edgeId=%d', ci, edgeId - 1);
                continue;
            end

            splineE{edgeId} = struct( ...
                'curveId', ci, ...
                'dir', dir, ...
                't0', t0, ...
                't1', t1, ...
                'xpp', xpp, ...
                'ypp', ypp ...
            );
        end
    end
end
