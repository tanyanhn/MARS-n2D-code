function sampleE = split_curve_samples(curves, E)
% SPLIT_CURVE_SAMPLES 将曲线采样拆分回逐 Edge 的采样
% 输入:
%   curves: constructSmoothCurves 的输出 (含 fields: edges, samples, sampleBreaks)
%   E:      原始边列表 (1xN cell)
% 输出:
%   sampleE: 1xN cell，sampleE{i} 为该边的采样点 (Nx2)，方向为 E{i} 的 v1->v2。

    numEdges = numel(E);
    sampleE = cell(1, numEdges);

    for ci = 1:numel(curves)
        c = curves{ci};
        segStart = c.sampleBreaks; % 每段 Bezier 的起始索引
        segStart = segStart(:)';    % 行向量
        segStart(end + 1) = size(c.samples, 1); % 便于取终点

        segIdx = 1; % 全局段索引
        for ei = 1:numel(c.edges)
            eInfo = c.edges(ei);
            k = eInfo.segCount;
            if segIdx + k - 1 >= numel(segStart)
                warning('sampleBreaks 与 segCount 不一致，曲线 %d 边序 %d', ci, ei);
                break;
            end

            % 当前边覆盖的采样范围
            sStart = segStart(segIdx);
            sEnd   = segStart(segIdx + k);
            segIdx = segIdx + k;

            pts = c.samples(sStart:sEnd, :);

            % 若曲线方向与 E 定义相反，则翻转采样点以保持端点一致 (v1->v2)
            if eInfo.dir == -1
                pts = flipud(pts);
            end

            edgeId = eInfo.edgeId + 1; % MATLAB 索引
            % 若同一边被多条曲线引用，后出现的覆盖前者
            sampleE{edgeId} = pts;
        end
    end
end
