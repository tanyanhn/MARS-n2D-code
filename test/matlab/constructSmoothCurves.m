function curves = constructSmoothCurves(V, E, S, opts)
% CONSTRUCTSMOOTHCURVES 按光滑约束将边分组并采样
% 输入:
%   V   2xM 顶点矩阵（仅用于潜在扩展，当前未直接使用）
%   E   1xN Cell，每单元 struct(ctrl 8xk, k, v1, v2, idx)
%   S   1xM Cell，空或 struct(pairs, withFlag)，withFlag=[eA,fA,eB,fB]
%   opts.maxStep   相邻采样点距离上界（默认 Inf）
%   opts.minSamples 每段 Bezier 至少采样点数，含端点（默认 2）
% 输出:
%   curves Cell，每单元 struct(edges, ctrl, isClosed, samples, sampleBreaks)

    if nargin < 4
        opts = struct;
    end
    if ~isfield(opts, 'minSamples'); opts.minSamples = 2; end

    numEdges = numel(E);
    minSegLen = compute_min_seg_length(E);
    if ~isfield(opts, 'maxStep') || isempty(opts.maxStep)
        opts.maxStep = minSegLen;
    end

    % adjacency: nextEdge(edge+1, flag+1) -> edgeId (0-based), flag
    nextEdge = nan(numEdges, 2);
    nextFlag = nan(numEdges, 2);

    % 构建端点连接映射
    for v = 1:numel(S)
        entry = S{v};
        if isempty(entry)
            continue;
        end
        if isstruct(entry)
            wf = entry.withFlag;
        else
            wf = entry;
        end
        if isempty(wf)
            continue;
        end
        if size(wf, 2) ~= 4
            warning('光滑数据格式异常: 顶点 %d', v - 1);
            continue;
        end
        for r = 1:size(wf, 1)
            eA = wf(r, 1); fA = wf(r, 2); eB = wf(r, 3); fB = wf(r, 4);
            validate_flags(eA, fA, numEdges);
            validate_flags(eB, fB, numEdges);
            [nextEdge, nextFlag] = set_conn(nextEdge, nextFlag, eA, fA, eB, fB);
            [nextEdge, nextFlag] = set_conn(nextEdge, nextFlag, eB, fB, eA, fA);
        end
    end

    visited = false(1, numEdges);
    curves = {};

    for startId = 0:(numEdges - 1)
        if visited(startId + 1)
            continue;
        end
        visited(startId + 1) = true;

        % 初始边按正向放入
        edgeList = [];
        ctrlList = {};
        [edgeList, ctrlList] = append_edge(edgeList, ctrlList, E, startId, +1);

        % 向前（尾部）延伸
        currEdge = startId; currDir = +1; forwardOpen = false;
        while true
            outFlag = (currDir == 1); % 0/1
            nEdge = nextEdge(currEdge + 1, outFlag + 1);
            nFlag = nextFlag(currEdge + 1, outFlag + 1);
            if isnan(nEdge)
                forwardOpen = true;
                break;
            end
            if visited(nEdge + 1)
                break;
            end
            visited(nEdge + 1) = true;
            nDir = flag_to_dir(nFlag);
            [edgeList, ctrlList] = append_edge(edgeList, ctrlList, E, nEdge, nDir);
            currEdge = nEdge; currDir = nDir;
        end

        % 向后（头部）延伸
        currEdge = startId; currDir = +1; backwardOpen = false;
        while true
            outFlag = (currDir == 1); % 当前方向的起点 flag
            outFlag = 1 - outFlag;   % 取另一端
            nEdge = nextEdge(currEdge + 1, outFlag + 1);
            nFlag = nextFlag(currEdge + 1, outFlag + 1);
            if isnan(nEdge)
                backwardOpen = true;
                break;
            end
            if visited(nEdge + 1)
                break;
            end
            visited(nEdge + 1) = true;
            nDir = flag_to_dir(nFlag);
            [edgeList, ctrlList] = prepend_edge(edgeList, ctrlList, E, nEdge, nDir);
            currEdge = nEdge; currDir = nDir;
        end

        % 拼接 ctrl
        fullCtrl = cat(2, ctrlList{:});

        % 采样
        [samples, breaks] = sample_chain(ctrlList, opts.maxStep, opts.minSamples);

        % 判定闭合
        headDir = edgeList(1).dir; headId = edgeList(1).edgeId;
        tailDir = edgeList(end).dir; tailId = edgeList(end).edgeId;
        headFlag = (headDir == 1) * 0 + (headDir == -1) * 1;
        tailFlag = (tailDir == 1) * 1 + (tailDir == -1) * 0;
        isClosed = ~forwardOpen && ~backwardOpen;
        if isClosed
            % 严格检查首尾连通
            nf1 = nextEdge(headId + 1, headFlag + 1);
            nf2 = nextEdge(tailId + 1, tailFlag + 1);
            isClosed = ~isnan(nf1) && ~isnan(nf2);
        end

        c = struct( ...
            'edges', edgeList, ...
            'ctrl', fullCtrl, ...
            'isClosed', isClosed, ...
            'samples', samples, ...
            'sampleBreaks', breaks ...
        );
        curves{end + 1} = c; %#ok<AGROW>
    end
end

%% 辅助
function [nextEdge, nextFlag] = set_conn(nextEdge, nextFlag, eA, fA, eB, fB)
    if ~isnan(nextEdge(eA + 1, fA + 1))
        warning('端点重复光滑连接: edge=%d flag=%d', eA, fA);
    end
    nextEdge(eA + 1, fA + 1) = eB;
    nextFlag(eA + 1, fA + 1) = fB;
end

function validate_flags(e, f, numEdges)
    if e < 0 || e >= numEdges || (f ~= 0 && f ~= 1)
        error('光滑引用越界: edge=%d flag=%d', e, f);
    end
end

function dir = flag_to_dir(enterFlag)
% 进入 enterFlag 端，遍历方向
    if enterFlag == 0
        dir = +1; % v1->v2
    else
        dir = -1; % v2->v1
    end
end

function [edgeList, ctrlList] = append_edge(edgeList, ctrlList, E, edgeId, dir)
    ctrl = E{edgeId + 1}.ctrl;
    if dir == -1
        ctrl = reverse_ctrl(ctrl);
    end
    edgeList = [edgeList; struct('edgeId', edgeId, 'dir', dir, 'segCount', size(ctrl, 2))];
    ctrlList{end + 1} = ctrl; %#ok<AGROW>
end

function [edgeList, ctrlList] = prepend_edge(edgeList, ctrlList, E, edgeId, dir)
    ctrl = E{edgeId + 1}.ctrl;
    if dir == -1
        ctrl = reverse_ctrl(ctrl);
    end
    edgeList = [struct('edgeId', edgeId, 'dir', dir, 'segCount', size(ctrl, 2)); edgeList];
    ctrlList = [{ctrl}, ctrlList];
end

function revMat = reverse_ctrl(mat)
% mat 8 x K，反向列并交换 P0<->P3, P1<->P2
    mat = fliplr(mat);
    x_rows = mat(1:4, :);
    y_rows = mat(5:8, :);
    x_rows = x_rows([4, 3, 2, 1], :);
    y_rows = y_rows([4, 3, 2, 1], :);
    revMat = [x_rows; y_rows];
end

function [samples, breaks] = sample_chain(ctrlList, maxStep, minSamples)
    samples = [];
    breaks = [];
    segCounter = 0;
    for ci = 1:numel(ctrlList)
        ctrl = ctrlList{ci};
        k = size(ctrl, 2);
        for s = 1:k
            segCounter = segCounter + 1;
            pts = ctrl(:, s);
            pxy = [pts(1:4), pts(5:8)]; % 4x2 (rows are control points)
            segPts = sample_bezier(pxy, maxStep, minSamples);
            startIdx = size(samples, 1) + 1;
            if ~isempty(samples)
                segPts = segPts(2:end, :); % 去掉重复端点
            end
            samples = [samples; segPts]; %#ok<AGROW>
            breaks(segCounter, 1) = startIdx;
        end
    end
    % 闭合时去掉重复首尾
    if size(samples, 1) > 1 && norm(samples(end, :) - samples(1, :)) < 1e-9
        samples = samples(1:end-1, :);
    end
end

function pts = sample_bezier(p, maxStep, minSamples)
% p: 4x2 控制点
    chordLen = sum(vecnorm(diff(p, 1, 1), 2, 2));
    intervals = 1;
    if isfinite(maxStep)
        intervals = max(intervals, ceil(chordLen / maxStep));
    end
    intervals = max(intervals, max(minSamples - 1, 1));
    t = linspace(0, 1, intervals + 1);
    pts = zeros(numel(t), 2);
    % Bernstein形式
    for i = 1:numel(t)
        tt = t(i);
        b0 = (1 - tt)^3;
        b1 = 3 * (1 - tt)^2 * tt;
        b2 = 3 * (1 - tt) * tt^2;
        b3 = tt^3;
        pts(i, :) = b0 * p(1, :) + b1 * p(2, :) + b2 * p(3, :) + b3 * p(4, :);
    end
end

function minLen = compute_min_seg_length(E)
    minLen = Inf;
    for i = 1:numel(E)
        ctrl = E{i}.ctrl;
        k = size(ctrl, 2);
        for s = 1:k
            pts = ctrl(:, s);
            pxy = [pts(1:4), pts(5:8)]; % 4x2 (rows are control points)
            chordLen = sum(vecnorm(diff(pxy, 1, 1), 2, 2));
            if chordLen < minLen
                minLen = chordLen;
            end
        end
    end
    if ~isfinite(minLen) || minLen <= 0
        minLen = Inf;
    end
end
