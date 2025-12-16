function newBoundaries = merge_boundaries(boundaries, parent, id, sampleE, tol)
% MERGE_BOUNDARIES 将子树所有边合并、抵消反向重复，并重建闭合边界
% 输入:
%   boundaries: build_boundaries_from_F 生成的 Cell
%   parent: build_containment_tree 的输出
%   id: 作为外边界的根下标，其子树视为洞（方向取反）
%   sampleE: split_curve_samples 的输出
%   tol: 容差（默认 1e-7）
% 输出:
%   newBoundaries: Cell，包含所有重建后的闭合环（edgeId 唯一，无反向重叠）

    if nargin < 5 || isempty(tol)
        tol = 1e-7;
    end

    idxSet = collect_subtree(parent, id);
    if length(idxSet) < 2
        newBoundaries = boundaries(id);
        return;
    end
    idxSet = idxSet(:).';

    % 聚合边，洞取反方向，并做反向抵消（累加方向，净为0则删除）
    dirSum = containers.Map('KeyType','double','ValueType','double');
    for b = 1:numel(idxSet)
        src = boundaries{idxSet(b)};
        edges = src.edges;
        if b > 1
            for k = 1:numel(edges)
                edges(k).dir = -edges(k).dir;
            end
        end
        for k = 1:numel(edges)
            eId = edges(k).edgeId;
            d = edges(k).dir;
            if isKey(dirSum, eId)
                dirSum(eId) = dirSum(eId) + d;
            else
                dirSum(eId) = d;
            end
        end
    end

    % 构造剩余边列表，方向取符号
    remEdges = {};
    keysE = dirSum.keys;
    for kk = 1:numel(keysE)
        eId = keysE{kk};
        dsum = dirSum(eId);
        if abs(dsum) < eps
            continue; % 完全抵消
        end
        dir = sign(dsum);
        remEdges{end+1} = struct('edgeId', eId, 'dir', dir); %#ok<AGROW>
    end

    % 用剩余边重建闭合环
    loops = build_loops(remEdges, sampleE, tol);

    % 颜色继承根节点
    col = boundaries{idxSet(1)}.color;
    newBoundaries = {};
    for i = 1:numel(loops)
        newBoundaries{end+1} = struct( ... %#ok<AGROW>
            'color', col, ...
            'edges', loops{i}, ...
            'vertexLoop', [] ...
        );
    end
end

function loops = build_loops(edges, sampleE, tol)
    loops = {};
    if isempty(edges), return; end
    remaining = edges; % cell of struct
    while ~isempty(remaining)
        curr = remaining{1};
        remaining(1) = [];
        [startPt, currEnd] = edge_endpoints(curr, sampleE);
        loop = curr;
        closed = norm(startPt - currEnd) <= tol;
        while ~isempty(remaining) && ~closed
            idx = find_next(remaining, currEnd, sampleE, tol);
            if isempty(idx)
                break;
            end
            nextEdge = remaining{idx};
            remaining(idx) = [];
            loop(end+1) = nextEdge; %#ok<AGROW>
            [~, currEnd] = edge_endpoints(nextEdge, sampleE);
            closed = norm(currEnd - startPt) <= tol;
        end
        if ~closed
            warning('未能闭合一个环，剩余边可能不连续');
        end
        loops{end+1} = [loop(:)]; %#ok<AGROW> % 转为 struct 数组
    end
end

function idx = find_next(candidates, currEnd, sampleE, tol)
    idx = [];
    for i = 1:numel(candidates)
        [sPt, ~] = edge_endpoints(candidates{i}, sampleE);
        if norm(sPt - currEnd) <= tol
            idx = i;
            return;
        end
    end
end

function [sPt, ePt] = edge_endpoints(e, sampleE)
    sPt = [NaN, NaN]; ePt = [NaN, NaN];
    edgeId = e.edgeId + 1;
    if edgeId > numel(sampleE) || isempty(sampleE{edgeId})
        return;
    end
    pts = sampleE{edgeId};
    if e.dir == -1
        pts = flipud(pts);
    end
    sPt = pts(1, :);
    ePt = pts(end, :);
end
