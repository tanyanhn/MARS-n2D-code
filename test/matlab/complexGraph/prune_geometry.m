function [V2, E2, S2, F2, maps] = prune_geometry(V, E, S, F, removeVerts)
% PRUNE_GEOMETRY 删除指定顶点及其关联边，调整 S/F 下标
% 输入:
%   V, E, S, F: parse_geometry_json 的输出（0-based 语义）
%   removeVerts: 要删除的顶点下标数组（0-based）
% 输出:
%   V2, E2, S2, F2: 过滤后的结构
%   maps: struct，包含 oldV2new, oldE2new 映射

    if nargin < 5
        error('缺少 removeVerts');
    end
    if isempty(removeVerts)
        V2 = V; E2 = E; S2 = S; F2= F; maps = {};
        return;
    end
    removeVerts = unique(removeVerts(:)');
    numV = size(V, 2);
    keepV = true(1, numV);
    keepV(removeVerts + 1) = false;

    % 顶点映射 old -> new
    oldV2new = -ones(1, numV);
    V2 = [];
    cntV = 0;
    for i = 1:numV
        if ~keepV(i), continue; end
        cntV = cntV + 1;
        oldV2new(i) = cntV - 1; % 0-based
        V2(:, cntV) = V(:, i);
    end

    % 边过滤
    numE = numel(E);
    keepE = false(1, numE);
    oldE2new = -ones(1, numE);
    E2 = {};
    cntE = 0;
    for i = 1:numE
        e = E{i};
        if ~keepV(e.v1 + 1) || ~keepV(e.v2 + 1)
            continue;
        end
        cntE = cntE + 1;
        keepE(i) = true;
        oldE2new(i) = cntE - 1; % 0-based
        e.v1 = oldV2new(e.v1 + 1);
        e.v2 = oldV2new(e.v2 + 1);
        E2{cntE} = e; %#ok<AGROW>
    end

    % S 过滤
    S2 = cell(1, cntV);
    for oldV = 0:(numV - 1)
        newVid = oldV2new(oldV + 1);
        if newVid < 0, continue; end
        entry = S{oldV + 1};
        if isempty(entry)
            S2{newVid + 1} = [];
            continue;
        end
        if isstruct(entry)
            wf = entry.withFlag;
        else
            wf = entry;
        end
        % 过滤被删边
        mask = keepE(wf(:, 1) + 1) & keepE(wf(:, 3) + 1);
        wf = wf(mask, :);
        if isempty(wf)
            S2{newVid + 1} = [];
            continue;
        end
        % 重映射 edgeId
        wf(:, 1) = oldE2new(wf(:, 1) + 1);
        wf(:, 3) = oldE2new(wf(:, 3) + 1);
        S2{newVid + 1} = struct( ...
            'pairs', wf(:, [1, 3]), ...
            'withFlag', wf ...
        );
    end

    % F 过滤
    F2 = [];
    if ~isempty(F)
        F2 = repmat(struct('color', [], 'vertexSet', [], 'edgeSet', [], 'boundaries', []), 1, 0);
        for fi = 1:numel(F)
            cs = F(fi);
            newBoundaries = {};
            newVset = [];
            newEset = [];
            for bi = 1:numel(cs.boundaries)
                b = cs.boundaries{bi};
                vloop = b.vertexLoop;
                if any(~keepV(vloop + 1))
                    continue;
                end
                eloop = b.edgeLoop;
                edgeIds = [eloop.edgeId] + 1;
                if any(~keepE(edgeIds))
                    continue;
                end
                % remap
                vloop = oldV2new(vloop + 1);
                for k = 1:numel(eloop)
                    eloop(k).edgeId = oldE2new(eloop(k).edgeId + 1);
                    % segCount 可保持原值
                end
                newBoundaries{end+1} = struct( ... %#ok<AGROW>
                    'vertexLoop', vloop, ...
                    'edgeLoop', eloop, ...
                    'edgeIdxRaw', b.edgeIdxRaw ...
                );
                newVset = union(newVset, vloop);
                for k = 1:numel(eloop)
                    eid = eloop(k).edgeId;
                    if eid >= 0
                        eMeta = E2{eid + 1};
                        newEset = [newEset; [eid, eMeta.v1, eMeta.v2, eMeta.idx, eMeta.k]]; %#ok<AGROW>
                    end
                end
            end
            if isempty(newBoundaries)
                continue;
            end
            cs2 = struct( ...
                'color', cs.color, ...
                'vertexSet', unique(newVset), ...
                'edgeSet', unique_edgeSet(newEset), ...
                'boundaries', {newBoundaries} ...
            );
            F2(end+1) = cs2; %#ok<AGROW>
        end
    end

    maps = struct('oldV2new', oldV2new, 'oldE2new', oldE2new);
end

function edgeSet = unique_edgeSet(edgeSet)
    if isempty(edgeSet), return; end
    [~, ia] = unique(edgeSet(:, 1), 'stable');
    edgeSet = edgeSet(ia, :);
end
