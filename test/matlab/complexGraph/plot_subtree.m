function plot_subtree(boundaries, parent, targetIdx, sampleE, V)
% PLOT_SUBTREE 绘制指定节点及其子树对应的边界曲线
% 输入:
%   boundaries: build_boundaries_from_F 的输出（含 edges/vertexLoop）
%   parent: build_containment_tree 的输出
%   targetIdx: 要绘制的根节点下标 (1-based)
%   sampleE: split_curve_samples 的输出 (每条边的采样点，方向 v1->v2)
%   V: 2xM 顶点坐标（仅在缺采样时作为回退）

    if targetIdx < 1 || targetIdx > numel(boundaries)
        error('targetIdx 超界');
    end
    % 收集子树节点
    nodes = collect_subtree(parent, targetIdx);
    colors = lines(numel(nodes));

    holdState = ishold;
    hold on;
    axis equal; grid on;

    for k = 1:numel(nodes)
        bi = nodes(k);
        col = colors(k, :);
        edgeLoop = boundaries{bi}.edges;
        ptsAll = [];
        for eidx = 1:numel(edgeLoop)
            edgeId = edgeLoop(eidx).edgeId + 1;
            if edgeId > numel(sampleE) || isempty(sampleE{edgeId})
                continue;
            end
            pts = sampleE{edgeId};
            if edgeLoop(eidx).dir == -1
                pts = flipud(pts);
            end
            if isempty(ptsAll)
                ptsAll = pts;
            else
                if norm(ptsAll(end, :) - pts(1, :)) < 1e-9
                    ptsAll = [ptsAll; pts(2:end, :)]; %#ok<AGROW>
                else
                    ptsAll = [ptsAll; pts]; %#ok<AGROW>
                end
            end
        end
        if isempty(ptsAll)
            loop = boundaries{bi}.vertexLoop;
            ptsAll = V(:, loop + 1).';
        end
        plot(ptsAll(:, 1), ptsAll(:, 2), '-', 'Color', col, 'LineWidth', 1.5);
    end

    if ~holdState
        hold off;
    end
end

function nodes = collect_subtree(parent, root)
    nodes = root;
    children = find(parent == root);
    for c = children(:).'
        nodes = [nodes, collect_subtree(parent, c)]; %#ok<AGROW>
    end
end
