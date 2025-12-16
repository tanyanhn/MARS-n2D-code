function plot_subtree(boundaries, parent, targetIdx, sampleE)
% PLOT_SUBTREE 绘制指定节点及其子树对应的边界曲线
% 输入:
%   boundaries: build_boundaries_from_F 的输出（含 edges/vertexLoop）
%   parent: build_containment_tree 的输出
%   targetIdx: 要绘制的根节点下标 (1-based)
%   sampleE: split_curve_samples 的输出 (每条边的采样点，方向 v1->v2)
%   V: 2xM 顶点坐标（仅在缺采样时作为回退）

    nodes = collect_subtree(parent, targetIdx);
    colors = lines(numel(nodes));

    for i = 1:numel(nodes)
        newBoundaries{i} = boundaries{nodes(i)};
    end

    % for k = 1:numel(nodes)
    %     bi = nodes(k);
        plot_bound(newBoundaries, sampleE);
    % end

end
