function boundaries = build_boundaries_from_F(F)
% BUILD_BOUNDARIES_FROM_F 将 F 结构展开为 boundary Cell 列表
% 输入:
%   F: parse_geometry_json 的颜色分组结构，含 fields: color, boundaries.edgeLoop/vertexLoop
% 输出:
%   boundaries: 1 x NB cell，每个单元为 struct:
%       color: [r g b] (uint8)
%       edges: struct 数组，元素含 edgeId (0-based), dir (+1/-1)
%       vertexLoop: 顶点序列 (0-based)

    boundaries = {};
    if isempty(F)
        return;
    end
    for fi = 1:numel(F)
        col = uint8(F(fi).color(:)).';
        bnds = F(fi).boundaries;
        for bi = 1:numel(bnds)
            el = bnds{bi}.edgeLoop;
            edges = arrayfun(@(e) struct('edgeId', e.edgeId, 'dir', e.dir), ...
                             el, 'UniformOutput', true);
            boundaries{end+1} = struct( ... %#ok<AGROW>
                'color', col, ...
                'edges', edges, ...
                'vertexLoop', bnds{bi}.vertexLoop ...
            );
        end
    end
end
