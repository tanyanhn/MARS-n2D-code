function [V, E, S, F] = parse_geometry_json(filepath)
% PARSE_GEOMETRY_JSON 读取几何 JSON 并返回 V/E/S/F
% V : 2 x m 顶点矩阵
% E : 1 x n Cell，每个单元为 struct(ctrl 8xk, k, v1, v2, idx)
% S : 1 x m Cell；非空时为 struct(pairs kx2, withFlag kx4)
% F : 颜色聚合的边界结构（struct 数组）

    if ~isfile(filepath)
        error('文件不存在: %s', filepath);
    end

    raw_text = fileread(filepath);
    try
        data = jsondecode(raw_text);
    catch ME
        error('JSON 解析失败: %s', ME.message);
    end

    raw_vertices = data.vertex;
    raw_edges = data.edge;
    raw_phases = [];
    if isfield(data, 'phases')
        raw_phases = data.phases;
    end

    num_vertices = numel(raw_vertices);
    num_edges = numel(raw_edges);

    %% 构建 E 和边映射
    E = cell(1, num_edges);
    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'double'); % key -> 1-based edge idx
    edges_meta = repmat(struct('v1', 0, 'v2', 0, 'idx', 0, 'k', 0), 1, num_edges);

    for i = 1:num_edges
        item = raw_edges(i);

        if item.v1 == item.v2
            warning('闭合边检测：edgeId=%d (idx=%d) v1==v2==%d', i-1, item.idx, item.v1);
        end

        segments = normalize_bezier(item.bezier);
        num_segs = size(segments, 1);
        ctrl_mat = zeros(8, num_segs);
        for s = 1:num_segs
            pts = squeeze(segments(s, :, :)); % 4 x 2
            ctrl_mat(1:4, s) = pts(:, 1);
            ctrl_mat(5:8, s) = pts(:, 2);
        end

        E{i} = struct( ...
            'ctrl', ctrl_mat, ...
            'k', num_segs, ...
            'v1', item.v1, ...
            'v2', item.v2, ...
            'idx', item.idx ...
        );

        key = edge_key(item.v1, item.v2, item.idx);
        edge_map(key) = i; % 1-based
        edges_meta(i).v1 = item.v1;
        edges_meta(i).v2 = item.v2;
        edges_meta(i).idx = item.idx;
        edges_meta(i).k = num_segs;
    end

    %% 构建 V 和 S
    V = zeros(2, num_vertices);
    S = cell(1, num_vertices);

    for i = 1:num_vertices
        v_item = raw_vertices(i);
        V(1, i) = v_item.p(1);
        V(2, i) = v_item.p(2);

        raw_smooth = v_item.smoothnessIndicator;
        if isempty(raw_smooth)
            S{i} = [];
            continue;
        end

        pairs_struct = normalize_pairs(raw_smooth);
        if isempty(pairs_struct)
            S{i} = [];
            continue;
        end

        num_pairs = size(pairs_struct, 1);
        with_flag = zeros(num_pairs, 4);
        valid_count = 0;
        current_vertex_idx = i - 1; % 0-based

        for r = 1:num_pairs
            ref1 = pairs_struct(r, 1);
            ref2 = pairs_struct(r, 2);
            try
                info1 = get_edge_info(ref1, current_vertex_idx, edge_map);
                info2 = get_edge_info(ref2, current_vertex_idx, edge_map);
                valid_count = valid_count + 1;
                with_flag(valid_count, :) = [info1(1), info1(2), info2(1), info2(2)];
            catch ME
                warning('平滑对解析失败: 顶点 %d, 行 %d, %s', current_vertex_idx, r, ME.message);
            end
        end

        if valid_count == 0
            S{i} = [];
        else
            with_flag = with_flag(1:valid_count, :);
            S{i} = struct( ...
                'pairs', with_flag(:, [1, 3]), ...
                'withFlag', with_flag ...
            );
        end
    end

    %% 构建 F（按颜色聚合）
    F = build_F(raw_phases, edge_map, edges_meta);
end

%% 辅助函数
function key = edge_key(v1, v2, idx)
    key = sprintf('%d_%d_%d', v1, v2, idx);
end

function segments = normalize_bezier(bz_data)
% 输出: segments (k x 4 x 2)
    segments = [];
    if isnumeric(bz_data)
        dims = size(bz_data);
        if numel(dims) == 2
            seg = ensure_4x2(bz_data);
            segments = reshape(seg, 1, 4, 2);
        else
            segments = bz_data;
            if size(segments, 2) == 2 && size(segments, 3) == 4
                segments = permute(segments, [1, 3, 2]); % k x 4 x 2
            end
        end
    elseif iscell(bz_data)
        k = numel(bz_data);
        segments = zeros(k, 4, 2);
        for i = 1:k
            segments(i, :, :) = ensure_4x2(bz_data{i});
        end
    else
        error('无法识别的 bezier 数据类型');
    end
end

function mat = ensure_4x2(seg)
    mat = seg;
    dims = size(mat);
    if isequal(dims, [2, 4])
        mat = mat';
    elseif ~isequal(dims, [4, 2])
        mat = reshape(mat, 4, 2);
    end
end

function pairs_struct = normalize_pairs(raw_smooth)
    pairs_struct = [];
    if isstruct(raw_smooth)
        if numel(raw_smooth) == 2
            pairs_struct = reshape(raw_smooth, 1, 2);
        else
            dims = size(raw_smooth);
            if numel(dims) == 2 && dims(2) == 2
                pairs_struct = raw_smooth;
            end
        end
    elseif iscell(raw_smooth)
        try
            tmp = [raw_smooth{:}];
            if mod(numel(tmp), 2) == 0
                pairs_struct = reshape(tmp, [], 2);
            end
        catch
            pairs_struct = [];
        end
    end
end

function info = get_edge_info(edge_ref, current_vertex_idx, edge_map)
    key = edge_key(edge_ref.v1, edge_ref.v2, edge_ref.idx);
    if ~isKey(edge_map, key)
        error('缺少边: %s', key);
    end
    edge_idx_1based = edge_map(key);
    if current_vertex_idx == edge_ref.v1
        flag = 0;
    elseif current_vertex_idx == edge_ref.v2
        flag = 1;
    else
        error('顶点 %d 不是边 %s 的端点', current_vertex_idx, key);
    end
    info = [edge_idx_1based - 1, flag]; % 返回 0-based edgeId
end

function F = build_F(raw_phases, edge_map, edges_meta)
    if isempty(raw_phases)
        F = [];
        return;
    end

    color_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(raw_phases)
        ph = raw_phases(i);
        verts = ph.boundaryEdge;
        idxs = ph.edgeIdx;
        color = ph.color;
        ckey = sprintf('%d_%d_%d', color(1), color(2), color(3));
        if ~isKey(color_map, ckey)
            color_map(ckey) = init_color_struct(color);
        end
        cs = color_map(ckey);

        cs.vertexSet = union(cs.vertexSet, verts);

        edgeLoop = [];
        num_seg = numel(verts) - 1;
        for s = 1:num_seg
            vA = verts(s);
            vB = verts(s + 1);
            idxVal = idxs(min(s, numel(idxs)));
            [edgeId, dir] = resolve_edge(vA, vB, idxVal, edge_map, edges_meta);
            edgeLoop = [edgeLoop; struct( ... %#ok<AGROW>
                'edgeId', edgeId, ...
                'dir', dir, ...
                'segCount', edges_meta(edgeId + 1).k ...
            )];

            cs.edgeSet = add_edge_set(cs.edgeSet, edgeId, edges_meta(edgeId + 1));
        end

        b = struct( ...
            'vertexLoop', verts, ...
            'edgeLoop', edgeLoop, ...
            'edgeIdxRaw', idxs ...
        );
        cs.boundaries{end + 1} = b;
        color_map(ckey) = cs;
    end

    % 输出 struct 数组
    keys = color_map.keys;
    F = repmat(struct('color', [], 'vertexSet', [], 'edgeSet', [], 'boundaries', []), 1, numel(keys));
    for i = 1:numel(keys)
        cs = color_map(keys{i});
        cs.vertexSet = unique(cs.vertexSet);
        cs.edgeSet = unique_edgeSet(cs.edgeSet);
        F(i) = cs;
    end
end

function cs = init_color_struct(color)
    cs = struct( ...
        'color', color, ...
        'vertexSet', [], ...
        'edgeSet', [], ...
        'boundaries', {{}});
end

function [edgeId, dir] = resolve_edge(vA, vB, idxVal, edge_map, edges_meta)
    key_forward = edge_key(vA, vB, idxVal);
    key_reverse = edge_key(vB, vA, idxVal);
    if isKey(edge_map, key_forward)
        edgeId = edge_map(key_forward) - 1; % 0-based
        dir = +1;
    elseif isKey(edge_map, key_reverse)
        edgeId = edge_map(key_reverse) - 1;
        dir = -1;
    else
        error('边缺失: (%d,%d,idx=%d)', vA, vB, idxVal);
    end
    % 保底检查元数据
    meta = edges_meta(edgeId + 1);
    if meta.idx ~= idxVal
        warning('edgeId=%d idx 不匹配，期望 %d, 实际 %d', edgeId, idxVal, meta.idx);
    end
end

function edgeSet = add_edge_set(edgeSet, edgeId, meta)
    entry = [edgeId, meta.v1, meta.v2, meta.idx, meta.k];
    edgeSet = [edgeSet; entry];
end

function edgeSet = unique_edgeSet(edgeSet)
    if isempty(edgeSet)
        return;
    end
    [~, ia] = unique(edgeSet(:, 1), 'stable'); % 按 edgeId 去重
    edgeSet = edgeSet(ia, :);
end
