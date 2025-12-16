# YinSet 生成流程设计（基于 F / boundaries / sampleE）

目标：从 `parse_geometry_json` 的 `F`（颜色分组面片）出发，构造：
1) `boundaries`：所有闭合曲线列表（含颜色、edgeId、dir）。
2) `parent`：闭曲线之间的直接包含关系树。
3) `newBoundaries`：消除重叠边后重建的闭合曲线集合。
4) `components`：每个连通分量（外边界 + 直接洞）在 `newBoundaries` 中的下标。
5) `YinSet`：按颜色聚合的连通分量集合。

## 依赖函数
- `build_boundaries_from_F(F) -> boundaries`
  - 展开 `F` 为闭合曲线列表（fields: `color`, `edges`(edgeId,dir), `vertexLoop`）。
- `split_curve_samples(curves, E) -> sampleE`
  - 将 `constructSmoothCurves` 的采样按边拆分。
- `boundary_contains(boundaries, V, sampleE) -> containMat`
  - 用边采样重建多边形，判断两两包含（面积打破对称包含）。
- `build_containment_tree(containMat) -> parent`
  - 生成直接父节点映射（无父则 0）。
- `merge_boundaries(boundaries, parent, id, sampleE, tol) -> newBoundaries`
  - 以 `id` 为根，聚合其子树边，抵消反向重叠，重建闭合环。

## 数据结构约定
- `boundaries{i}`: struct `{ color[1x3], edges[1xK struct(edgeId,dir)], vertexLoop }`
- `parent`: 1xN 向量，`parent(j)=i` 表示 i 直接包含 j；0 表示无父。
- `newBoundaries`: Cell，同 `boundaries` 格式，edge 已去重、闭合。
- `components`: Cell，每单元为 struct:
  - `color`: [r g b]
  - `indices`: `newBoundaries` 的下标，`indices(1)` 外边界，其余为洞。
- `YinSet`: Cell/struct 数组，按颜色聚合:
  - `color`: [r g b]
  - `components`: 该颜色下的 `components` 下标或直接列表。

## 推荐流程（伪代码）
```matlab
% 输入：F, V, curves, E
boundaries = build_boundaries_from_F(F);
sampleE    = split_curve_samples(curves, E);
containMat = boundary_contains(boundaries, V, sampleE);
parent     = build_containment_tree(containMat);

newBoundaries = {};
components    = {};
for i = 1:numel(boundaries)
    if parent(i) ~= 0, continue; end % 只从根开始
    % 收集根及其直接子洞
    child = find(parent == i);
    idxSet = [i, child(:)'];
    nb = merge_boundaries(boundaries, parent, i, sampleE); % 返回闭合环列表
    baseIdx = numel(newBoundaries);
    for k = 1:numel(nb)
        newBoundaries{end+1} = nb{k}; %#ok<AGROW>
    end
    components{end+1} = struct( ... %#ok<AGROW>
        'color', boundaries{i}.color, ...
        'indices', baseIdx + (1:numel(nb)) ...
    );
end

% 按颜色聚合到 YinSet
YinSet = {};
for c = 1:numel(components)
    col = components{c}.color;
    pos = find(cellfun(@(y) isequal(y.color, col), YinSet), 1);
    if isempty(pos)
        YinSet{end+1} = struct('color', col, 'components', {{components{c}}}); %#ok<AGROW>
    else
        YinSet{pos}.components{end+1} = components{c};
    end
end
```

## 说明
- `merge_boundaries` 内部已自动将洞方向取反、抵消重叠边，并用 `sampleE` 端点重新闭合。
- `components` 以“外边界 + 直接洞”为单位；若存在更深层嵌套，可按需递归多层合并。
- 聚合颜色时默认使用外边界的颜色；若存在颜色不一致的洞，可在合并前做检查/过滤。
