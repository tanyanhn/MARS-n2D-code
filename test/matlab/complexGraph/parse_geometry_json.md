# `parse_geometry_json` 接口说明（扩展 F）

基于 `test/data1/jsonInput.md` 的 JSON 结构约定，对 `test/matlab/parse_geometry_json.m`（当前实现入口为 `parse_geometry_json_v2`）的输入 / 输出做接口定义，并在保持 `V、E、S` 输出的前提下，新增 `F` 用于按颜色聚合同一大结构下的边界与顶点信息。本版将 `E` 设计为“每条逻辑边一个 Cell”，并在单元内展开该边的多段 Bezier。

## 函数签名
```matlab
[V, E, S, F] = parse_geometry_json(filepath)
```
- 兼容旧调用：若只需要 `V, E, S`，可以按原三输出调用；扩展场景下按四输出取 `F`。
- 0-based 索引语义仅体现在返回值的“ID”含义上（矩阵 / Cell 的物理下标仍是 MATLAB 1-based）。

## 输入
- `filepath`：几何 JSON 文件路径，遵循 `jsonInput.md` 中的 `vertex / edge / phases` 结构。

## 输出
### 1) `V` — 顶点矩阵
- 形状：`2 x m`。
- 第 1 行为 x，第 2 行为 y。列下标（0-based 语义）与 JSON 中 `vertex` 的顺序一致。

### 2) `E` — 曲线单元列表（Cell）
- 类型：`1 x n` Cell，其中 `n` 为 JSON `edge` 的项数。
- 每个单元与一条逻辑边一一对应（同一个 `(v1, v2, idx)`），内部可包含多段 Bezier。
- 建议单元结构（示例）：
  ```matlab
  E{i} = struct( ...
      'ctrl', ctrl_mat, ...   % 8 x k，列为该逻辑边的第 j 段 Bezier
      'k',    k, ...          % 段数，可选
      'v1',   v1, ...         % 0-based
      'v2',   v2, ...         % 0-based
      'idx',  idx ...         % 0-based
  );
  % ctrl_mat 行 1-4: P0.x..P3.x，行 5-8: P0.y..P3.y
  % 段顺序遵循 v1 -> v2，列 1 与 v1 相连，列 k 与 v2 相连。
  ```
- 逻辑边的 ID = `E` 的单元下标（0-based 语义，对应 MATLAB 下标 `i+1`）。

### 3) `S` — 光滑性约束
- 类型：`1 x m` Cell；`S{i+1}` 对应顶点 `i`。
- 数据结构（非空时为 struct，空约束则为 `[]`）：
    - `pairs`: `k x 2` 矩阵，每行一对光滑连接的逻辑边 ID `[e_i, e_j]`（0-based，指向 `E{e+1}`）。
    - `withFlag`: `k x 4` 矩阵，显式保存端标记 `[e_i, flag_i, e_j, flag_j]`，`flag = 0` 表示 `v1` 端，`flag = 1` 表示 `v2` 端。
- 端点映射规则与 JSON 的 `smoothnessIndicator` 对应：
    - 若当前顶点等于 `E{e+1}.v1`，使用该单元 `ctrl` 的 **第 1 列** 作为连接端（`flag = 0`）。
    - 若当前顶点等于 `E{e+1}.v2`，使用该单元 `ctrl` 的 **第 k 列** 作为连接端（`flag = 1`）。
- 允许 `e_i == e_j`（自光滑）。

### 4) `F` — 颜色聚合的面/边界结构（新增）
- 类型：建议为 struct 数组或 Cell，每个元素对应一种颜色（一个“大结构”）。
- 字段约定（推荐）：
    - `color`: `[r, g, b]`，与 `phases.color` 完全一致。
    - `vertexSet`: 本颜色下所有出现的顶点下标集合（0-based 语义，源自 `phases.boundaryEdge`）。
    - `edgeSet`: 本颜色下涉及的逻辑边集合，元素形如 `[edgeId, v1, v2, idx, k]`（`edgeId` 为 `E` 的单元下标，`k` 为段数，可选）。
    - `boundaries`: Cell 数组，每个单元代表一个闭合边界（对应单个 `phases` 条目）：
        - `vertexLoop`: `boundaryEdge` 原样存储（顶点序列，通常首尾重复闭合）。
        - `edgeLoop`: 与 `vertexLoop` 的相邻顶点对一一对应，元素包含：
            - `edgeId`: 逻辑边 ID（指向 `E{edgeId+1}`，包含 `v1, v2, idx` 元数据）。
            - `dir`: `+1` 按 `v1->v2` 使用 `ctrl`，`-1` 按 `v2->v1`（消费时可反转列顺序）。
            - `segCount`: 该逻辑边的段数（可直接取 `E{edgeId+1}.k` 或 `size(ctrl,2)`）。
        - `edgeIdxRaw`: 保留 JSON 中的 `edgeIdx` 序列以便调试。
- 一个 `F` 元素可以包含多个闭合边界（同色多个 `phases` 合并），但边界内部相互独立存储。

## 构建要点（实现建议）
- 解析顺序：`vertex` -> `edge` -> `phases`。先构建 `V`、`E` 及映射 `edge_idx_map((v1,v2,idx)) = edgeId`，再用它们驱动 `S` 与 `F`。
- `S` 和 `F` 都依赖逻辑边到 `E` 单元的映射；务必统一 `0-based` 语义。
- 若发现闭合边（`v1 == v2`），建议输出警告，便于定位数据质量问题。
- 方向处理：`phases.boundaryEdge` 给出行进方向。若相邻顶点顺序与逻辑边定义相反，`dir = -1`，消费时可按需反转 Bezier 顺序。
- 兼容性：若调用侧只关心旧输出，可忽略/不请求第 4 个返回值；新增 `F` 不改变 `V, E, S` 的含义与形状。
