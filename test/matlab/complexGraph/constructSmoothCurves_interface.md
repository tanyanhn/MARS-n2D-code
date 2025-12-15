# `constructSmoothCurves` 接口说明

基于 `parse_geometry_json` 的输出 (`V, E, S`)，将曲线按光滑约束组合成有序子集合，并在每个集合上生成采样点。

## 函数签名
```matlab
curves = constructSmoothCurves(V, E, S, opts)
```
- `opts` 可选，结构体；常用字段：
  - `maxStep` (默认 = 所有 Bezier 段控制折线长度的最小值): 相邻采样点的距离上界（欧氏距离）；如未传入则以最短段长度为上界，保证至少两点。
  - `minSamples` (默认 `2`): 每段 Bezier 至少采样的点数（包含端点）。

## 输入
- `V` (`2 x m`): 顶点矩阵，列下标的 0-based 语义与 JSON `vertex` 一致。
- `E` (`1 x n` Cell): 每个单元为 struct `{ctrl 8xk, k, v1, v2, idx}`，表示一条逻辑边及其多段 Bezier（列顺序为 v1→v2）。
- `S` (`1 x m` Cell):
  - 空约束用 `[]`。
  - 非空时为 struct:
    - `pairs` (`p x 2`): 成对光滑连接的逻辑边 ID（0-based，对应 `E{edgeId+1}`）。
    - `withFlag` (`p x 4`): `[edgeA, flagA, edgeB, flagB]`，`flag=0` 表示该边的 v1 端，`flag=1` 表示 v2 端。

## 子集合构造规则
- 基于 `withFlag`，每条边的每个端点最多与一条边光滑连接（可与自己连接形成闭环）。
- 将光滑连接视为无交叉的链/环，按如下方式生成有序集合：
  - 若某条边在某端无光滑连接，则从该端开始/结束，形成开放链。
  - 若所有端点均有光滑连接，则形成闭合环。
  - 构建过程中，若连接方向与边定义 (v1→v2) 相反，需要反转该边的 `ctrl` 列顺序以保证相邻端点一致。
- 每个集合附带一个标记：
  - `isClosed = true` 表示周期光滑曲线（首尾相连且首尾也在 `withFlag` 中出现）。
  - `isClosed = false` 表示普通曲线段（可能只有一条边）。

## 采样规则
- 采样覆盖集合中的所有 Bezier 段；每段至少包含其两个端点。
- 若设置 `opts.maxStep` 为有限值，则沿参数方向细分该段，确保相邻采样点间的弦长不超过该上界。
- 各段的采样点按照集合的有序方向串接，若集合闭合且首尾重合，则保留首点作为唯一重复点或根据实现选择去重。

## 输出 `curves`
- Cell 数组，每个单元对应一个光滑子集合，推荐结构：
  - `edges`: 有序列表，元素为 struct `{edgeId, dir, segCount}`，`dir=+1` 表示按 `E{edgeId+1}` 的方向，`dir=-1` 表示反向。
  - `ctrl`: 串接后的 `8 x K` Bezier 控制矩阵（已按方向调整列顺序）。
  - `isClosed`: 布尔，周期光滑标记。
  - `samples`: `N x 2` 采样点坐标，按照曲线行进顺序排列。
  - （可选）`sampleBreaks`: 采样点对应到每段 Bezier 的起始索引，便于做参数映射。

## 错误与警告
- 若 `withFlag` 引用不存在的边、或顶点不匹配、或检测到违反“一端最多一条光滑连接”的情况，应报错或警告（实现者可选策略）。
