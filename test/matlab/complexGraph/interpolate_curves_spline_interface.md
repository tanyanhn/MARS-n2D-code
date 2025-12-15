# `interpolate_curves_spline` 接口说明

使用已有的五次样条插值例程（`quinticSpline.m`(周期) / `construct_spline`(B样条) 思路）对 `constructSmoothCurves` 产生的 `curves` 逐条插值，并将得到的全局样条按 `E` 的边结构切分回每条逻辑边。

## 函数签名
```matlab
[curvesOut, splineE] = interpolate_curves_spline(curves, E)
```

## 输入
- `curves`: `constructSmoothCurves` 的输出，需包含：
  - `samples` (`N x 2`): 曲线采样点，顺序为曲线行进方向。
  - `sampleBreaks`: 每段 Bezier 的起始采样索引（与 `edges.segCount` 一一对应）。
  - `edges`: 有序列表，元素含 `edgeId` (0-based)、`dir` (±1)、`segCount`。
  - `isClosed`: 闭合标记。
- `E`: 原始边列表（1 x n Cell），用于获取边计数和方向语义（v1→v2）。

## 输出
- `curvesOut`: 在原有 `curves` 基础上附加样条信息：
  - `spline.xpp`, `spline.ypp`: 分别为 x/y 方向的 `pp` 样条（来自 `quinticSpline`）。
  - `spline.breaks`: 样条节点（0..1）。
- `splineE`: 1 x n Cell；`splineE{i}` 对应 `E{i}`，包含：
  - `curveId`: 源 `curves` 下标。
  - `dir`: 相对于 `E` 的方向，`+1` 表示 `v1->v2`，`-1` 表示反向。
  - `t0`, `t1`: 在 `curves{curveId}.spline` 参数域中的子区间（已按 `dir` 方向排列）。
  - `xpp`, `ypp`: 对应曲线的全局样条（评估时需用 `t = t0 + (t1 - t0)*u`，`u ∈ [0,1]`；若 `dir=-1` 则需反向映射）。
  - （二进制输出时，为避免冗余，可仅存 `curveId, dir, t0, t1`，读取侧用 `curveId` 查样条系数。）

## 处理策略
- 插值：直接调用 `quinticSpline` 对 `curves{i}.samples'` 进行五次样条插值；闭合/开曲线均可使用。若需要严格周期插值，可替换为 `quintic_periodic_spline` 思路（需返回 `pp`）。
- 参数化：使用采样点弦长归一化到 `[0,1]` 作为样条参数；首点为 0，末点为 1。零长度的曲线会给出常值参数并跳过插值。
- 切分：利用 `sampleBreaks` 和 `edges.segCount` 找到每条逻辑边在采样序列中的起止索引，转换为参数区间 `[t0, t1]`；若 `dir=-1` 则交换区间端点以保持与 `E` (v1→v2) 一致。
- 多次引用：如果同一 `edgeId` 在不同曲线/方向重复出现，后写覆盖前写（可按需修改为冲突检查）。

## 评估提示
- 计算边上点：给定局部参数 `u∈[0,1]`，先映射到全局参数 `t = t0 + (t1 - t0)*u`（`dir=-1` 则 `t0 > t1`），再用 `ppval(xpp,t)` / `ppval(ypp,t)` 获取坐标。
- 如果需要保持首尾点严格一致，可在生成样条前确保 `samples` 首尾重合（对闭合曲线）。***
