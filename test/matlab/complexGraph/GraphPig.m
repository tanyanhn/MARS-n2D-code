clear; close all; clc;

% 数据文件路径（可切换为 raccoonV2 / pigNew 等）
name = '../../data1/raccoonV2';
addpath ../

% 解析 JSON
[V, E, S, F] = parse_geometry_json(name + ".json");
% rect = [.2, .1, .8, .9];
rect = [.1, .1, .9, .9];
[V, E, S] = normalize_geometry(V, E, S, rect, pi);

% 构建光滑曲线集合
opts = struct; % 可设置 opts.maxStep / opts.minSamples
curves = constructSmoothCurves(V, E, S, opts);

sampleE = split_curve_samples(curves, E);

boundaries = build_boundaries_from_F(F);
containMat = boundary_contains(boundaries, V, sampleE);
parent = build_containment_tree(containMat);
targetIdx = [1];
plot_subtree(boundaries, parent, targetIdx, sampleE, V);

% [curves, splineE] = interpolate_curves_spline(curves, E);
% filepath = name + '.input';
% write_geometry_bin(filepath, V, E, S, curves, splineE);

% 可视化（已抽离为独立函数）
f1 = figure('Name', 'Smooth Curves Visualization');
axis([.1 .9 .1 .9]);
visualize_curves(curves, E, V);

% figure
% fillIdx = [12:12]; % 例: [1 3] 仅渲染 F(1), F(3) 的区域；空数组则不按 F 上色
% plot_boundaries(F, sampleE, fillIdx);
% axis([.1 .9 .1 .9]);