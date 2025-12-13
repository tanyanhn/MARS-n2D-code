clear; close all; clc;

% 数据文件路径（可切换为 raccoonNew.json / pigNew.json 等）
file = '../data1/raccoonNew.json';

% 解析 JSON
[V, E, S, ~] = parse_geometry_json(file);

% 构建光滑曲线集合
opts = struct; % 可设置 opts.maxStep / opts.minSamples
curves = constructSmoothCurves(V, E, S, opts);

% 可视化（已抽离为独立函数）
figure('Name', 'Smooth Curves Visualization');
visualize_curves(curves, E, V);
