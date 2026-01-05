clear; close all; clc;
addpath ../

% 数据文件路径（可切换为 Raccoon / Pig 等）
name = '../../data1/Pig';
removeV = [];
% removeV = [0, 1, 2, 3]; % Raccoon boundary vertex
removeV = [2, 9, 36, 40]; % Pig boundary vertex

% 解析 JSON
[V, E, S, F] = parse_geometry_json(name + ".json");
[V, E, S, F, maps] = prune_geometry(V, E, S, F, removeV);
rect = [.2 .8 .2 .8];
[V, E, S] = normalize_geometry(V, E, S, rect, pi);

% 构建光滑曲线集合
opts = struct; % 可设置 opts.maxStep / opts.minSamples
curves = constructSmoothCurves(V, E, S, opts);

result = build_yinset(F, V, curves, E);

[curves, splineE] = interpolate_curves_spline(curves, E);
% print_curve(curves);
filepath = name + ".inputbak";
write_geometry_bin(filepath, V, E, S, curves, splineE, result);

% 可视化（已抽离为独立函数）
f1 = figure('Name', 'Smooth Curves Visualization');
axis(rect);
visualize_curves(curves, E, V);

% figure
% fillIdx = [12:12]; % 例: [1 3] 仅渲染 F(1), F(3) 的区域；空数组则不按 F 上色
% plot_boundaries(F, sampleE, fillIdx);
% targetIdx = [3];
% plot_subtree(boundaries, parent, targetIdx, sampleE);
% axis([.1 .9 .1 .9]);


% 收集每个 YinSet 的颜色到 cell
% yinColors = cell(1, numel(result.YinSet));
% for yi = 1:numel(result.YinSet)
%     yinColors{yi} = result.YinSet{yi}.color;
% end
% fprintf('yinColors = {');
% for yi = 1:numel(yinColors)
%     c = yinColors{yi};
%     fprintf('[%d %d %d]', c(1), c(2), c(3));
%     if yi ~= numel(yinColors), fprintf(', '); end
% end
% fprintf('};\n');



function [] = print_curve(curves)
    for i = 1:numel(curves)
        spline = curves{i}.spline;
        xsp = spline.xpp;
        ysp = spline.ypp;
        t = spline.breaks;
        disp(["i = ", int2str(i), "\n"]);
        print_pp_derivatives(xsp, t);
        print_pp_derivatives(ysp, t);
    end
end