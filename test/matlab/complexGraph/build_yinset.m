function result = build_yinset(F, V, curves, E, tol)
% BUILD_YINSET 按 yinset_pipeline 生成闭合边界/组件/YinSet
% 输入:
%   F, V, curves, E: 解析/构造管线的输出
%   tol: 可选，端点匹配容差，默认 1e-7
% 输出 result 结构体:
%   boundaries     : 所有原始闭合曲线（build_boundaries_from_F）
%   sampleE        : 按边拆分的采样（split_curve_samples）
%   containMat     : 包含矩阵
%   parent         : 直接父映射
%   newBoundaries  : 去重闭合后的边界集合
%   components     : 外边界+洞的组件列表（含 color, indices）
%   YinSet         : 按颜色聚合的组件

    if nargin < 5 || isempty(tol)
        tol = 1e-7;
    end

    result = struct();

    % 原始边界
    boundaries = build_boundaries_from_F(F);
    sampleE = split_curve_samples(curves, E);
    containMat = boundary_contains(boundaries, V, sampleE);
    parent = build_containment_tree(containMat);

    result.boundaries = boundaries;
    result.sampleE = sampleE;
    result.containMat = containMat;
    result.parent = parent;

    % 重建闭合边界集合与组件
    newBoundaries = {};
    components = {};
    for i = 1:numel(boundaries)
        nb = merge_boundaries(boundaries, parent, i, sampleE, tol);
        if i >=1
            f = figure;
            plot_bound(nb, sampleE);
            axis([.1 .9 .1 .9]);
            close(f);
        end
        baseIdx = numel(newBoundaries);
        for k = 1:numel(nb)
            newBoundaries{end+1} = nb{k}; %#ok<AGROW>
        end
        components{end+1} = struct( ... %#ok<AGROW>
            'color', boundaries{i}.color, ...
            'indices', baseIdx + (1:numel(nb)) ...
        );
    end

    result.newBoundaries = newBoundaries;
    result.components = components;

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

    result.YinSet = YinSet;
end
