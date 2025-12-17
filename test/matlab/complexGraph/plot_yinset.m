function plot_yinset(result, yinIdx)
% PLOT_YINSET 绘制指定 YinSet 的所有边界
% 输入:
%   result: build_yinset 的返回结构（含 newBoundaries, YinSet）
%   sampleE: split_curve_samples 输出
%   yinIdx: 要绘制的 YinSet 下标（可为标量或向量，1-based）
    
      
    if nargin < 2 || isempty(yinIdx)
        yinIdx = 1:numel(result.YinSet);
    end
    sampleE = result.sampleE;

    holdState = ishold;
    hold on; axis equal; grid on;

    colors = lines(numel(yinIdx));
    for k = 1:numel(yinIdx)
        yi = yinIdx(k);
        ys = result.YinSet{yi};
        col = colors(k, :);
        % 聚合该 YinSet 下所有边界索引
        idxList = [];
        for ci = 1:numel(ys.components)
            idxList = [idxList, ys.components{ci}.indices]; %#ok<AGROW>
        end
        idxList = unique(idxList);
        for bi = idxList
            plot_bound(result.newBoundaries(bi), sampleE, true);
        end
    end

    if ~holdState
        hold off;
    end
end
