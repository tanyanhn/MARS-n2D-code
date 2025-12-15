function plot_boundaries(F, sampleE, regionIdx)
% PLOT_BOUNDARIES 使用 F.boundaries.edgeId 和 sampleE 绘制区域边界
% 输入:
%   F         parse_geometry_json 的 F 结构
%   sampleE   1xN cell，split_curve_samples 的输出，sampleE{i} 为 edge i 的采样点 (Nx2)，方向为 v1->v2
%   regionIdx 可选，F 的下标数组；缺省或空表示绘制全部

    if nargin < 3 || isempty(regionIdx)
        regionIdx = 1:numel(F);
    end
    regionIdx = regionIdx(regionIdx >= 1 & regionIdx <= numel(F));
    if isempty(regionIdx)
        warning('无可绘制的区域索引');
        return;
    end

    holdState = ishold;
    hold on;
    colors = lines(numel(regionIdx));

    for ii = 1:numel(regionIdx)
        fi = regionIdx(ii);
        faceColor = colors(ii, :);
        boundaries = F(fi).boundaries;
        for bi = 1:numel(boundaries)
            edgeLoop = boundaries{bi}.edgeLoop;
            ptsAll = [];
            for ej = 1:numel(edgeLoop)
                eInfo = edgeLoop(ej);
                edgeId = eInfo.edgeId + 1; % MATLAB 索引
                if edgeId > numel(sampleE) || isempty(sampleE{edgeId})
                    warning('区域 %d 边界 %d 的 edgeId=%d 缺少采样', fi, bi, edgeId - 1);
                    continue;
                end
                pts = sampleE{edgeId};
                if eInfo.dir == -1
                    pts = flipud(pts);
                end
                if isempty(ptsAll)
                    ptsAll = pts;
                else
                    % 避免重复连接点
                    if norm(ptsAll(end, :) - pts(1, :)) < 1e-9
                        ptsAll = [ptsAll; pts(2:end, :)]; %#ok<AGROW>
                    else
                        ptsAll = [ptsAll; pts]; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(ptsAll)
                % f = figure;
                plot(ptsAll(1:end-1, 1), ptsAll(1:end-1, 2), '-', 'Color', faceColor, 'LineWidth', 1.5); hold on;
                plot(ptsAll(1, 1), ptsAll(1, 2), 'o');
                % close(f);
            end
        end
    end

    hold on;
    axis equal;
    grid on;

    if ~holdState
        hold off;
    end
end
