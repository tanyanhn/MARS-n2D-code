function plot_bound(boundaries, sampleE)

holdState = ishold;
hold on;
axis equal; grid on;
colors = lines(numel(boundaries));
for k = 1:numel(boundaries)
    boundary = boundaries{k};
    edgeLoop = boundary.edges;
    ptsAll = [];
    for eidx = 1:numel(edgeLoop)
        edgeId = edgeLoop(eidx).edgeId + 1;
        if edgeId > numel(sampleE) || isempty(sampleE{edgeId})
            continue;
        end
        pts = sampleE{edgeId};
        if edgeLoop(eidx).dir == -1
            pts = flipud(pts);
        end
        if isempty(ptsAll)
            ptsAll = pts;
        else
            if norm(ptsAll(end, :) - pts(1, :)) < 1e-9
                ptsAll = [ptsAll; pts(2:end, :)]; %#ok<AGROW>
            else
                ptsAll = [ptsAll; pts]; %#ok<AGROW>
            end
        end
    end
    if nargin < 4
        col = 'k';
    end
    plot(ptsAll(:, 1), ptsAll(:, 2), '-', 'Color', colors(k, :), 'LineWidth', 1.5);
end

if ~holdState
    hold off;
end
end
