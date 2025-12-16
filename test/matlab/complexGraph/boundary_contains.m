function containMat = boundary_contains(boundaries, V, sampleE)
% BOUNDARY_CONTAINS 计算闭曲线两两之间的包含关系
% 输入:
%   boundaries: build_boundaries_from_F 输出的 Cell，每个含 edges (edgeId,dir)
%   V: 2xM 顶点坐标矩阵（备用）
%   sampleE: 1xN cell，split_curve_samples 的输出，edgeId 对应 E 下标
% 输出:
%   containMat: logical NxN，containMat(i,j)=true 表示 boundary i 包含 boundary j

    n = numel(boundaries);
    containMat = false(n, n);
    if n == 0
        return;
    end

    polys = cell(1, n);
    centers = cell(1, n);
    for i = 1:n
        % 用 sampleE 上的采样点拼接多边形，方向按照 edge.dir 调整
        ptsAll = [];
        ptsSearch = [];
        edgeLoop = boundaries{i}.edges;
        for k = 1:numel(edgeLoop)
            edgeId = edgeLoop(k).edgeId + 1;
            if edgeId > numel(sampleE) || isempty(sampleE{edgeId})
                continue;
            end
            pts = sampleE{edgeId};
            if edgeLoop(k).dir == -1
                pts = flipud(pts);
            end
            ptsSearch = [ptsSearch; pts(1:2, :)];
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
        if isempty(ptsAll)
            coords = V(:, boundaries{i}.vertexLoop + 1);
        else
            coords = ptsAll.';
        end
        try
            p = polyshape(coords(1, :), coords(2, :), 'Simplify', false);
        catch
            % 若构造失败，尝试简化
            p = polyshape(coords(1, :), coords(2, :));
        end
        polys{i} = p;
        centers{i} = ptsSearch(:, :);
    end

    for i = 1:n
        for j = 1:n
            if i == j, continue; end
            if isempty(polys{i}) || isempty(polys{j})
                continue;
            end
            pts = centers{j};
            containMat(i, j) = -1;
            for k = 1:size(pts, 1)
                pt = pts(k, :);
                [in, on] = isinterior(polys{i}, pt(1), pt(2));
                if on 
                    continue;
                end
                containMat(i, j) = in;
            end
            if containMat(i, j) == -1
                error("fail to search boundaries " + int2str(i) + " and " + int2str(j) + " relateion\n");
            end
            % containMat(i, j) = isinterior(polys{i}, pt(1), pt(2));
            if (containMat(i, j) == 1 && containMat(j, i) == 1)
                error("boundary i and j contain each other.");
            end
        end
    end
end
