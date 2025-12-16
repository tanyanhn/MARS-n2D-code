function write_geometry_bin(filepath, V, E, S, curves, splineE, yinsetResult)
% WRITE_GEOMETRY_BIN 将 (V,E,S) 及样条数据写出为紧凑二进制 (v4)
% 曲线以 Curve<2,6>::dump 兼容格式存储。
% 格式 (小端):
%   magic[4] = 'GEOB'
%   uint32 version = 4
%   uint32 numV, numE, numCurves, numNB, numYin
%   V: double[2*numV] (先 x 行后 y 行按列写)
%   Edge 元数据 (numE 次):
%       int32 v1, v2, idx
%       int32 curveId   (对应 curves，-1 表示缺失)
%       int8  dir       (+1 v1->v2, -1 反向)
%       double t0, t1   (曲线参数区间，方向对应 v1->v2)
%   Curves (numCurves 次): Curve<2,6>::dump 格式
%   S 光滑性 (numV 次):
%       uint32 pairCount
%       重复 pairCount 次: int32 eA, int8 fA, int32 eB, int8 fB
%   newBoundaries (numNB 次):
%       uint8 color[3]
%       uint32 edgeCount
%       重复 edgeCount 次: int32 edgeId, int8 dir
%   YinSet (numYin 次):
%       uint8 color[3]
%       uint32 idxCount
%       uint32 boundaryIdx[idxCount] (指向 newBoundaries，0-based)

    fid = fopen(filepath, 'w');
    if fid < 0
        error('无法打开文件: %s', filepath);
    end
    c = onCleanup(@() fclose(fid));

    numV = size(V, 2);
    numE = numel(E);
    numCurves = numel(curves);
    if nargin < 7 || isempty(yinsetResult)
        error('缺少 yinsetResult，请先运行 build_yinset');
    end
    numNB = numel(yinsetResult.newBoundaries);
    numYin = numel(yinsetResult.YinSet);

    % header
    fwrite(fid, 'GEOB', 'uint8');
    fwrite(fid, uint32(4), 'uint32');
    fwrite(fid, uint32([numV, numE, numCurves, numNB, numYin]), 'uint32');

    % V
    fwrite(fid, V, 'double');

    % Edge meta + splineE 映射
    for i = 1:numE
        meta = E{i};
        if i <= numel(splineE) && ~isempty(splineE{i})
            se = splineE{i};
            curveId = int32(se.curveId - 1); % 0-based
            dir = int8(se.dir);
            t0 = se.t0; t1 = se.t1;
        else
            curveId = int32(-1);
            dir = int8(0);
            t0 = 0; t1 = 0;
        end
        fwrite(fid, int32([meta.v1, meta.v2, meta.idx]), 'int32');
        fwrite(fid, curveId, 'int32');
        fwrite(fid, dir, 'int8');
        fwrite(fid, double([t0, t1]), 'double');
    end

    % Curves 存为 Curve<2,6>::dump 兼容格式
    for ci = 1:numCurves
        if ~isfield(curves{ci}, 'spline')
            error('曲线 %d 缺少 spline 字段，请先运行 interpolate_curves_spline', ci);
        end
        sp = curves{ci}.spline;
        writepp(fid, sp.xpp, sp.ypp);
    end

    % S 平滑性
    for vi = 1:numV
        entry = S{vi};
        if isempty(entry)
            fwrite(fid, uint32(0), 'uint32');
            continue;
        end
        if isstruct(entry)
            wf = entry.withFlag;
        else
            wf = entry;
        end
        pairCount = size(wf, 1);
        fwrite(fid, uint32(pairCount), 'uint32');
        if pairCount > 0
            fwrite(fid, int32(wf(:, [1, 3]) - 0), 'int32'); % edges 0-based already
            fwrite(fid, int8(wf(:, [2, 4])), 'int8');
        end
    end

    % newBoundaries
    for bi = 1:numNB
        nb = yinsetResult.newBoundaries{bi};
        fwrite(fid, uint8(nb.color(:)), 'uint8');
        edgeCount = numel(nb.edges);
        fwrite(fid, uint32(edgeCount), 'uint32');
        for k = 1:edgeCount
            fwrite(fid, int32(nb.edges(k).edgeId), 'int32');
            fwrite(fid, int8(nb.edges(k).dir), 'int8');
        end
    end
    % YinSet
    for yi = 1:numYin
        ys = yinsetResult.YinSet{yi};
        fwrite(fid, uint8(ys.color(:)), 'uint8');
        idxList = [];
        for cci = 1:numel(ys.components)
            idxList = [idxList, ys.components{cci}.indices]; %#ok<AGROW>
        end
        idxList = unique(idxList);
        fwrite(fid, uint32(numel(idxList)), 'uint32');
        fwrite(fid, uint32(idxList - 1), 'uint32'); % 0-based
    end
end
