function write_geometry_bin(filepath, V, E, S, curves, splineE, F)
% WRITE_GEOMETRY_BIN 将 (V,E,S) 及样条数据写出为紧凑二进制
% 格式 (小端):
%   magic[4] = 'GEOB'
%   uint32 version = 2
%   uint32 numV, numE, numCurves, numRegions
%   V: double[2*numV] (先 x 行后 y 行按列写)
%   Edge 元数据 (numE 次):
%       int32 v1, v2, idx
%       int32 curveId   (对应 curves/spline，-1 表示缺失)
%       int8  dir       (+1 v1->v2, -1 反向)
%       double t0, t1   (在 curves{curveId}.spline 参数域的区间)
%   Curves 样条 (numCurves 次):
%       uint32 order            (应为 6)
%       uint32 nSeg             (样条段数)
%       double breaks[nSeg+1]
%       double coefsX[nSeg*order] (行主序，段优先)
%       double coefsY[nSeg*order]
%   S 光滑性 (numV 次):
%       uint32 pairCount
%       重复 pairCount 次: int32 eA, int8 fA, int32 eB, int8 fB
%   F 区域 (numRegions 次):
%       uint32 numBoundaries
%       重复 numBoundaries 次:
%           uint32 edgeCount
%           重复 edgeCount 次: int32 edgeId, int8 dir (+1 v1->v2, -1 反向)
%
% 读取时可重建:
%   V: 顶点坐标
%   E: 每条边记录 (v1,v2,idx, curveRef=curveId, dir, [t0,t1])
%   S: 顶点平滑对
%   F: 区域由边序列给出
%
% 依赖: splineE 需由 interpolate_curves_spline 生成 (避免冗余存 xpp/ypp).

    fid = fopen(filepath, 'w');
    if fid < 0
        error('无法打开文件: %s', filepath);
    end
    c = onCleanup(@() fclose(fid));

    numV = size(V, 2);
    numE = numel(E);
    numCurves = numel(curves);
    if nargin < 7
        F = [];
    end
    numRegions = numel(F);

    % header
    fwrite(fid, 'GEOB', 'uint8');
    fwrite(fid, uint32(2), 'uint32');
    fwrite(fid, uint32([numV, numE, numCurves, numRegions]), 'uint32');

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

    % Curves spline
    for ci = 1:numCurves
        if ~isfield(curves{ci}, 'spline')
            error('曲线 %d 缺少 spline 字段，请先运行 interpolate_curves_spline', ci);
        end
        sp = curves{ci}.spline;
        xpp = sp.xpp; ypp = sp.ypp;
        order = uint32(size(xpp.coefs, 2));
        nSeg = uint32(size(xpp.coefs, 1));
        fwrite(fid, order, 'uint32');
        fwrite(fid, nSeg, 'uint32');
        fwrite(fid, xpp.breaks, 'double'); % 长度 nSeg+1
        fwrite(fid, xpp.coefs.', 'double'); % 转置使行主序写入
        fwrite(fid, ypp.coefs.', 'double');
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

    % F 区域信息（仅存 edgeId 和 dir）
    for fi = 1:numRegions
        bnds = F(fi).boundaries;
        numB = numel(bnds);
        fwrite(fid, uint32(numB), 'uint32');
        for bi = 1:numB
            loop = bnds{bi}.edgeLoop;
            edgeCount = numel(loop);
            fwrite(fid, uint32(edgeCount), 'uint32');
            for ej = 1:edgeCount
                fwrite(fid, int32(loop(ej).edgeId), 'int32'); % 0-based
                fwrite(fid, int8(loop(ej).dir), 'int8');
            end
        end
    end
end
