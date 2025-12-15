function parent = build_containment_tree(containMat)
% BUILD_CONTAINMENT_TREE 根据包含矩阵构造父节点映射
% 输入:
%   containMat: NxN 逻辑矩阵，containMat(i,j)=true 表示 i 包含 j
% 输出:
%   parent: 1xN 向量，parent(j)=i 表示 i 是 j 的最近父级；无父则为 0

    n = size(containMat, 1);
    parent = zeros(1, n);
    for j = 1:n
        candidates = find(containMat(:, j))'; % 可能的父
        if isempty(candidates), continue; end
        % 选择最小包含者（不被其他候选包含）
        isMinimal = true(1, numel(candidates));
        for a = 1:numel(candidates)
            for b = 1:numel(candidates)
                if a == b, continue; end
                if containMat(candidates(b), candidates(a))
                    isMinimal(b) = false;
                    break;
                end
            end
        end
        mins = candidates(isMinimal);
        if ~isempty(mins)
            parent(j) = mins(1); % 若有多个等价最小父，取首个
        end
    end
end
