function pp = sympoly2pp(polyCell, breaks, order)
% sympoly2pp 将符号多项式转换为pp格式的分段多项式（加速版）
% 输入：
%   polyCell - 符号多项式单元数组，每个元素对应一个区间的多项式
%   breaks   - 区间断点，格式为[a0, a1, ..., an]
%   order    - 样条的阶数（例如，4表示三次样条）
% 输出：
%   pp - pp结构，包含分段多项式信息

n = numel(polyCell);  % 分段数
coefs = zeros(n, order);  % 初始化系数矩阵

for i = 1:n
    p = polyCell{i};
    coeff = sym2poly(p);  % 获取多项式系数（按降幂排列）
    numCoeff = numel(coeff);
    
    % 左侧补零至指定阶数
    if numCoeff < order
        coeff = [zeros(1, order - numCoeff), coeff];
    end
    
    coefs(i, :) = coeff(1:order);  % 确保长度匹配
end

% 构建pp结构
pp = struct(...
    'form', 'pp',...
    'breaks', breaks,...
    'coefs', coefs,...
    'pieces', n,...
    'order', order,...
    'dim', 1);
end
