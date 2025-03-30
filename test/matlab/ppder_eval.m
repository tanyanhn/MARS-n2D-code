function d = ppder_eval(t, k, pp)
% PPDER_EVAL 计算pp样条的k阶导数在指定点的值
% 输入参数：
%   t: 求值点（标量或向量）
%   k: 导数阶数（非负整数）
%   pp: pp样条结构体
% 输出参数：
%   d: 导数在t处的值（与t同维度）

% 检查k是否为非负整数
if k < 0 || fix(k) ~= k
    error('k必须是非负整数');
end

% 处理k=0的情况（直接返回原函数值）
if k == 0
    d = ppval(pp, t);
    return;
end

% 提取pp结构信息
coefs = pp.coefs;
order = pp.order;

% 处理k超过可导次数的情况（返回全零）
if k >= order
    d = zeros(size(t));
    return;
end

% 对系数进行k次求导处理
for j = 1:k
    [~, m] = size(coefs);
    if m == 1  % 已为常数项，导数后全零
        coefs = zeros(size(coefs, 1), 0);
        break;
    end
    factors = (m-1:-1:1);  % 各阶次数因子
    coefs = coefs(:, 1:end-1) .* factors;  % 系数乘以次数并截断
end

% 构造导数后的pp结构
pp_der = struct(...
    'form', 'pp',...
    'breaks', pp.breaks,...
    'coefs', coefs,...
    'pieces', pp.pieces,...
    'order', order - k,...
    'dim', pp.dim);

% 计算导数在t处的值
d = ppval(pp_der, t);
end
