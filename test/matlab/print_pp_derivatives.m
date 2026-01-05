function print_pp_derivatives(pp, t0, max_order)
% PRINT_PP_DERIVATIVES 打印样条函数在 t0 处的各阶导数值
%
% 输入:
%   pp        - 样条结构体 (例如通过 spline, mkpp, pchip 等生成)
%   t0        - 需要计算导数的点
%   max_order - (可选) 最大导数阶数，默认为 4
%
% 依赖:
%   此函数依赖于你工作区中已有的 ppder_eval 函数

% use case
% print_pp_derivatives(xsp1, t0);

    % 如果未提供 max_order，则默认为 4
    if nargin < 3
        max_order = 4;
    end

    % fprintf('--- 在 t0 = %g 处的导数分析 ---\n', t0);
    
    for k = 0:max_order
        % 调用你的自定义函数 ppder_eval
        val = ppder_eval(t0, k, pp);
        
        % 打印结果 (使用 num2str 保持和你源代码一致的格式，或者用 %g)
        disp([num2str(k), '阶导数最大: ', num2str(max(abs(val)))]);
    end
    
    fprintf('---------------------------\n');
end