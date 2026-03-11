function [A, rhs] = cubic_spline_matrix(h, bc_type, y)
% CUBIC_SPLINE_MATRIX 构建三次样条插值的线性方程组 Ax = rhs
%
% 输入:
%   h       : 步长向量 (h_i = x_{i+1} - x_i). 
%             对于 n 个点，not-a-knot 需输入 n-1 个 h.
%             对于 periodic，需输入 n-1 个 h (最后一段 wrap around 由逻辑处理).
%   bc_type : 边界条件字符串 ('periodic' 或 'not-a-knot')
%   y       : 样本点函数值向量 (长度为 n)
%
% 输出:
%   A       : 系数矩阵
%   rhs     : 方程右端项 (Right Hand Side)
%
% 说明:
%   求解出的 x = A \ rhs 对应的是样本点处的二阶导数 M (Moments).

    % 1. 输入处理与预检查
    h = h(:); % 确保列向量
    y = y(:); 
    n = length(y);
    
    if length(h) ~= n - 1
        error('Input h must have length n-1 (where n is length of y).');
    end

    % 2. 初始化
    % -----------------------------------------------------------
    
    if strcmp(bc_type, 'periodic')
        % === 周期性边界 (Periodic) ===
        % 假设: y(1) == y(end). 
        % 未知数数量: n-1 (因为 M_n = M_1)
        % 方程组大小: (n-1) * (n-1)
        
        N = n - 1; 
        A = zeros(N, N);
        rhs = zeros(N, 1);
        
        % 周期性处理：我们需要 h 和 y 的环绕索引
        % h 向量通常是 [h1, h2, ..., h_{n-1}]
        % y 向量是 [y1, y2, ..., y_{n-1}, y_n] (其中 y_n = y1)
        
        for i = 1:N
            % 确定索引 (Wrap around logic)
            % i-1: 如果 i=1, 前一个点是 N
            idx_prev = mod(i-2, N) + 1; 
            % i+1: 如果 i=N, 后一个点是 1
            idx_next = mod(i, N) + 1;
            
            % 获取步长
            % h_i 对应当前点 i 到 i+1 的距离
            % h_{i-1} 对应前一个点 i-1 到 i 的距离
            hi = h(i);
            hi_prev = h(idx_prev);
            
            % 构造矩阵行: h_{i-1}*M_{i-1} + 2*(h_{i-1}+h_i)*M_i + h_i*M_{i+1}
            A(i, idx_prev) = hi_prev;
            A(i, i)        = 2 * (hi_prev + hi);
            A(i, idx_next) = hi;
            
            % 构造右端项: 6 * (差商后 - 差商前)
            % 注意 y 的索引：y(1)...y(N) 是独立的，y(N+1) 是 y(1)
            % 为了方便，我们直接取 y 向量，注意 y(idx_next) 当 idx_next=1 时应取 y(n) 或 y(1)
            % 但为了差商计算的一致性，最好直接用传入的 y
            
            val_curr = y(i);
            val_prev = y(idx_prev); 
            val_next = y(idx_next);
            
            % 特殊情况修正：当涉及到边界跨越时
            if i == 1
                val_prev = y(N); % 前一个点实际上是倒数第二个 y
            end
            if i == N
                val_next = y(1); % 后一个点是第一个 y
                % 或者直接使用 y(n) 因为 y(n)==y(1)
                val_next = y(n);
            end

            diff_first = (val_next - val_curr) / hi;
            diff_second = (val_curr - val_prev) / hi_prev;
            
            rhs(i) = 6 * (diff_first - diff_second);
        end

    elseif strcmp(bc_type, 'not-a-knot')
        % === 非扭结边界 (Not-a-Knot) ===
        % 未知数数量: n (M_1 ... M_n)
        % 方程组大小: n * n
        
        N = n;
        A = zeros(N, N);
        rhs = zeros(N, 1);
        
        % A. 内部节点方程 (i = 2 到 n-1)
        % 标准连续性方程
        for i = 2:N-1
            hi = h(i);       % x_{i+1} - x_i
            hi_prev = h(i-1); % x_i - x_{i-1}
            
            A(i, i-1) = hi_prev / (hi_prev + hi);
            A(i, i)   = 2;
            A(i, i+1) = hi / (hi_prev + hi);
            
            % 差商
            d1 = (y(i+1) - y(i)) / hi;
            d2 = (y(i) - y(i-1)) / hi_prev;
            rhs(i) = 6 * (d1 - d2) / (hi_prev + hi);
        end
        
        % B. 边界条件方程
        % Not-a-knot 条件意味着 S''' 在 x_2 和 x_{n-1} 处连续
        % 这等价于前两段是同一多项式，后两段是同一多项式
        
        % Row 1: h2*M1 - (h1+h2)*M2 + h1*M3 = 0
        h1 = h(1); 
        h2 = h(2);
        A(1, 1) = h2 / (h1 + h2);
        A(1, 2) = -1;
        A(1, 3) = h1 / (h1 + h2);
        rhs(1) = 0;
        
        % Row N: h_{n-1}*M_{n-2} - (h_{n-2}+h_{n-1})*M_{n-1} + h_{n-2}*M_n = 0
        hn_1 = h(end);   % h_{n-1}
        hn_2 = h(end-1); % h_{n-2}
        A(N, N-2) = hn_1 / (hn_2 + hn_1);
        A(N, N-1) = -1;
        A(N, N)   = hn_2 / (hn_2 + hn_1);
        rhs(N) = 0;
        
    else
        error('Unknown bc_type. Use "periodic" or "not-a-knot".');
    end
end