function [A, condA] = not_a_knot_interpolation(n, r_tiny)
    % n: 插值节点数-1
    % r_tiny: 小间距比 (l_i = r_tiny * h)
    
    % Step 1: 设置 l_0, l_1, ..., l_{n-1}
    % 构造一个均匀间隔的插值节点
    h = 1;  % 基础间距
    l = zeros(1, n);
    
    % 手动设置节点，保证 l_2 - l_1 = r_tiny * h
    rh = r_tiny * h;
    rrh = r_tiny * rh;
    % l(1) = h;
    % l(2) = r_tiny * h;
    
    for i = 1:1:n
%         l(i) = rh;  % 其它间距为 h
        l(i) = (10 ^ (rand() * log10(rh)));
%         l(i + 1) = h;
        % l(i + 2) = rh;
        % l(i + 1) = rand() * (h - rh) + rh;
    end
    l = l(1:n);

%     l(1) = h;
%     l(2) = rh;
%     l(end - 1) = h;
%     l(end) = rh;

    % Step 2: 构造 not-a-knot 插值矩阵 A
    M = n + 1;
%     M = n;
    A = zeros(M, M);
    
        for i = 2:M-1
%         for i = 1:M
            % 对于开放曲线:
            % h 索引: h(1)是P1->P2, h(i-1)是Pi-1->Pi, h(i)是Pi->Pi+1
            prei = mod(i - 2 + M, M) + 1;  
            posi = mod(i, M) + 1;
            h_prev_val = l(prei);
            h_curr_val = l(i);
            
            sum_h = h_prev_val + h_curr_val;
            
            mu     = h_prev_val / sum_h;
            lambda = h_curr_val / sum_h;
            
            A(i, prei) = mu;
            A(i, i)   = 2;
            A(i, posi) = lambda;
        end

        % --- 第一个点 (Row 1) ---
        h1 = l(1);
        h2 = l(2);
        A(1, 1) = h2 / (h1 + h2);
        A(1, 2) = -1;
        A(1, 3) = h1/ (h1 + h2);
        
        % --- 最后一个点 (Row M) ---
        % 方程对称: h_{end} * M_{end-2} - (...) * M_{end-1} + h_{end-1} * M_{end} = 0
        h_end_1 = l(end-1); % h_{M-2}
        h_end   = l(end);   % h_{M-1}
        
        A(M, M-2) = h_end/ (h_end_1 + h_end);
        A(M, M-1) = -1;
        A(M, M)   = h_end_1 / (h_end_1 + h_end);
    
    % Step 3: 计算矩阵 A 的条件数
    condA = cond(A, 'inf');
%     condA = cond(A, 2);

    % if condA > 10
    %     fprintf('矩阵 A 的条件数为: %.4f\n', condA);
    % end
    
    % 显示结果
    % fprintf('矩阵 A 的条件数为: %.4f\n', cond_A);
end
