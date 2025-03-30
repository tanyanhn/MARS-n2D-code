
function m = solve_periodic_fft(n, b)
A = zeros(n, n);


% 构建循环卷积核c
c = zeros(n, 1);
c(1) = 66;        % 中心点 (位移0)
c(2) = 26;        % 右侧第一个点 (位移+1)
c(3) = 1;         % 右侧第二个点 (位移+2)
c(end) = 26;      % 左侧第一个点 (位移-1)
c(end-1) = 1;     % 左侧第二个点 (位移-2)

% % FFT快速求解
% lambda = fft(c);           % 特征值谱
% B = fft(b(:));             % 右侧项FFT
% m = real(ifft(B ./ lambda)); % 频域相除后逆变换

A = gallery('circul', c); % 生成周期五对角矩阵
m = A \ b; % 反斜杠运算符自动调用稀疏求解器


% 验证解的正确性 (可选)
if max(abs(imag(m))) > 1e-10
    warning('解存在显著虚部，请检查输入条件');
end
end