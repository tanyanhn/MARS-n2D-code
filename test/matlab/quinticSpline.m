function [xsp, ysp, breaks] = quinticSpline(pts)

syms t;
oneMinusT = 1 - t;
onePlus3T = 1 + 3 * t;
onePlus3TPlus6T = 6 * t^2 + 3 * t + 1;

% 计算多项式乘积的系数向量
oneMinusT3 = oneMinusT * oneMinusT * oneMinusT;
A = oneMinusT3 * onePlus3TPlus6T;
B = subs(A, t, (1 - t));
C = t * oneMinusT3 * onePlus3T;
D = -subs(C, t, (1 - t));
E = 1 / 2 * t^2 * oneMinusT3;
F = subs(E, t, (1 - t));


n = size(pts, 2) - 1;
breaks = linspace(0, 1, n + 1);
A = subs(A, t, n * t);
B = subs(B, t, n * t);
C = subs(C, t, n * t);
D = subs(D, t, n * t);
E = subs(E, t, n * t);
F = subs(F, t, n * t);

xsp = fitCoordinate(pts(1, :), breaks, A, B, C, D, E, F);
ysp = fitCoordinate(pts(2, :), breaks, A, B, C, D, E, F);

end


function sp = fitCoordinate(f, breaks, A, B, C, D, E, F)
n = size(f, 2) - 1;
h = 1 / n;
b = zeros(n, 1);
syms t;


for i = 1:n
    % 处理周期性索引 (i+2, i+1, i-1, i-2)
    indices = mod(i - 1 + [2, 1, 0, -1, -2] + n, n) + 1;
    term = f(indices(1)) + 10*f(indices(2))...
          -10*f(indices(4)) - f(indices(5));
    b(i) = 5 * n * term;
end
lambda = solve_periodic_fft(n, b);
for i = 1:n
    indices = mod(i - 1 + [2, 1, 0, -1, -2] + n, n) + 1;
    term = f(indices(1)) + 2*f(indices(2)) - 6*f(indices(3))...
          + 2*f(indices(4)) + f(indices(5));
    b(i) = 20 * n * n * term;
end
gamma = solve_periodic_fft(n, b);

P = cell(n, 1);
for i = 1:n
    indices = mod(i - 1 + [2, 1, 0, -1, -2] + n, n) + 1;
    P{indices(4)} = expand(f(indices(4)) * A + f(indices(3)) * B + ...
        lambda(indices(4)) * C / n + lambda(indices(3)) * D / n + ...
        gamma(indices(4)) * E / (n * n) + gamma(indices(3)) * F / (n * n));
end


% for i = 1:n
%     indices = mod(i - 1 + [2, 1, 0, -1, -2] + n, n) + 1;
%     dP0 = P{indices(3)};
%     dP1 = P{indices(2)};
%     for k = 0:2
%         l = subs(dP0, t, h) - subs(dP1, t, 0);
%         if l > 1e-14
%             throw("fit not corect in " + num2str(k) + "k_" + num2str(i));
%         end
%         dP0 = diff(dP0, t);
%         dP1 = diff(dP1, t);
%     end
% end

sp = sympoly2pp(P, breaks, 6);
end