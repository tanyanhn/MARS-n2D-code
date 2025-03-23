% 绘制多项式曲线
% 输入多项式系数和参数范围，绘制曲线

clc;
close all;

% 输入多项式系数
% coeff = input('请输入多项式系数（按降幂排列，例如[3, 0, 2]对应3x² + 2）: ');
% x_start = input('请输入x的起始值：');
% x_end = input('请输入x的结束值：');
coeff = [0.000870393697073973209542 -3.27447776442540955779E-7 3.56190587759540292362E-13 -1.11022302462515654042E-16];
xrange = [0 2.7253688180545188E-7];
orig_x1 = xrange(1);
orig_x2 = xrange(2);
exp_x1 = orig_x1 - 1e-7;
exp_x2 = orig_x2 + 1e-7;

% 生成扩展范围数据
x_exp = linspace(exp_x1, exp_x2, 1000);
y_exp = polyval(coeff, x_exp);

% 创建图形窗口
figure;

% 先绘制背景区域
y_range = [min(y_exp), max(y_exp)];  % 自动适应y轴范围
if diff(y_range) == 0  % 处理常函数情况
    y_range = y_range + [-1, 1];
end
rectangle('Position', [orig_x1 y_range(1) orig_x2-orig_x1 y_range(2)-y_range(1)],...
          'FaceColor', [0.9 0.8 0.8], 'EdgeColor','none');
hold on;

% 绘制扩展曲线
plot(x_exp, y_exp, 'LineWidth', 2, 'Color', [0 0.447 0.741]);
xlabel('x');
ylabel('y');
title('多项式曲线（蓝色为扩展范围，粉色为原始区间）');
grid on;

% 绘制原始区间边界线
plot([orig_x1 orig_x1], ylim, '--', 'Color', [1 0 0], 'LineWidth', 1);
plot([orig_x2 orig_x2], ylim, '--', 'Color', [1 0 0], 'LineWidth', 1);

% 生成多项式表达式字符串
terms = {};
for i = 1:length(coeff)
    power = length(coeff) - i;
    coeff_value = coeff(i);
    
    if coeff_value ~= 0
        % 处理不同次数项
        if power == 0
            term = format_coeff(coeff_value, true);
        else
            term = format_coeff(coeff_value, false);
            if power > 1
                term = [term 'x^{' num2str(power) '}'];
            elseif power == 1
                term = [term 'x'];
            end
        end
        
        % 处理符号显示
        if i > 1
            if coeff_value > 0
                terms{end+1} = ['+ ' term];
            else
                terms{end+1} = ['- ' term(2:end)];
            end
        else
            terms{end+1} = term;
        end
    end
end

% 生成最终多项式字符串
if isempty(terms)
    poly_str = '0';
else
    poly_str = strjoin(terms, ' ');
end

% 添加复合图例
legend_items = {
    ['y = ' poly_str], 
    '原始参数范围', 
    '区间边界'
};
legend(legend_items, 'Location', 'best');

%% 辅助函数：格式化系数显示
function str = format_coeff(value, is_constant)
    if value == 0
        str = '';
        return;
    end
    
    if is_constant
        if value == floor(value)
            str = sprintf('%d', value);
        else
            str = sprintf('%.2f', value);
        end
    else
        if abs(value) == 1
            str = '';
            if value == -1
                str = '-';
            end
        else
            if value == floor(value)
                str = sprintf('%d', value);
            else
                str = sprintf('%.2f', value);
            end
        end
    end
end
