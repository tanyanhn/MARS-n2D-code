close all 
clear

figure 
max_step = 1e-2;
move = [-1.7; -5.26];
scale = [1 / 14];

controlpoints1 =...
   [9.,       15.95;
    2.,       16.15;
    2.,       9.6;
    8,        9.];
controlpoints2 =...
   [8,        9.;
    12,       8.6;
    8.,      11.;
    10,       13.] ;
controlpoints3 =... 
   [10,       13;
    11,       14;
    12.5,     15.85;
    9.,       15.95];

adjust = 0.196445;
a = 11. -adjust;
b= 13 +adjust;
controlpoints4 =...
   [13.5,     14.75;
    10.5,     16;
    a,        14;
    12.,      13];
adjust = -0.423635;
c = 12.8 -adjust;
d= 8.8 +adjust;
controlpoints5 = ...
   [12.,      13;
    b,        12;
    c,        12.75;
    10.8,     10.75];
controlpoints6 =...
   [10.8,    10.75;
    d,        8.75;
    11.25,    8.75;
    12.75,    10.25];
controlpoints7 =...
   [12.75,   10.25;
    14.25,    11.75;
    13.75,    12.5;
    13.5,     14.75];


edge1 = bezierCurve(controlpoints1, max_step)';
edge2 = bezierCurve(controlpoints2, max_step)';
edge3 = bezierCurve(controlpoints3, max_step)';

edge4 = bezierCurve(controlpoints4, max_step)';
edge5 = bezierCurve(controlpoints5, max_step)';
edge6 = bezierCurve(controlpoints6, max_step)';
edge7 = bezierCurve(controlpoints7, max_step)';



edge1 = moveScale(edge1, move, scale);
edge2 = moveScale(edge2, move, scale);
edge3 = moveScale(edge3, move, scale);
edge4 = moveScale(edge4, move, scale);
edge5 = moveScale(edge5, move, scale);
edge6 = moveScale(edge6, move, scale);
edge7 = moveScale(edge7, move, scale);
[~, k1] = min(abs(edge2(2,:) - 0.305636));
[~, k4] = min(abs(edge6(2,:) - 0.305507980538836));
[~, k2] = min(abs(edge3(2,:) - 0.672228));
[~, k3] = min(abs(edge4(2,:) - 0.672226));
[edge2, edge6] = attachPoint(edge2, edge6, k1, k4);
[edge3, edge4] = attachPoint(edge3, edge4, k2, k3);
plot(edge2(1, k1), edge2(2, k1), '*'); hold on
plot(edge2(1, k1 + 1), edge2(2, k1 + 1), '*'); hold on
plot(edge2(1, k1 - 1), edge2(2, k1 - 1), '*'); hold on
plot(edge4(1, k3), edge4(2, k3), '*'); hold on
plot(edge4(1, k3 + 1), edge4(2, k3 + 1), '*'); hold on
plot(edge4(1, k3 - 1), edge4(2, k3 - 1), '*'); hold on


% p = [7.25  9; 
%     13.25 14.5];
% p = [7.25 8.25;
%     11 11.5];
p = [11.1 ;
    11.85];
p = moveScale(p, move, scale);

k1 = k1 - 1 + size(edge1, 2) - 1;
k2 = k2 - 1 + size(edge1, 2) - 1 + size(edge2, 2) - 1;
k3 = k3 - 1;
k4 = k4 - 1 + size(edge4, 2) - 1 + size(edge5, 2) - 1;
jordanCurve1 = [edge1(:, 1:end-1), edge2(:, 1:end-1), edge3(:, 1:end)];
jordanCurve2 = [edge4(:, 1:end-1), edge5(:, 1:end-1), edge6(:, 1:end-1), edge7(:, 1:end)];

jordanCurve1(:, 1) - jordanCurve1(:, end);
jordanCurve2(:, 1) - jordanCurve2(:, end);

[xsp2, ysp2, t2] = construct_spline(jordanCurve2, 6);
t_total2 = xsp2.breaks(end);
t_eval2 = linspace(0, t_total2, 10000);
t_eval2 = sort([t_eval2, t2(k3 + 1), t2(k4 + 1)]);
x_eval2 = ppval(xsp2, t_eval2);
y_eval2 = ppval(ysp2, t_eval2);
plot(x_eval2(1, :), y_eval2(1,:)); hold on;

x3 = ppval(xsp2, t2(k3 + 1));
y3 = ppval(ysp2, t2(k3 + 1));
x4 = ppval(xsp2, t2(k4 + 1));
y4 = ppval(ysp2, t2(k4 + 1));
check(3) = norm(jordanCurve2(:, k3 + 1) - [x3; y3]);
check(4) = norm(jordanCurve2(:, k4 + 1) - [x4; y4]);
plot(x3, y3, 'o'); hold on;

[xsp1, ysp1, t1] = quinticSpline(jordanCurve1);
t_total1 = xsp1.breaks(end);
t_eval1 = linspace(0, t_total1, 10000);
t_eval1 = sort([t_eval1, t1(k1 + 1), t1(k2 + 1)]);
x_eval1 = ppval(xsp1, t_eval1);
y_eval1 = ppval(ysp1, t_eval1);
plot(x_eval1(1, :), y_eval1(1,:)); hold on;

x1 = ppval(xsp1, t1(k1 + 1));
y1 = ppval(ysp1, t1(k1 + 1));
x2 = ppval(xsp1, t1(k2 + 1));
y2 = ppval(ysp1, t1(k2 + 1));
check(1) = norm(jordanCurve1(:, k1 + 1) - [x1; y1]);
check(2) = norm(jordanCurve1(:, k2 + 1) - [x2; y2]);
plot(x1, y1, 'o'); hold on;
% t0 = t1(1);
% tend = t1(end);
% pp = xsp1;
% disp(['0阶导数: ', num2str(ppder_eval(t0, 0, pp))])  % 0.125
% disp(['1阶导数: ', num2str(ppder_eval(t0, 1, pp))])  % 0.75
% disp(['2阶导数: ', num2str(ppder_eval(t0, 2, pp))])  % 3
% disp(['3阶导数: ', num2str(ppder_eval(t0, 3, pp))])  % 6
% disp(['4阶导数: ', num2str(ppder_eval(t0, 4, pp))])  % 0

% plot(edge1(1,:), edge1(2,:)); hold on
% plot(edge2(1,:), edge2(2,:)); hold on
% plot(edge3(1,:), edge3(2,:)); hold on
% 
% plot(edge4(1,:), edge4(2,:)); hold on
% plot(edge5(1,:), edge5(2,:)); hold on
% plot(edge6(1,:), edge6(2,:)); hold on
% plot(edge7(1,:), edge7(2,:)); hold on

axis equal


filename = "../data1/Graph41.input";
hd = fopen(filename, 'w');
writepp(hd, xsp1, xsp2);
writepp(hd, xsp2, ysp2);

% dumpJordanCurve(jordanCurve1, hd);
% dumpJordanCurve(jordanCurve2, hd);
fwrite(hd, [t1(k1 + 1), t1(k2 + 1), t2(k3 + 1), t2(k4 + 1)], 'double');
fclose(hd);


% hd = fopen(filename, 'r');
% a = fread(hd, 1, "int");
% b = fread(hd, 2, "double");
% fclose(hd);


function [pts] = moveScale(pts, move, scale) 
pts = (pts + move) * scale;
end

function dumpJordanCurve(jordanCurve, hd)
num = size(jordanCurve, 2);
fwrite(hd, num, 'int');
fwrite(hd, jordanCurve, "double");
end

function [pts] = reverse(pts)
pts = pts(end:-1:1, :);
end

function [pts1, pts2] = attachPoint(pts1, pts2, i1, i2) 
mid = (pts1(:, i1) + pts2(:, i2)) / 2;
pts1(:, i1) = mid;
pts2(:, i2) = mid;
% p1 = pts1(:, i1);
% p2 = pts1(:, i1 + 1);
% p3 = pts2(:, i2);
% p4 = pts2(:, i2 + 1);
% mid = (p1 + p2 + p3 + p4) / 4;
% pts1 = [pts1(:, 1:i1); mid; pts1(:, i1 + 1:end)];
% pts2 = [pts2(:, 1:i2); mid; pts2(:, i2 + 1:end)];
end