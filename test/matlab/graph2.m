close all 

figure 
max_step = 1e-2;
move = [-1.7; -5.2];
scale = [1 / 14];

controlpoints1 = [9, 2.75, 2.75, 8.5;...
                 15.95, 15.85, 9.5, 9]'; 
controlpoints2 = [10, 8.5, 12, 8.5;...
                 13., 10.5, 9, 9]' ;
controlpoints2 = reverse(controlpoints2);
controlpoints3 = [9, 13.75, 10.5, 10;...
                 15.95, 15.75, 14, 13]';
controlpoints3 = reverse(controlpoints3);

controlpoints4 = [13.5, 10.523, 11.5, 12.25;...
                 14.75, 16, 14, 13]' ;
controlpoints5 = [12.25, 13., 12.75, 11.75;...
                 13., 12., 11.75, 11.25]' ;
controlpoints6 = [11.75, 8.826, 10.25, 12.75;...
                 11.25, 9.5, 8.25, 10.25]' ;
controlpoints7 = [13.5, 13.75, 14.25, 12.75;...
                 14.75, 12.5, 11.75, 10.25]' ;
controlpoints7 = reverse(controlpoints7);




edge1 = bezierCurve(controlpoints1, max_step)';
edge2 = bezierCurve(controlpoints2, max_step)';
edge3 = bezierCurve(controlpoints3, max_step)';

edge4 = bezierCurve(controlpoints4, max_step)';
edge5 = bezierCurve(controlpoints5, max_step)';
edge6 = bezierCurve(controlpoints6, max_step)';
edge7 = bezierCurve(controlpoints7, max_step)';

k1 = 272; k2 = 315;
k3 = 301; k4 = 326;
[edge2, edge6] = attachPoint(edge2, edge6, k1, k4);
[edge3, edge4] = attachPoint(edge3, edge4, k2, k3);
k1 = k1 - 1 + size(edge1, 2) - 1;
k2 = k2 - 1 + size(edge1, 2) - 1 + size(edge2, 2) - 1;
k3 = k3 - 1;
k4 = k4 - 1 + size(edge4, 2) - 1 + size(edge5, 2) - 1;


edge1 = moveScale(edge1, move, scale);
edge2 = moveScale(edge2, move, scale);
edge3 = moveScale(edge3, move, scale);
edge4 = moveScale(edge4, move, scale);
edge5 = moveScale(edge5, move, scale);
edge6 = moveScale(edge6, move, scale);
edge7 = moveScale(edge7, move, scale);

jordanCurve1 = [edge1(:, 1:end-1), edge2(:, 1:end-1), edge3(:, 1:end)];
jordanCurve2 = [edge4(:, 1:end-1), edge5(:, 1:end-1), edge6(:, 1:end-1), edge7(:, 1:end)];

jordanCurve1(:, 1) - jordanCurve1(:, end)
jordanCurve2(:, 1) - jordanCurve2(:, end)



plot(edge1(1,:), edge1(2,:)); hold on
plot(edge2(1,:), edge2(2,:)); hold on
plot(edge3(1,:), edge3(2,:)); hold on

plot(edge4(1,:), edge4(2,:)); hold on
plot(edge5(1,:), edge5(2,:)); hold on
plot(edge6(1,:), edge6(2,:)); hold on
plot(edge7(1,:), edge7(2,:)); hold on

axis equal


filename = "../data1/Graph4_1.input";
hd = fopen(filename, 'w');

dumpJordanCurve(jordanCurve1, hd);
dumpJordanCurve(jordanCurve2, hd);
fwrite(hd, [k1, k2, k3, k4], 'int');
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