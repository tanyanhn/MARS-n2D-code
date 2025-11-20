close all
clear

% Shape = "Disk4";
% Shape = "Disk5";
Shape = "Graph41";
vel = "Vortex";
% vel = "Deformation";
T = 16;
Order = 4;
dir = "../../results/TrackInterface/";
filedir = sprintf('%s%s%sT%dOrder%d/', dir, Shape, vel, T, Order);
% filedir = "../../results/TrackInterface/Disk5DeformationT4Order4/";
% filedir = "../../results/TrackInterface/Graph41DeformationT4Order4/";
% filedir = "../../results/TrackInterface/Graph41VortexT16Order4/";
filename1 = "1Order4_grid32_00marksHistory.dat";
filename2 = "1Order4_grid32_00LengthHistory.dat";
% ylim_range = [0 9];

hd1 = fopen(filedir + filename1);
hd2 = fopen(filedir + filename2);

rows1 = fread(hd1, 1, 'int');
cols1 = fread(hd1, 1, 'int');
rows2 = fread(hd2, 1, 'int');
cols2 = fread(hd2, 1, 'int');

A = fread(hd1, [rows1, cols1], "int");
B = fread(hd2, [rows2, cols2], "double");

A = A ./ A(:, 1);
B = B ./ B(:, 1);

% fillColor = 'kmgrcb';
fillColor = getColor();
fillColor{rows1} = 'k';
t = linspace(0, 1, cols1);
fig = figure;
ax1 = axes('FontSize', 23); % 创建坐标轴时直接指定属性
hold(ax1, 'on');
for i = 1:rows1
    xx = 0:0.00001:1; 
    yy = spline(t, A(i, :), xx);
    h = plot(ax1, xx, yy, 'Color', fillColor{i}, 'LineWidth', 5); 
    set(h, 'LineJoin', 'round');
    set(gcf, 'Renderer', 'opengl');
    set(h, 'LineStyle', '-'); % 确保是实线 
    hold on;
end
mkdir('tmp');
filename = 'tmp/markers.png';
% ylim(ax1, ylim_range);
set(gcf, 'Color', 'w');
export_fig('-q101', '-m4', fig, filename);
hold(ax1, 'off');
fig = figure;
ax2 = axes('FontSize', 23); % 创建坐标轴时直接指定属性
hold(ax2, 'on');
for i = 1:rows2
    xx = 0:0.00001:1; 
    yy = spline(t, B(i, :), xx);
    plot(ax2, xx, yy, 'Color', fillColor{i}, 'LineWidth', 5); 
    set(h, 'LineJoin', 'round');
    set(gcf, 'Renderer', 'opengl');
    set(h, 'LineStyle', '-'); % 确保是实线 
    hold on;
end
filename = 'tmp/lengths.png';
% ylim(ax2, ylim_range);
set(gcf, 'Color', 'w');
export_fig('-q101', '-m4', fig, filename);
hold(ax2, 'off');

fclose(hd1);
fclose(hd2);