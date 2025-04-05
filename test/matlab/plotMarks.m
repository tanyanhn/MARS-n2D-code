close all

% filedir = "../../results/TrackInterface/Disk5DeformationT4Order4/";
filedir = "../../results/TrackInterface/Graph41DeformationT4Order4/";
filename1 = "00marksHistory.dat";
filename2 = "00LengthHistory.dat";

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
figure
ax1 = axes('FontSize', 14); % 创建坐标轴时直接指定属性
hold(ax1, 'on');
for i = 1:rows1
    plot(ax1, t, A(i, :), 'Color', fillColor{i}, 'LineWidth', 5); hold on;
end
hold(ax1, 'off');
figure
ax2 = axes('FontSize', 14); % 创建坐标轴时直接指定属性
hold(ax2, 'on');
for i = 1:rows2
    plot(ax2, t, B(i, :), 'Color', fillColor{i}, 'LineWidth', 5); hold on;
end
hold(ax2, 'off');

fclose(hd1);
fclose(hd2);