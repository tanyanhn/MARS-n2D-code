%% MARS plot
close all
fclose('all');
clear
addpath 'export_fig-3.54';
% hold on
lineColor = '';
% dir = "results/";
dir = "../../results/TrackInterface/";
Shape = "Disk";
% Shape = "Graph";
% Shape = "Raccoon";
% Shape = "Pig";
% fillColor = getColor(Shape);
fillColor = {'r', 'g'};
part = 2;
% vel = 'Vortex';
vel = 'Deformation';
dynamic = "";
% dynamic = "Dynamic";
T = 2;
Order = 4;
filedir = dir + Shape + num2str(part) + vel + "T" + num2str(T) + "Order" + num2str(Order)+ dynamic + "/";
grid = 32;
steps = [0, 2, 4];
stepSize = 128;
type = "_c.dat"; plotF = false;
% type = "_f.dat"; plotF = true;
% figHandles = [];
plotInaSinglegraph = false;
if plotInaSinglegraph
    f = figure('Position', [100, 200, 1200, 350], 'Color', 'w');
end
for iStep = 1:length(steps)
    plotI = steps(iStep);
%     plotI = 1;
    step = plotI * stepSize;
    % if plotI == steps
    %     step = T * 256;
    % end
    filename = "1Order4_grid" + num2str(grid) + "_Step" + num2str(step) + type;
    if plotInaSinglegraph
        subplot(1, 3, iStep);
    else
        f = figure;
    end
    if part == 41
        part = 5;
    end
    hd = fopen(filedir + filename);
    tensor = true;
    volume = false;
    if volume
        dat = readCellVolume(hd);
    else
        for k = 1:1:(part + 1)

            if ~tensor
                sf = readYinSet(hd);
                plotYinSet(sf, lineColor, fillColor{k}); hold on
                clear sf;
            else
                vecSf = readVecYinSet(hd);
                if k == part + 1
                    continue;
                end
                if plotF 
                    for i=length(vecSf)-1:-1:1
                        plotYinSet(vecSf{i}, lineColor, fillColor{i}); hold on
                    end
                else 
                    for i=1:1:length(vecSf)
                        plotYinSet(vecSf{i}, lineColor, fillColor{k}); hold on
                    end
                end
                clear vecSf;
            end
        end
    end
    if plotInaSinglegraph
        hold on;
        box on;
        axis([0 1 0 1]);
        axis square; % 强制坐标轴为正方形
        set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1, 'LineWidth', 1);
    else
        saveFiguresAsFrames(vel, [f], plotI);
        close(f);
    end
    fclose(hd);
end
if plotInaSinglegraph
    saveFiguresAsFrames(vel, [f], plotI, 'Tight', 0);
end
% saveFiguresAsFrames(vel, figHandles, 0, 'Tight', 1);

