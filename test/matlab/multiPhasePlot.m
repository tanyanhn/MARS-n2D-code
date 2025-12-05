%% MARS plot
close all
fclose('all');
clear
addpath 'export_fig-3.54';
% hold on
lineColor = '';
fillColor = getColor();
% dir = "results/";
dir = "../../results/TrackInterface/";
Shape = "Disk";
% Shape = "Graph";
part = 4;
% part = 41;
vel = 'Vortex';
% vel = 'Deformation';
dynamic = "";
% dynamic = "Dynamic";
T = 4;
Order = 4;
filedir = dir + Shape + num2str(part) + vel + "T" + num2str(T) + "Order" + num2str(Order)+ dynamic + "/";
grid = 32;
steps = 4;
stepSize = 256;
type = "_c.dat"; plotF = false;
% type = "_f.dat"; plotF = true;
figHandles = [];
for plotI = 0:1:steps
%     plotI = 1;
    step = plotI * stepSize;
    if plotI == steps
        step = T * 256;
    end
    filename = "1Order4_grid" + num2str(grid) + "_Step" + num2str(step) + type;
    f = figure;
    if part == 41
        part = 5;
    end
    % filedir = "results/CutCell/";
    % filedir = "results/InterfaceGraph/";
    % filedir = "../../results/TrackInterface/Disk4VortexT8Order4/";
    % filedir = "../../results/TrackInterface/Disk2DeformationT2Order4/";
    % filedir = "../../results/TrackInterface/Graph41VortexT8Order4/";
    % filedir = "../../test/data1/";
    % round = "No4_";

    % filename = "spadjor-13.input.dat";
    % filename = "localVolumes.dat";
    % filename = "localYinset.dat";
    % filename = "4Circle_grid32_Step256_c.dat";
    % filename = "4Circle_grid32.dat";
    % filename = "00test.dat";
    % filename = "test.dat";
    % Order = 4;
    % N = 4;
    % Shape = "Rose";
    % Shape = "Rectangle";
    % filename = num2str(Order) + Shape + num2str(N) + "_" + num2str(k) + ".dat";
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
    saveFiguresAsFrames(vel, [f], plotI);
    close(f);
    fclose(hd);
end
% saveFiguresAsFrames(vel, figHandles, 0, 'Tight', 1);

