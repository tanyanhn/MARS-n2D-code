%% MARS plot
close all
% clear
% hold on
figure
lineColor = '';
fillColor = getColor();
dir = "../../results/TrackInterface/";
% Shape = "Disk";
Shape = "Graph";
part = 41;
% vel = "Vortex";
vel = "Deformation";
T = 4;
Order = 4;
filedir = dir + Shape + num2str(part) + vel + "T" + num2str(T) + "Order" + num2str(Order) + "/";
grid = 32;
plotT = 1 / 4;
step = T * grid * 8 * plotT;
type = "_c.dat";
filename = "4Circle_grid" + num2str(grid) + "_Step" + num2str(step) + type;
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
    for k = (part + 1):-1:1

        if ~tensor
            sf = readYinSet(hd);
            plotYinSet(sf, lineColor, fillColor{k}); hold on
        else
            vecSf = readVecYinSet(hd);
            if k == part + 1
                continue;
            end
            for i=1:1:length(vecSf)
                plotYinSet(vecSf{i}, lineColor, fillColor{k}); hold on
            end
        end
    end
end

axis equal