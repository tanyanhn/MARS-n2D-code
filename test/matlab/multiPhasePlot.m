%% MARS plot
close all
% clear
% hold on
figure
lineColor = '';
fillColor = getColor();
% filedir = "results/CutCell/";
% filedir = "results/InterfaceGraph/";
% filedir = "../../results/TrackInterface/Disk4Vortex16/";
filedir = "../../results/TrackInterface/Disk5DeformationT4Order8/";
% filedir = "../../results/TrackInterface/Graph41VortexT8Order4/";
% filedir = "../../test/data1/";
% round = "No4_";

% filename = "spadjor-13.input.dat";
% filename = "localVolumes.dat";
% filename = "localYinset.dat";
filename = "4Circle_grid32_Step512_f.dat";
% filename = "4Circle_grid32.dat";
% filename = "00test.dat";
% filename = "test.dat";
Order = 4;
N = 4;
Shape = "Rose";
% Shape = "Rectangle";
% filename = num2str(Order) + Shape + num2str(N) + "_" + num2str(k) + ".dat";
hd = fopen(filedir + filename);
tensor = true;
volume = false;
if volume
    dat = readCellVolume(hd);
else
    for k = 2:1:2

        if ~tensor
            sf = readYinSet(hd);
            plotYinSet(sf, lineColor, fillColor{k}); hold on
        else
            vecSf = readVecYinSet(hd);
            if k == 1 
                continue;
            end
            for i=2:1:2
                plotYinSet(vecSf{i}, lineColor, fillColor{i}); hold on
            end
        end
    end
end