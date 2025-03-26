%% MARS plot
close all
% clear
% hold on
figure
lineColor = '';
fillColor = 'bcrgmy';
% filedir = "results/CutCell/";
% filedir = "results/InterfaceGraph/";
% filedir = "../../results/TrackInterface/Disk4Vortex16/";
% filedir = "../../results/TrackInterface/Disk5Deformation4/";
filedir = "../../results/TrackInterface/Graph41VortexT4Order4/";
% filedir = "../../test/data1/";
% round = "No4_";

% filename = "spadjor-13.input.dat";
% filename = "localVolumes.dat";
% filename = "localYinset.dat";
% filename = "4Circle_grid32_Step2048_4_n.dat";
% filename = "4Circle_grid32.dat";
filename = "test.dat";
Order = 4;
N = 4;
Shape = "Rose";
% Shape = "Rectangle";
% filename = num2str(Order) + Shape + num2str(N) + "_" + num2str(k) + ".dat";
hd = fopen(filedir + filename);
tensor = false;
volume = false;
if volume
    dat = readCellVolume(hd);
else
    for k = 4:1:5

        if ~tensor
            sf = readYinSet(hd);
            plotYinSet(sf, lineColor, fillColor(k)); hold on
        else
            vecSf = readVecYinSet(hd);
            for i=1:length(vecSf)
                plotYinSet(vecSf{i}, lineColor, fillColor(k)); hold on
            end
        end
    end
end