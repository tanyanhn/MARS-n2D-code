%% MARS plot
close all
clear
lineColor = '';
fillColor = 'bcrgmy';
% filedir = "results/CutCell/";
% filedir = "results/InterfaceGraph/";
filedir = "../../results/TrackInterface/Disk4Vortex8/";
% round = "No4_";

% filename = round + "Start.dat";
% filename = round + "Step160.dat";
% filename = "spadjor-13.input.dat";
% filename = "localVolumes.dat";
filename = "localYinset.dat";

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
    figure
    for k = 3:1:3

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