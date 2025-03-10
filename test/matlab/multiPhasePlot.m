%% MARS plot
close all
lineColor = '';
fillColor = 'bcrgmy';
% filedir = "results/CutCell/";
filedir = "results/InterfaceGraph/";
% round = "No4_";
% filename = round + "Start.dat";
% filename = round + "Step160.dat";
% filename = "spadjor-13.input.dat";
Order = 4;
N = 4;
Shape = "Circle";
% Shape = "Rectangle";
tensor = true;
figure
for k = 1:1:6
filename = num2str(Order) + Shape + num2str(N) + "_" + num2str(k) + ".dat";

hd = fopen(filedir + filename);

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