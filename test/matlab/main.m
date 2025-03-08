%% MARS plot
% filedir = "results/Vortex/";
% filedir = "results/Deformation/";
% filedir = "results/CircleShrink/";
% filedir = "test/data/";
filedir = "results/CutCell/";
% round = "No4_";
% filename = round + "Start.dat";
% filename = round + "Step160.dat";
% filename = "spadjor-13.input.dat";
Order = 2;
N = 512;
% Shape = "Circle";
Shape = "Rectangle";
filename = num2str(Order) + Shape + num2str(N) + ".dat";

hd = fopen(filedir + filename);

lineColor = '';
fillColor = 'bcrgm';
figure

% sf = readYinSet(hd);
% plotYinSet(sf, lineColor);

vecSf = readVecYinSet(hd);
for i=1:length(vecSf)
    plotYinSet(vecSf{i}, lineColor, fillColor(1)); hold on
end