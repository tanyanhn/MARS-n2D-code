%% MARS plot
% filedir = "results/Vortex/";
% filedir = "results/Deformation/";
filedir = "results/CircleShrink/";
round = "No4_";
filename = round + "Start.dat";
% filename = round + "Step160.dat";

hd = fopen(filedir + filename);

lineColor = '';

sf = readYinSet(hd);
figure
plotYinSet(sf, lineColor);