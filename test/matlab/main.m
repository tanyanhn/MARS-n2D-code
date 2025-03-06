%% MARS plot
% filedir = "results/Vortex/";
% filedir = "results/Deformation/";
% filedir = "results/CircleShrink/";
filedir = "test/data/";
% round = "No4_";
% filename = round + "Start.dat";
% filename = round + "Step160.dat";
filename = "spadjor-13.input.dat";

hd = fopen(filedir + filename);

lineColor = '';

sf = readYinSet(hd);
% figure
plotYinSet(sf, lineColor);