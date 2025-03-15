function dat = readCellVolume(hd)

tmp = fread(hd, 2, 'int');

dat = fread(hd, [tmp(1) tmp(2)], 'double');


end