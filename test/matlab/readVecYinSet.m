function vecSf = readVecYinSet(hd)
num = fread(hd, 1, 'int');
vecSf = cell(num);
for i = 1:num
    vecSf{i} = readYinSet(hd);
end
end