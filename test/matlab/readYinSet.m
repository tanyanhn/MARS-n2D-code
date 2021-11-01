function sf = readYinSet(hd)
    nl = fread(hd, 1, 'int');
    sf = cell(1,nl);
    for i=1:nl
        sf{i} = readpp(hd);
    end
end