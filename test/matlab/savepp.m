function savepp(handle, pp)
% save the pp form to the file 'outfile'
    fwrite(handle, [pp.dim pp.order pp.pieces], 'int');
    fwrite(handle, pp.breaks, 'double');
    tmp = pp.coefs(:,end:-1:1);
    fwrite(handle, tmp', 'double');
end

