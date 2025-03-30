function writepp(handle, xsp, ysp)
    tmp = [xsp.dim, xsp.order, xsp.pieces];
    fwrite(handle, tmp, 'int');
    fwrite(handle, xsp.breaks, 'double');
    for k= 1:xsp.pieces
        for d = 1:xsp.order
            fwrite(handle, xsp.coefs(k, :), 'double');
            fwrite(handle, ysp.coefs(k, :), 'double');
        end
    end
end