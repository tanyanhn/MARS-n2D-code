function writepp(handle, xsp, ysp)
    tmp = [xsp.dim, xsp.order, xsp.pieces];
    fwrite(handle, tmp, 'int');
    fwrite(handle, xsp.breaks, 'double');
    xco = xsp.coefs';
    xco = xco(end:-1:1, :);
    yco = ysp.coefs';
    yco = yco(end:-1:1, :);
    for k= 1:xsp.pieces
        fwrite(handle, xco(:, k), 'double');
        fwrite(handle, yco(:, k), 'double');
    end
end