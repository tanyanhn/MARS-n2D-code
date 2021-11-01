function pp = ppFromPolygon(polyvert)
% polyvert = [ X Y ] with (xn,yn) ~= (x1,y1)
    n = size(polyvert,1);
    
    pp.form = 'pp';
    pp.pieces = n;
    pp.order = 2;
    pp.dim = 2;
    
    polyvert = [polyvert; polyvert(1,:)];
    
    dt = diff(polyvert, 1, 1);
    dl = sqrt(dt(:,1).^2 + dt(:,2).^2);
    pp.breaks = [0 cumsum(dl)'];
    
    pp.coefs = zeros(n*2, 2);
    tmp = polyvert';
    pp.coefs(:,2) = tmp(1:end-2);
    dt = dt ./ dl;
    dt = dt';
    pp.coefs(:,1) = dt(:);
    
end