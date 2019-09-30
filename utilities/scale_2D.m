function out = scale_2D(pnts,maxx,minx,maxy,miny,numx,numy)
    n_pnts = size(pnts,1);
    out = zeros(n_pnts,2);
    out(:,1) = floor(numx/(maxx-minx)*(pnts(:,1)-minx));
    out(:,2) = floor(numy/(maxy-miny)*(pnts(:,2)-miny));
end