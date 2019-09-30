function out = bw_adjacency_matrix(img)

    % indices of all nonzero pixels
    p = find(img);
    [pr,pc] = ind2sub(size(img),p);
    n_pnts = size(pr,1);
    
    % compute adjacency for each point
    out = zeros(n_pnts);
    for i=1:n_pnts
        
        %choose row & column of pixel
        r = pr(i);
        c = pc(i);
        
        %find neighborhoor around pixel
        grid = img(r-1:r+1,c-1:c+1);
        
        %find nonzero pixels in neighborhood
        [nr,nc] = ind2sub(size(grid),find(grid));
        nr = nr-2;
        nc = nc-2;
        
        %convert from neighborhood to full matrix indices
        pr_i = nr+r;
        pc_i = nc+c;
        temp_inds = sub2ind(size(img),pr_i,pc_i);
        
        inds = find(ismember(p,temp_inds));
        out(i,inds) = 1;
        out(inds,i) = 1;
    end
end