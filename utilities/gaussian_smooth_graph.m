function v_smooth = gaussian_smooth_graph(G,v,s)

    %rows of v get smoothed, each row is like a trial or something
    
    v_smooth = zeros(size(v));
    
    %only compute for rows that are not all zero
    rows = find(sum(v,2));
    
    kernel = 0:(2*s);
    kernel = (2*pi*s^2)^-.5*exp(-((kernel/s).^2)/2);
    
    for i=1:size(v,2)
        nodes = neighbors_of_neighbors(G,i,2*s);
        nodes = vertcat([i,0],nodes);
        weights = kernel(nodes(:,2)+1);
        v_smooth(rows,i) = weights*v(rows,nodes(:,1))';
    end
    
    v_smooth = v_smooth.*(sum(v,2)./rep_zero(sum(v_smooth,2)));
end