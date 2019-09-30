function out = neighbors_of_neighbors(G,n,s)

neighbs = [];
for i=1:s
    neighbs = vertcat(neighbs,G.nearest(n,i));
end

nodes = unique(neighbs);
out = zeros(length(nodes),2);
for j=1:length(nodes)
    out(j,1) = nodes(j);
    tmp = sum(neighbs==nodes(j));
    out(j,2) = -tmp + s + 1;
end
end