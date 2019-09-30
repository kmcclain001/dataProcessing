function pos_norm = transform_to_trial_frame(l,map,distance_matrix)

pos_norm = zeros(size(l));

for i = 1:length(l)
    
    if ismember(l(i),map)
        pos_norm(i) = find(map==l(i))/length(map);
    else
        [~,ind] = min(distance_matrix(l(i),map));
        pos_norm(i) = ind/length(map);
    end
end
end