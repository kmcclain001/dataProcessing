function out = project_to_track(pos,track)
    out = zeros(size(pos,1),1);
    for i=1:size(pos,1)
        out(i) = find_nearest_pixel(pos(i,:),track);
    end
end