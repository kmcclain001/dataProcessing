
function out = find_nearest_pixel(p,p_list)
    if sum(size(p)==[1,2]) ~= 2
        error('point input needs to be in the form (x,y)')
    elseif size(p_list,2) ~= 2
        error('list of points needs to be in the form (x_array,y_array)')
    end
    
    [~,out] = min(vecnorm(p_list-p,2,2));
end