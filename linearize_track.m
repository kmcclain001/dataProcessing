% Compute linear track from set of trajectories
%   Skeletonize occupancy map to get a graph representation of track
%   inputs:
%       filename - name of behavior file (buzcode format)
%       plot info - boolean of whether you want track plotted
%
%   outputs:
%       f - function that takes in a set of x, y positions and outputs
%       linear positions
%       G - a graph that represents the track as a set of connected
%       nodeswhere each node is a spatial bins
%       fig - plotted linear representations
%
%   Author: Kathryn McClain

function [f, G, fig] = linearize_track(filename,plot_info)

    % Load Data
    load(filename)
    n_trials = size(behavior.events.trials,2);
    trials = behavior.events.trials;
    x = cellfun(@(S) S.x, trials,'uni',false);
    y = cellfun(@(S) S.y, trials,'uni',false);
    x_all = vertcat(x{:});
    y_all = vertcat(y{:});

    % Transform to Standard Map
    west = min(x_all);
    east = max(x_all);
    north = max(y_all);
    south = min(y_all);
    grain = 100;
    scale = @(p) scale_2D(p,east,west,north,south,grain-1,grain-1);

    t = scale([x_all,y_all]);
    u = t(:,1);
    v = t(:,2);

    % Compute Heat Map over 2D grid
    occ = zeros(grain);
    for i=1:grain
        for j=1:grain
            occ(i,j) = sum(v==i-1&u==j-1);
        end
    end

    % Rectify
    %thresh = .3*mean((occ(occ>0)));
    c = occ.^2;
    thresh = median(c(c>0));
    occ_rect = occ;
    occ_rect(c<thresh) = 0;
    
    % scale
    occ_sc = occ_rect.^.25;
    
    % Embedd
    border = floor(grain/4);
    occ_bord = zeros(grain+(2*border));
    occ_bord((border+1):(grain+border),(border+1):(grain+border)) = occ_sc;
    
    % Smooth Heat Map
    occ_smooth = imgaussfilt(occ_bord,floor(border/6));
    occ_flat = occ_smooth.^.25;

    % Threshold
    %thresh = .25*max(max(occ_smooth));
    thresh=mean(occ_flat(occ_flat>0));
    occ_flat(occ_flat<=thresh) = 0;
    occ_flat(occ_flat>thresh) = 1;

    % Remove Blobs
    temp = bwareaopen(occ_flat,floor(grain/2));
    surface = double(~bwareaopen(~temp,floor(grain/2)));
    
    % Skeletonize and remove border
    skel_img_bord = bwskel(logical(surface),'MinBranchLength', floor(grain/10));
    skel_img = skel_img_bord((border+1):(grain+border),(border+1):(grain+border));

    % List of track points
    p = find(skel_img);
    [px,py] = ind2sub(size(skel_img),p);
    f = @(p) project_to_track(scale(p)+1,[px,py]);

    % Make Graph
    adj = bw_adjacency_matrix(skel_img);
    G = graph(adj,'OmitSelfLoops');

    % Order Points in sensible way
    start = f([x_all(1),y_all(1)]);
    v = dfsearch(G,start);
%     px = px(v);
%     py = py(v);
%     f = @(p) project_to_track(scale(p)+1,[px,py]);
    adj = adj(v,v);
    G = graph(adj,'OmitSelfLoops');
    
    % plot info
    if plot_info
        fig = figure();hold on
        subplot(1,3,1);
        imagesc(occ);
        
        subplot(1,3,2);
        imagesc(skel_img);
%         inds = 1:length(v);
%         colored_img = zeros(size(skel_img));
%         for i = 1:length(v)
%             colored_img(py(i),px(i)) = inds(i);
%         end
%         imagesc(colored_img);
        %image(skel_img,'CDataMapping','scaled')
        
        subplot(1,3,3);
        l = f([x_all,y_all]);
        scatter(x_all,y_all,[],l,'.');
    else
        fig = [];
    end
        
end



