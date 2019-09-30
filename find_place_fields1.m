% Function for finding place fields on graph representation of environment
%
% Inputs:
%   G = graph of environment (found in .linear.behavior.mat files)
%   FR = matrix of firing rates, column corresponds to node index in G, row
%       corresponds to single trial
%   max_frac = fraction of peak height FR must be above
%   area_thresh = number of contiguous bins required
%   FR_thresh = minimum FR required for peak of each putative field
%   spatial_co_thres = spatial coherence required for field (explained in
%                      spatial_coherence.m)
%
% Outputs:
%   node_inds = list of node indices in environment graph that are part of
%               a place field that fits these standards
%   field_labels = list of labels for which field each node is part of,
%                  index corresponds to index in node_inds
%   ex: node_inds = [4,5,6,8,9], field_labels = [1,1,1,2,2] corresponds to
%       nodes 4, 5 and 6 being in field 1 and nodes 8 and 9 in field 2
%
% Criterion based on Mizuseki & Buzsaki 2013 Cell Reports
% Author: Kathryn McClain
%

function [node_inds,field_labels] = find_place_fields1(G,FR,varargin)

    %Parse inputs
    p = inputParser;
    p.addParameter('max_frac',.1,@isfloat);
    p.addParameter('area_thresh',10,@isfloat);
    p.addParameter('FR_thresh',2,@isfloat);
    p.addParameter('spatial_co_thresh',0.7,@isfloat);
    p.addParameter('area_lim', 5/8,@isfloat);
    p.addParameter('smoothed',[],@isfloat);
    p.addParameter('min_trial_prop',1/5,@isfloat);
    
    parse(p,varargin{:});
    max_frac = p.Results.max_frac;
    area_thresh = p.Results.area_thresh;
    FR_thresh = p.Results.FR_thresh;
    spatial_co_thresh = p.Results.spatial_co_thresh;
    area_lim = p.Results.area_lim;
    smoothed = p.Results.smoothed;
    min_trial_prop = p.Results.min_trial_prop;
    
    %Find spots with FR above threshold
    FR_raw = FR;
    if ~isempty(smoothed)
        FR = smoothed;
    else
        FR = gaussian_smooth_graph(G,FR,4);
    end
    max_FR = max(FR);
    above_thresh = FR > (max_frac*max_FR);
    possible_node_inds = find(above_thresh);
    
    %Find spatially contiguous groups they belong to
    G_small = G.subgraph(possible_node_inds);
    groups = G_small.conncomp;
    
    %Evaluate each group
    fields = unique(groups);
    good = ones(length(fields),1);

    for i = fields
        %Check number of bins
        n_bins = sum(groups==i);
        if n_bins < area_thresh | n_bins > area_lim*G.numnodes
            good(i) = 0;
            continue
        end

        %Check peak FR
        m = max(FR(possible_node_inds(groups==i)));
        if m < FR_thresh
            good(i) = 0;
            continue
        end

        %Check spatial coherence 
        if spatial_co_thresh >-1
            sp_co = spatial_coherence(G,possible_node_inds(groups==i),mean(FR_raw,1),'radius',1);
            if sp_co < spatial_co_thresh
                good(i) = 0;
                continue
            end
        end
        
        %Make sure field is present for certain proportion of trials
        n_trials = size(FR_raw,1);
        rate_in_bins = sum(FR_raw(:,possible_node_inds(groups==i)),2);
        if sum(rate_in_bins>0)< min_trial_prop*n_trials
            good(i) = 0;
        end
        
    end
    
    %Collect final nodes and groups
    n = 0;
    node_inds = [];
    field_labels = [];
    for j = 1:length(good)
        if good(j) == 1
            n = n+1;
            this_group = possible_node_inds(groups==j);
            node_inds = [node_inds this_group];
            field_labels = [field_labels n*ones(1,length(this_group))];
        end
    end
    
end