%% Add usable trials, smooth occupancy & spike count, add fr

function Tuning = fr_maps_session(varargin)

% Parse Inputs

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'filename',[],@isstr);

parse(p,varargin{:});

basepath = p.Results.basepath;
filename = p.Results.filename;

% load shit
behav_file = [basepath filesep filename '.linear.behavior.mat'];
tuning_file = [basepath filesep filename '.Tuning.cellinfo.mat'];
load(behav_file)
load(tuning_file)

% good trial conditions
usable_types = [];
for j = 1:max(behavior.events.trialConditions)
    if sum(behavior.events.trialConditions==j)>=10
        usable_types = [usable_types j];
    end
end
behavior.usableTypes = usable_types;
Tuning.usableTypes = usable_types;

% smooth occupancy
Tuning.occupancySmooth = gaussian_smooth_graph(behavior.trackGraph,Tuning.occupancy,4);

% smooth spike counts
smoothed_spikes = zeros(size(Tuning.spikeCount));
for j = 1:Tuning.nCells
    if strcmp(Tuning.region{j},'hpc')
        smoothed_spikes(:,:,j) = gaussian_smooth_graph(behavior.trackGraph,Tuning.spikeCount(:,:,j),4);
    end
end
Tuning.spikeCountSmoothed = smoothed_spikes;

% add fr
Tuning.fr = Tuning.spikeCountSmoothed./rep_zero(Tuning.occupancySmooth);

% add rate maps
for j = 1:max(Tuning.trialType)
    trial_inds = find(Tuning.trialType==j);
    pos_inds = behavior.events.mapLinear{j};
    rate_maps{j} = zeros(Tuning.nCells,length(pos_inds));
    se_maps{j} = zeros(Tuning.nCells,length(pos_inds));
    for k = 1:Tuning.nCells
        rate_maps{j}(k,:) = mean(Tuning.fr(trial_inds,pos_inds,k),1);
        se_maps{j}(k,:) = std(Tuning.fr(trial_inds,pos_inds,k),1)/(length(trial_inds).^.5);
    end
end
Tuning.rateMaps = rate_maps;
Tuning.seMaps = se_maps;
    
% save shit
save(behav_file,'behavior');
save(tuning_file,'Tuning');

end