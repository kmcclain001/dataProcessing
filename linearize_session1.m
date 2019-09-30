% create behavior file with relevent linearized info
% and speed info

function behavior = linearize_session1(varargin)

% Parse Inputs

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'speedSmoothing',5/6,@isfloat);
addParameter(p,'filename',[],@isstr);

parse(p,varargin{:});

basepath = p.Results.basepath;
force_reload = p.Results.forceReload;
smoothing_param = p.Results.speedSmoothing;
filename = p.Results.filename;

savefile = [basepath filesep filename '.linear.behavior.mat'];
oldfile = [basepath filesep filename '.behavior.mat'];

% Check what files exist
if ~exist(oldfile,'file')
    disp([filename ' has no original behavior file']);
    behavior = [];
    return
end

if exist(savefile,'file') && ~force_reload
    disp([filename ' already had a linearized behavior file']);
    load(savefile)
    return
end

% Linearize track
[f,G,fig] = linearize_track(oldfile,true);
saveas(fig,[basepath filesep 'linearization.fig']);
saveas(fig,[basepath filesep 'linearization.png']);

% Load original behavior file
load(oldfile);
n_trials = size(behavior.events.trials,2);
behavior.trackGraph = G;
behavior.trackDistanceMatrix = distances(G);

% Add description of behavior
if ~isfield(behavior,'description')
    behavior.description = behavior.behaviorinfo.description;
end

% SD of gaussian speed filter
sigma = smoothing_param*behavior.samplingRate;
kern = gausswin(round(sigma));
kern = kern./sum(kern);

for i=1:n_trials
    
    %Compute linearization for trial
    x = behavior.events.trials{i}.x;
    y = behavior.events.trials{i}.y;
    behavior.events.trials{i}.l = f([x,y]);
    
    %Compute speed for trial
    dx = diff(x);
    dy = diff(y);
    ds_raw = (dx.^2+dy.^2).^.5;
    ds = conv(ds_raw,kern,'same');
    behavior.events.trials{i}.speedArb = interp1(1:length(ds),ds,0:length(ds),'cubic');
end

subsamp = round(behavior.samplingRate/4);

% if map has already been computed just linearize that
if isfield(behavior.events,'map')
    for i = 1:max(behavior.events.trialConditions)
        l = f([behavior.events.map{i}.x,behavior.events.map{i}.y]);
        lsub = [l(1:subsamp:end); l(end);l(1)];
        P = [];
        for j = 1:length(lsub)-1
            P = [P shortestpath(behavior.trackGraph,lsub(j),lsub(j+1))];
        end
        P = unique(P,'stable');
        behavior.events.mapLinear{i} = P;
    end
else
    
    % subsample each trial trajectory then average, then find linearized
    % position and do graph search from piont to point
    trial_pos_subsamps = zeros(n_trials,10,2);
    for i = 1:n_trials
        if isempty(behavior.events.trials{i}.x)
            continue
        end
        sub_inds = round(linspace(1,length(behavior.events.trials{i}.x),10));
        trial_pos_subsamps(i,:,1) = behavior.events.trials{i}.x(sub_inds);
        trial_pos_subsamps(i,:,2) = behavior.events.trials{i}.y(sub_inds);
    end
    for j = 1:max(behavior.events.trialConditions)
        trial_inds = find(behavior.events.trialConditions==j);
        mean_sub_traj = mean(trial_pos_subsamps(trial_inds,:,:),1);
        lin_sub_traj = f(reshape(mean_sub_traj,10,2));
        P = [];
        for k = 1:length(lin_sub_traj)-1
            P = [P shortestpath(behavior.trackGraph,lin_sub_traj(k),lin_sub_traj(k+1))];
        end
        P = unique(P,'stable');
        behavior.events.mapLinear{j} = P;
    end
end
            
behavior.linearFunction = f;

save(savefile,'behavior');

end




