% Compute firing maps by trial
% each cell has a matrix of nPos X nTrials
% and each entry is the number of spikes at that position
% then theres an overarching occupancy map that is
% nPos x nTrials
% Could also be helpful to have location at which each spike occured

function Tuning = compute_tuning1(varargin)

% Parse Inputs

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'filename',[],@isstr);
addParameter(p,'spikeFile',[],@isstr);

parse(p,varargin{:});

basepath = p.Results.basepath;
force_reload = p.Results.forceReload;
filename = p.Results.filename;
spikeFile = p.Results.spikeFile;

if isempty(spikeFile)
    spikeFile = basepath;
end

if isempty(filename)
    sessionInfo = bz_getSessionInfo(basepath,'noprompts',true);
    filename = sessionInfo.FileName;
end

% Load spikes and behavior

savefile = [basepath filesep filename '.Tuning.cellinfo.mat'];
behavfile = [basepath filesep filename '.linear.behavior.mat'];

if ~exist(behavfile, 'file')
    disp('There is not a linear behavior file in basepath.')
    Tuning = [];
    return
end

if exist(savefile,'file') && ~force_reload
    disp([filename ' already had a tuning file']);
    load(savefile)
else

spikes = bz_GetSpikes('basepath',spikeFile);
load(behavfile);

Tuning.UID = spikes.UID;
Tuning.nTrials = length(behavior.events.trials);
Tuning.nCells = length(spikes.UID);
Tuning.nPos = height(behavior.trackGraph.Nodes);
Tuning.occupancy = zeros(Tuning.nTrials,Tuning.nPos);
Tuning.spikeCount = zeros(Tuning.nTrials,Tuning.nPos,Tuning.nCells);
Tuning.trialType = behavior.events.trialConditions;
Tuning.region = spikes.region;
 
for i = 1:Tuning.nTrials
    traj = behavior.events.trials{i}.l;
    for j = 1:length(traj)
        Tuning.occupancy(i,traj(j)) = Tuning.occupancy(i,traj(j))+1;
        t1 = behavior.events.trials{i}.timestamps(j);
        t2 = t1 + 1/behavior.samplingRate;
        for k = 1:Tuning.nCells
            nSpikes = sum(spikes.times{k}>=t1&spikes.times{k}<t2);
            Tuning.spikeCount(i,traj(j),k) = Tuning.spikeCount(i,traj(j),k) + nSpikes;
        end
    end
end

% Convert occupancy to seconds
Tuning.occupancy = Tuning.occupancy./behavior.samplingRate;

% Add description of behavior
Tuning.description = behavior.description;

save(savefile,'Tuning')
end

end
    