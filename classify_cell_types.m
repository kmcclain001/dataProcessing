% Classify CA1 cells using clusters found previously (fast pyr, slow pyr,
% fast int, slow int). Uses fitpyrint for double exponential fitting of
% ACG. The log of the upward and downward slope of the fit, as well as log
% of firing rate and average value of second half of waveform are used as
% the classification features. 
%
% Inputs:
%   spikes = spikes from recording session (use bz_GetSpikes)
%   clusters = centroids of clusters used (saved in
%               CellClassificationParams.mat)
%   feature_means = means from each classification parameter from original
%                   clustering (for normalization)
%   feature_stds = standard deviation from each classification parameter
%                   from original clustering (for normalization)
%
% Outputs:
%   CellClass = cellinfo object with classification parameters and
%               classifications for each object
%
% Author: Kathryn McClain
%

function CellClass = classify_cell_types(spikes, clusters, feature_means, feature_stds)

nCells = length(spikes.times);

for i = 1:nCells
    
    if ~strcmp(spikes.region{i},'hpc')
            CellClass.params{i}.ACG = [];
            CellClass.params{i}.waveform = [];
            CellClass.params{i}.Rate = length(spikes.times{i})/(spikes.times{i}(end)-spikes.times{i}(1));
            CellClass.params{i}.doubexpo = [];
            CellClass.type{i} = 'N/A';
            CellClass.PyrInt(i) = 0;
            CellClass.Cluster(i) = 0;
    else
        cluster_params = zeros(1,4);
        
        % Find wave info
        wave = spikes.rawWaveform{i}';
        wave_scaled = wave-min(wave);
        wave_scaled = wave_scaled./max(wave_scaled);
        [~,ind] = min(wave_scaled);
        last_half = sum(wave_scaled(ind:end))/rep_zero(32-ind);
        cluster_params(4) = last_half;
        
        % Find fr info
        Rate = length(spikes.times{i})/(spikes.times{i}(end)-spikes.times{i}(1));
        cluster_params(1) = log(Rate);
        
        % Find ACG params
        ACGex = CrossCorr(spikes.times{i},spikes.times{i},.001,100);
        ACGex = 10000*ACGex/length(spikes.times{i});
        ACGex(51) = 0;
        [fmodel,ydata,xdata,ACGparams] = fitpyrint(ACGex',0:50,0,10);
        ACG = ACGex(51:100);
        ACG_scaled = ACG/rep_zero(max(ACG));
        cluster_params(2) = log(ACGparams(1));
        cluster_params(3) = log(ACGparams(3));
        
        % Add params to cell info
        CellClass.params{i}.ACG = ACG;
        CellClass.params{i}.waveform = wave;
        CellClass.params{i}.Rate = Rate;
        CellClass.params{i}.doubexpo = ACGparams;
        
        % Find nearest cluster
        features = (cluster_params-feature_means)./feature_stds;
        centroids = (clusters-feature_means)./feature_stds;
        d = sum((centroids-features).^2,2);
        [~,c] = min(d);
        CellClass.Cluster(i) = c;
        
        switch c
            case {5,7}
                CellClass.type{i} = 'fPyr';
                CellClass.PyrInt(i) = 1;
            case {1,3,4,12,11,15}
                CellClass.type{i} = 'sPyr';
                CellClass.PyrInt(i) = 1;
            case {2,8}
                CellClass.type{i} = 'fInt';
                CellClass.PyrInt(i) = -1;
            case {6,10,13}
                CellClass.type{i} = 'sInt';
                CellClass.PyrInt(i) = -1;
            case {9,14}
                CellClass.type{i} = 'Garb';
                CellClass.PyrInt(i) = 0;
        end
    end
end

CellClass.UID = spikes.UID;
CellClass.Pyr = CellClass.PyrInt ==1;
CellClass.Int = CellClass.PyrInt ==-1;


end
