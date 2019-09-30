%% collect relevant info for fitting model for each place field
%   need position, phase and spike count arrays, could have timestamp &
%   speed too
%   for each field, have cell of trials, for each trial* have vector of each
%   *make sure only include visits to field that are real
%
% UPDATED from scratch_st_model_object with smoothing of position

basename = [DATASET DIRECTORY HERE]; %assumes datasets are organized by session

%load list of experimental sessions used
load('session_names.mat');

%%% load relavent stuff %%%
for sesh = 1:length(used_names)

filename = used_names{sesh};

%get channel number
try
load([basename filesep filename filesep filename '.sessionInfo.mat']);
chan = sessionInfo.ca1;

%load tuning 
load([basename filesep filename filesep filename '.consTuning.cellinfo.mat'])

%load cell type info
load([basename filesep filename filesep filename '.CellClass.cellinfo.mat']);

%get behavior interval
load([basename filesep filename filesep filename '.linear.behavior.mat']);

%load theta info
load([basename filesep filename filesep filename '.thetaInfo.mat'])

%get spikes
spikes = bz_GetSpikes('basepath',[basename filesep filename]);

%make filter
passband = [4 15];
lfp_rate = sessionInfo.rates.lfp;
[b,a] = butter(4,[passband(1)/(lfp_rate/2) passband(2)/(lfp_rate/2)],'bandpass');

%new object
clear stModel
directory = [];

catch
    disp([filename ' no go'])
    continue
    
end

count = 1;

%%% iterate through cells
for cell = 1:Tuning.nCells
    
    % is it not in ca1, pyramidal, or has field, skip
    if ~(strcmp(Tuning.region{cell},'hpc') && CellClass.PyrInt(cell)==1 && Tuning.hasField(cell)==1)
        continue
    end
        
    nTypes = max(Tuning.trialType);
    
    %%% iterating through trial types
    for tr_type = 1:nTypes
        
        trial_inds = find(Tuning.trialType==tr_type);
        track = behavior.events.mapLinear{tr_type};
        
        if ~ismember(tr_type,Tuning.usableTypes)
            continue
        end
        
        fields = Tuning.placeFields{cell}.trialType{tr_type};
        all_field_inds = fields.fieldInds;
        all_field_labels = fields.fieldLabel;
        
        %if cell doesn't have field in this trial type, skip
        if isempty(all_field_inds)
            continue
        end
        
        %%% iterate through fields in this trial type
        for label = 1:max(all_field_labels)
            
            
            %location of this field on track
            field_inds = all_field_inds(all_field_labels==label);
            field_inds_orig = track(field_inds);
            
            % fast time scale values
            timestamps = [];
            spike_count = [];
            position = [];
            theta_phase = [];
            speed = [];
            
            % slow time scale values
            behav_time_stamps = [];
            behav_pos = [];
            behav_speed = [];
            
            data = [];
            %%% iterate through releveant trials
            for t = 1:length(trial_inds)
                
                trial_num = trial_inds(t);
                         
                %trim position on track
                pos_norm = transform_to_trial_frame(behavior.events.trials{trial_num}.l,track,behavior.trackDistanceMatrix);
                start_ind = min(find(pos_norm>(min(pos_norm)+.001)&pos_norm<.5));
                end_ind = max(find(pos_norm<(max(pos_norm)-.001)&pos_norm>.5));
                good_trial_inds = start_ind:end_ind;
                pos_on_track = behavior.events.trials{trial_num}.l(good_trial_inds);
                         
                %find points when in field
                in_field = ismember(pos_on_track,field_inds_orig);
                enter_field = find(diff([0;in_field])==1);
                exit_field = find(diff([in_field;0])==-1);
                field_int = zeros(length(enter_field),2);
                field_int(:,1) = behavior.events.trials{trial_num}.timestamps(enter_field);
                field_int(:,2) = behavior.events.trials{trial_num}.timestamps(exit_field);
                         
                %find largest interval in field
                [~,max_field_visit] = max(diff(field_int,1,2));
                idx1 = enter_field(max_field_visit);
                idx2 = exit_field(max_field_visit);
                pos_in_field = pos_on_track(idx1:idx2);
                         
                %ignore this visit if only in the field very breifly
                if sum(ismember(field_inds_orig,pos_in_field))<.4*length(field_inds_orig)
                    continue
                end
                         
                %get slow sampling variables
                slow_idx = good_trial_inds(idx1:idx2); %index of original timestamp
                slow_timestamps = behavior.events.trials{trial_num}.timestamps(slow_idx);
                slow_pos = transform_to_trial_frame(pos_in_field,field_inds_orig,behavior.trackDistanceMatrix);
                kern = gausswin(5,.5); %added
                kern = kern/sum(kern); %this
                slow_pos =  conv([0;0; slow_pos; 1;1],kern,'valid'); %part
                slow_speed = behavior.events.trials{trial_num}.speed(slow_idx);
                         
                %get phase
                interval = [slow_timestamps(1),slow_timestamps(end)];
                phase_data = thetaInfo.trial{trial_num}.data;
                in_field = find(phase_data(:,1)>=interval(1)&phase_data(:,1)<=interval(2));
                phase = phase_data(in_field,5); %using peak to peak phase measure
                         
                %get fast sampling variables
                fast_timestamps = phase_data(in_field,1);
                fast_pos = interp1(slow_timestamps,slow_pos,fast_timestamps,'spline');
                fast_speed = interp1(slow_timestamps,slow_speed,fast_timestamps,'spline');
                         
                %get spike count
                spike_count = zeros(size(fast_timestamps));
                rel_spikes = spikes.times{cell}(spikes.times{cell}>interval(1)&spikes.times{cell}<=interval(2));
                for s = 1:length(rel_spikes)
                    [~,ind] = min(abs(fast_timestamps-rel_spikes(s)));
                    spike_count(ind) = spike_count(ind)+1;
                end
                         
                %assemble data to matrix
                data = zeros(length(fast_timestamps),5);
                data(:,1) = fast_timestamps;
                data(:,2) = spike_count;
                data(:,3) = fast_pos;
                data(:,4) = phase;
                data(:,5) = fast_speed;
                         
                %save everything I want about the field
                stModel.field{count}.cell = cell;
                stModel.field{count}.trialType = tr_type;
                stModel.field{count}.label = label;
                stModel.field{count}.trial{t} = data;
                stModel.field{count}.fieldInds = field_inds;
                stModel.field{count}.fieldIndsOrig = field_inds_orig;
            end
            
            if ~isempty(data)
            
                count = count + 1;
                
                %record cell, trial type, label of this field
                directory = vertcat(directory, [cell tr_type label]);
            end
        end
    end
end

stModel.directory = directory;
stModel.nField = size(directory,1);
if isempty(directory)
    disp([filename ' has no place cells'])
else
    stModel.dt = diff(fast_timestamps(1:2));
end
save([basename filesep filename filesep filename '.stModel.fieldinfo.mat'],'stModel');
disp(sesh)

end


