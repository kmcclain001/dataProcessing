function Tuning = compute_place_fields_session1(varargin)

% Parse Inputs

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'filename',[],@isstr);

parse(p,varargin{:});

basepath = p.Results.basepath;
filename = p.Results.filename;

% Load tuning cell type and behavior
tuningfile = [basepath filesep filename '.Tuning.cellinfo.mat'];
constuningfile = [basepath filesep filename '.consTuning.cellinfo.mat'];
behavfile = [basepath filesep filename '.linear.behavior.mat'];
cellfile = [basepath filesep filename '.CellClass.cellinfo.mat'];

if ~exist(tuningfile,'file')
    disp('There is not a tuning file in basepath.')
    Tuning = [];
    return
end

if ~exist(behavfile, 'file')
    disp('There is not a linear behavior file in basepath.')
    Tuning = [];
    return
end

if ~exist(cellfile,'file')
    disp('There is not a cell type file in basepath.')
    Tuning = [];
    return
end

load(tuningfile)
load(behavfile);
load(cellfile);

G = behavior.trackGraph;
hasField = zeros(Tuning.nCells,1);

% erase old place field identification is it exists
if isfield(Tuning,'placeFields')
    Tuning = rmfield(Tuning,'placeFields');
end

% For each trial type
for i=1:max(Tuning.trialType)
    
    if ismember(i,Tuning.usableTypes)
        
        % compute firing rate
        position_type_inds = behavior.events.mapLinear{i};
        trial_type_inds = Tuning.trialType==i;
        subG = G.subgraph(position_type_inds);
        FR = Tuning.fr(trial_type_inds,position_type_inds,:);
        
        % for each cell
         for j=1:Tuning.nCells
            
            if strcmp(Tuning.region{j},'hpc') & CellClass.Pyr(j)==1
                
                % find fields
                [node_inds,field_labels] = find_place_fields1(subG,FR(:,:,j),...
                    'max_frac',.2,'FR_thresh',5,'spatial_co_thresh',-1,'smoothed',Tuning.rateMaps{i}(j,:));
                % This part of structure is kind of gross:
                % Tuning > cell index > trial group > fields and labels
                Tuning.placeFields{j}.trialType{i}.fieldInds = node_inds;
                Tuning.placeFields{j}.trialType{i}.fieldLabel = field_labels;
            else
                Tuning.placeFields{j}.trialType{i}.fieldInds = [];
                Tuning.placeFields{j}.trialType{i}.fieldLabel = [];
            end
        end
    else
        for j=1:Tuning.nCells
            Tuning.placeFields{j}.trialType{i}.fieldInds = [];
            Tuning.placeFields{j}.trialType{i}.fieldLabel = [];
        end
    end
end

% make list of cells with fields
for j = 1:Tuning.nCells
    field = 0;
    for k = 1:max(Tuning.trialType)
        if ~isempty(Tuning.placeFields{j}.trialType{k}.fieldInds)
            field = 1;
        end
    end
    Tuning.placeFields{j}.hasField = field;
    hasField(j) = field;
end
Tuning.hasField = hasField;

save(constuningfile,'Tuning')
end