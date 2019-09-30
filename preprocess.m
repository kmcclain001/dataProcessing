%% Master preprocessing script

function preprocess(varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'filename',[],@isstr);

parse(p,varargin{:});

basepath = p.Results.basepath;
filename = p.Results.filename;
spikeFile = p.Results.basepath;

% Linearize session
disp('linearize')
linearize_session1('basepath',basepath,'filename',filename,'forceReload',true);

% Convert speed units
disp('convert speed')
convert_speed_units1('basepath',basepath,'filename',filename);

% Compute tuning of cells
disp('tuning')
compute_tuning1('basepath',basepath,'filename',filename,'forceReload',true,'spikeFile',spikeFile);

% Compute rate maps
disp('rate maps')
fr_maps_session('basepath',basepath,'filename',filename);

% Find place fields
disp('place fields')
compute_place_fields_session1('basepath',basepath,'filename',filename);

end

