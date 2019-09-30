% Classify cell types across all sessions

load('CellClassificationParams.mat');

feature_means = mean(Data);
feature_stds = std(Data);

%load list of experimental sessions used
load('session_names.mat');

for j = 1:length(used_names)
    
    disp(j)
    
    basepath = [DATASET DIRECTORY HERE]; %assumes datasets are organized by session
    cell_class_file = [basepath filesep used_names{j} '.CellClass.cellinfo.mat'];
    old_file = [basepath filesep used_names{j} '.old1CellClass.cellinfo.mat'];
    if exist(cell_class_file,'file')
        load(cell_class_file);
        save(old_file,'CellClass');
    end
        
    spikes = bz_GetSpikes('basepath',basepath);
    CellClass = classify_cell_types(spikes,centroids,feature_means,feature_stds);
    save(cell_class_file,'CellClass')
end