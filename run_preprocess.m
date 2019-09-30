%Runs master script PREPROCESS over all datasets

load('session_names.mat')

for i= 1:length(used_names)
    disp(used_names{i})
    basepath = [DATASET DIRECTORY HERE];
    preprocess('basepath',basepath,'filename',used_names{i});
end