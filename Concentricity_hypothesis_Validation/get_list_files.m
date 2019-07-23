% A very simple function to get the names of all files in a directory
function out=get_list_files(path,type)
% Get all files
list_dir=dir(fullfile(path,type));
% Extract just the names of the files.
out={list_dir.name};
