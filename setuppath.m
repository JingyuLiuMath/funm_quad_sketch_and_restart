function setuppath()

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));
addpath(genpath(fullfile(root, 'test')));
addpath(genpath(fullfile(root, 'experiments')));

end
