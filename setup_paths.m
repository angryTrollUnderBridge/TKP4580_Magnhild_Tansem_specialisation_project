function setup_paths()
% setup_paths  Put the project's code on the MATLAB path (run once per session).
    root = fileparts(mfilename('fullpath'));
    addpath(root);
    addpath(genpath(fullfile(root, 'src')));
end
