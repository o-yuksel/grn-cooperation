%STARTUP  Put the release code on the MATLAB path. Run once per session.
%
%   Adds the code folders (model, shared, experiments, the six analysis subfolders, figures) and
%   the release root (grn_data_file, publication_params). Result files are
%   found through grn_data_file: set the environment variable GRN_DATA to
%   the folder holding the data archive, or unpack it into <release root>/data.

releaseRoot = fileparts(mfilename('fullpath'));
addpath(releaseRoot);
addpath(fullfile(releaseRoot, 'model'));
addpath(fullfile(releaseRoot, 'shared'));
addpath(fullfile(releaseRoot, 'experiments'));
analysisDirs = {'methods', 's1_congruence_identity', 's2_share_of_members', 's3_one_mutation', ...
    's4_losing_cooperation', 's5_tracking_task'};
for k = 1:numel(analysisDirs)
    addpath(fullfile(releaseRoot, 'analysis', analysisDirs{k}));
end
addpath(fullfile(releaseRoot, 'figures'));
clear releaseRoot analysisDirs k
