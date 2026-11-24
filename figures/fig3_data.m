function varargout = fig3_data(what, varargin)
% FIG3_DATA  The data of Fig. 3: the success rates and shares of members of the signaling
% architectures and of the spatial settings.
%
% Implements : the data of Fig. 3 (request 'data', called by fig3_draw)
% Inputs     : what   'data'; cfg (default figure_config()), whose dataDir (the data root)
%                     holds:
%              publication/css_predicts_sensing_gen500.mat   (condition_names, success_empirical,
%                reach_occupancy = the share of members; written by gen500_share_of_members)
%              publication/css_predicts_spatial.mat   (written by gen500_share_of_members_spatial)
% Outputs    : OD: the success rates (percent) and shares of members of the architectures (sens)
%              and the spatial settings (spat), and the architecture palette; prints the counts
% Seed       : none (no random draws)
%
% fig3_draw adds the group points and computes r over all points.

    switch what
        case 'data', [varargout{1:nargout}] = share_data(varargin{:});
        otherwise, error('fig3_data:request', 'unknown request: %s', what);
    end
end


function OD = share_data(cfg)
%SHARE_DATA  the two stored inputs and the palette.
    if nargin < 1 || isempty(cfg), cfg = figure_config(); end
    pubDir = fullfile(cfg.dataDir, 'publication');

    OD.sens = load(fullfile(pubDir, 'css_predicts_sensing_gen500.mat'), ...
        'condition_names', 'success_empirical', 'reach_occupancy');
    OD.spat = load(fullfile(pubDir, 'css_predicts_spatial.mat'), ...
        'reach_occupancy', 'success_empirical');
    OD.palette = figure_utils.architecture_palette();
    fprintf('fig3_data: %d architectures, %d spatial settings\n', ...
        numel(OD.sens.reach_occupancy), numel(OD.spat.reach_occupancy));
end
