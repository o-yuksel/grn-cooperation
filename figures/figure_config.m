function cfg = figure_config()
% FIGURE_CONFIG  The folders and the export resolution that several figure and analysis scripts
% take from one place.
%
% Implements : the paths of fig2_draw and fig2_data (Fig. 2), fig3_draw and fig3_data (Fig. 3),
%              fig2_signaling_architecture (Fig. 2) and tableS6_per_mutation_shares (Table S6)
% Inputs     : none
% Outputs    : cfg struct: projectRoot (the release root), analysisDir (this file's folder,
%              figures/), dataDir (the data root, from grn_data_file), figDir (output/figures under
%              the release root; the figure scripts create it when they export), dpi (300)
% Seed       : none (no random draws)

    % Folder of this file (figures/ in the release)
    thisDir = fileparts(mfilename('fullpath'));
    if isempty(thisDir)
        thisDir = pwd;
    end

    cfg.projectRoot = fullfile(thisDir, '..');
    cfg.analysisDir = thisDir;
    cfg.dataDir     = fileparts(grn_data_file('publication'));
    cfg.figDir      = fullfile(cfg.projectRoot, 'output', 'figures');

    % Export settings
    cfg.dpi = 300;
end
