function fig3_success_vs_share_of_members()
% FIG3_SUCCESS_VS_SHARE_OF_MEMBERS  Draws Fig. 3: success against the generation-500 share of
% members in the 11 settings (four signaling architectures, five spatial settings, and groups
% split by payoff or at random).
%
% Produces : Fig. 3
% Inputs   : read by fig3_draw and fig3_data:
%            publication/css_predicts_sensing_gen500.mat   (written by gen500_share_of_members)
%            publication/css_predicts_spatial.mat   (written by gen500_share_of_members_spatial)
%            publication/group_two_level/group_fission_analysis_N1000.mat   (written by
%              group_selection_arms)
% Outputs  : output/figures/Fig3.pdf (vector), .png (300 dpi) and .tif (600 dpi, LZW), under the
%            release root; prints each point with its 95% interval
% Run      : fig3_success_vs_share_of_members
% Seed     : none (no random draws)
%
% The drawing is fig3_draw's; this file sets the output path.

    root = fileparts(fileparts(mfilename('fullpath')));
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outBase = fullfile(outDir, 'Fig3');

    fig3_draw(struct('outStem', outBase));
    fprintf('fig3_success_vs_share_of_members: wrote %s.{pdf,png,tif}\n', outBase);
end
