function figS2_moving_target(pubDir)
% FIGS2_MOVING_TARGET  Fig. S2: with the target moving, do populations whose responsiveness was
% lowered restore it? Responsiveness relative to the control over the 150 generations after the
% split (generation 350), for the control and the lowered arm (medians, interquartile bands).
%
% Produces : Fig. S2
% Inputs   : pubDir (optional; default the publication folder under the data root)
%            publication/tracking_moving_target/moving_target.mat and moving_target.txt (written by
%              run_tracking_moving_target), and the two held-still runs fig6_figS2_data reads
%              and checks before it draws (tracking_cut_after_150/design1, design2)
% Outputs  : output/figures/FigS2.pdf (vector, 8.7 cm wide), .png (300 dpi) and .tif (1200 dpi,
%              LZW), under the release root; the checks, printed
% Run      : figS2_moving_target
% Seed     : none (no random draws)
%
% The drawing is fig6_figS2_data's style 'si_moving', written to a temporary folder that this file
% copies from and then deletes.

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1, pubDir = fullfile(grn_data_file('publication')); end
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    scratch = tempname; mkdir(scratch);
    cleanup = onCleanup(@() rmdir(scratch, 's')); %#ok<NASGU>

    fig6_figS2_data(struct('pubDir', pubDir, 'figDir', scratch, 'style', 'si_moving'));
    outBase = fullfile(outDir, 'FigS2');
    copyfile(fullfile(scratch, 'fig6_figS2_data_moving.pdf'), [outBase '.pdf']);
    copyfile(fullfile(scratch, 'fig6_figS2_data_moving.png'), [outBase '.png']);
    copyfile(fullfile(scratch, 'fig6_figS2_data_moving.tif'), [outBase '.tif']);
    fprintf('figS2_moving_target: wrote %s.{pdf,png,tif}\n', outBase);
end
