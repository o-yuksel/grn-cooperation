function figS1_repair_test(pubDir)
% FIGS1_REPAIR_TEST  Draws Fig. S1: in 65 populations about to fail, restarted just before the
% fall, the share of runs that collapse within 100 generations by arm (observed, sham repair,
% repair), with and without the defector that later took over added (A); and the cumulative
% collapse over those generations with the defector added (B).
%
% Produces : Fig. S1
% Inputs   : pubDir (optional; default the publication folder under the data root)
%            publication/repair_matched/repair_matched_build_targets.mat and repair_matched_main.mat
%              (written by run_repair_test)
%            publication/repair_test_observed/repair_test_observed.mat   (written by
%              run_repair_test_observed)
%            and the inputs of repair_test_summary('main'), which figS1_draw calls
% Outputs  : output/figures/FigS1.pdf (vector, 17.8 x 5.4 cm), .png (300 dpi) and .tif (1200 dpi,
%            LZW), under the release root; prints the levels, differences
%            and checks; repair_test_summary rewrites publication/repair_matched/repair_matched.txt
% Run      : figS1_repair_test   (about 1 min, the summary's run time)
% Seed     : rng(20261013, 'twister'), then randi(65, 65, 10000): the summary's resamples of the
%            65 populations (set in figS1_draw and in repair_test_summary)
%
% The drawing is figS1_draw's, written to a temporary folder, copied from there, then deleted.

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1, pubDir = fullfile(grn_data_file('publication')); end
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    scratch = tempname; mkdir(scratch);
    cleanup = onCleanup(@() rmdir(scratch, 's')); %#ok<NASGU>

    figS1_draw(struct('pubDir', pubDir, 'figDir', scratch));
    outBase = fullfile(outDir, 'FigS1');
    copyfile(fullfile(scratch, 'figS1_draw.pdf'), [outBase '.pdf']);
    copyfile(fullfile(scratch, 'figS1_draw.png'), [outBase '.png']);
    copyfile(fullfile(scratch, 'figS1_draw.tif'), [outBase '.tif']);
    fprintf('figS1_repair_test: wrote %s.{pdf,png,tif}\n', outBase);
end
