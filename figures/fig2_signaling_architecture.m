function fig2_signaling_architecture(pubDir)
% FIG2_SIGNALING_ARCHITECTURE  Draws Fig. 2: for each signaling architecture, mean cooperation
% over 500 generations with, above it, the direction of selection on cooperation (A); and the
% histogram of generation-500 cooperation across the 100 populations and how many succeed (B).
%
% Produces : Fig. 2
% Inputs   : pubDir (optional; default the publication folder under the data root)
%            publication/<architecture>_100reps.mat   (written by run_signaling_architectures;
%              read by fig2_data)
%            publication/selection_over_time/selection_z_<arch>_r<NN>.mat, selection_sign_<arch>_r<NN>.mat
%              and lineage_records_<arch>_r<NN>.mat, replicates 1 to 20   (written by run_selection_z,
%              run_selection_sign and run_lineage_records)
%            publication/selection_over_time/summary_selection_inherited.txt   (written by
%              selection_inherited; read for the checks, from the same publication folder)
% Outputs  : output/figures/Fig2.pdf (vector), .png (150 dpi) and .tif (600 dpi, LZW), under the
%            release root; prints the checks and the PDF width
% Run      : fig2_signaling_architecture
% Seed     : rng(20261022), default generator, for the bootstrap of the selection strip (fig2_draw)
%
% The drawing and its checks are fig2_draw's; this file sets the output path and stops unless
% the PDF is 17.8 cm wide (two PNAS columns; MediaBox, within 2 pt).

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1, pubDir = fullfile(grn_data_file('publication')); end
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outBase = fullfile(outDir, 'Fig2');
    widthCm = 17.8;

    cfg = figure_config();
    figure_utils.setup_defaults();
    cfg.dataDir = pubDir;
    cfg.outStem = outBase;
    fig = fig2_draw(cfg);
    close(fig);
    figure_utils.check_pdf_width([outBase '.pdf'], widthCm);
    fprintf('fig2_signaling_architecture: wrote %s.{pdf,png,tif}\n', outBase);
end
