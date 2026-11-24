function fig6_tracking_load(pubDir, opts)
% FIG6_TRACKING_LOAD  Fig. 6: (A) measured against predicted per-mutation load in the tracking
% task, target moving filled and held still open; (B) with the target held still, the mutation
% load with responsiveness lowered to 0.70, relative to the control's, against mutation size.
%
% Produces : Fig. 6
% Inputs   : pubDir (optional; default the publication folder under the data root); opts.outBase
%              (optional) replaces the output path stem
%            panel A, through fig6_figS2_data('data') and fig6_load_cells: the per-mutation loads
%              (publication/tracking_runs/tracking_runs_N1000_load.mat and .txt,
%              written by tracking_load_vs_identity) and the other inputs fig6_figS2_data checks
%            panel B: publication/tracking_load_regime/cut_vs_size_summary.txt (the "CUT VS SIZE FINAL" table;
%              written by run_tracking_cut_vs_mutation_size) and
%              each point's main.mat and result text: tracking_load_regime/design1_sd0p006 to design1_sd0p100
%              (run_tracking_cut_vs_mutation_size), tracking_load_regime/design1_sd0p060
%              (run_tracking_cut_vs_shift), tracking_cut_after_150/design1
%              (tracking_cut_after_150) and tracking_mutation_regime/size_10x (run_tracking_mutation_regime)
% Outputs  : output/figures/Fig6.pdf (vector, 8.7 x 10.4 cm), .png (300 dpi) and .tif (1200 dpi,
%              LZW), under the release root; panel B's values and the checks, printed
% Run      : fig6_tracking_load
% Seed     : none here; each bootstrap of fig6_load_cells uses its own fixed RandStream (twister)
%
% Panel B: design 1 (output self-activation a = 0.70) at seven mutation sizes, 0.3 to 10 times the
% standard sd 0.020; y = 1 - cut/100, with the cut and its paired-bootstrap 95% interval from
% cut_vs_size_summary.txt. The run stops unless each table row matches its run's own block and main.mat
% (cut and interval to 0.05 percentage points, sd exactly) and the plotted values equal 1 - cut/100.
% Colors: no cooperation in this task, so ink, greys and white only; panel B's points are open ink
% circles, panel A's mark for "target held still" (every panel B run holds the target still).

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(pubDir), pubDir = fullfile(grn_data_file('publication')); end
    if nargin < 2, opts = struct(); end
    outBase = fullfile(root, 'output', 'figures', 'Fig6');
    if isfield(opts, 'outBase'), outBase = opts.outBase; end
    if ~isempty(fileparts(outBase)) && ~exist(fileparts(outBase), 'dir'), mkdir(fileparts(outBase)); end

    runs = {'tracking_load_regime/design1_sd0p060', 0.70; 'tracking_load_regime/design1_sd0p006', 0.70; ...
        'tracking_load_regime/design1_sd0p010', 0.70; 'tracking_load_regime/design1_sd0p040', 0.70; 'tracking_load_regime/design1_sd0p100', 0.70};
    R = fig6_figS2_data(struct('pubDir', pubDir, 'style', 'data', 'bigGroup', true, 'extraHeldStill', {runs}));
    regime = fullfile(pubDir, 'tracking_load_regime');
    B = curve_data(fullfile(regime, 'cut_vs_size_summary.txt'), R, pubDir);
    draw(R, B, outBase);
end


% ===================== data + controls ===============================
function B = curve_data(file, R, pubDir)
    txt = fileread(file);
    fin = txt(strfind(txt, '==== CUT VS SIZE FINAL'):end);
    B.name = {'0.3x (E)', '0.5x (F)', '1x (d1 sd0.020)', '2x (G)', '3x (d1 sd0.060)', '5x (H)', '10x (d1 sd0.200)'};
    sdWant = [0.006 0.010 0.020 0.040 0.060 0.100 0.200]; sdStd = 0.020;
    runDir = {'tracking_load_regime/design1_sd0p006', 'tracking_load_regime/design1_sd0p010', 'tracking_cut_after_150/design1', ...
        'tracking_load_regime/design1_sd0p040', 'tracking_load_regime/design1_sd0p060', 'tracking_load_regime/design1_sd0p100', 'tracking_mutation_regime/size_10x'};
    run = {'E', 'F', '', 'G', '', 'H', ''};
    % each point's own result (main.mat, unrounded d, through fig6_figS2_data's held_still_design):
    % R.extra = design1_sd0p060, design1_sd0p006, F, G, H (extraHeldStill order); the standard size is the tail-150
    % design-1 run, sd 0.200 the big run
    own = {R.extra(2), R.extra(3), R.held_still(1), R.extra(4), R.extra(1), R.extra(5), R.big};
    n = numel(B.name); B.sd = nan(1, n); B.x = nan(1, n); B.cut = nan(1, n); B.ci = nan(n, 2);
    ok = ~isempty(regexp(fin, '\ndesign 1, target held still, responsiveness x0\.70;', 'once'));   % the table is design 1
    okOwn = true; okX = true;
    for i = 1:n
        pat = ['^' regexptranslate('escape', B.name{i}) '\s{2,}([\d.]+)\s+([\d.]+) \[([\d.]+), ([\d.]+)\]\s+' ...
            '[\d.]+ \[[\d.]+, [\d.]+\]\s+(\d+)\s+\d+\s'];
        tok = regexp(fin, pat, 'tokens', 'lineanchors');
        ok = ok && numel(tok) == 1;
        if numel(tok) ~= 1, continue; end
        v = str2double(tok{1});
        ok = ok && v(5) == 20;                                     % 20 responsiveness-lowered replicates
        B.x(i) = v(1); B.cut(i) = v(2); B.ci(i, :) = v(3:4); B.sd(i) = sdWant(i);
        if ~isempty(run{i})                                        % the run's own block
            blk = regexp(txt, ['==== RUN ' run{i} ' \(.*?(?=\n====|$)'], 'match', 'once');
            hd = regexp(blk, '^==== RUN \w+ \(design 1, a = 0\.70\) \| mutation distribution p 0\.080 sd ([\d.]+) ', 'tokens', 'once');
            t2 = str2double(regexp(blk, 'ALPHA-DOWN cut ([\d.]+)% \[([\d.]+), ([\d.]+)\]', 'tokens', 'once'));
            ok = ok && numel(hd) == 1 && abs(str2double(hd{1}) - sdWant(i)) < 1e-12 ...
                && numel(t2) == 3 && all(abs(t2 - [B.cut(i) B.ci(i, :)]) < 1e-9);
        end
        F = own{i};                                                % the run's own result, 0.05 pp
        okOwn = okOwn && strcmp(F.name, runDir{i}) && F.a == 0.70 ...
            && all(abs(-[F.pct(2) fliplr(F.pctCI(2, :))] - [B.cut(i) B.ci(i, :)]) < 0.05 + 1e-9);
        S = load(fullfile(pubDir, runDir{i}, 'main.mat'), 'out');  % the run's own mutation sd, exact
        if isfield(S.out.cfg, 'mutSd'), sdRun = S.out.cfg.mutSd; else, sdRun = S.out.P.mutation_sd; end
        okX = okX && sdRun == sdWant(i) && abs(B.x(i) * sdStd - sdWant(i)) < 1e-12;
        fprintf('  %s: run sd %.6f, table size %.4g x 0.020 = %.6f\n', runDir{i}, sdRun, B.x(i), B.x(i) * sdStd);
    end
    B.y = 1 - B.cut / 100; B.ylo = 1 - B.ci(:, 2)' / 100; B.yhi = 1 - B.ci(:, 1)' / 100;
    fprintf('PANEL B (design 1, mutation load with responsiveness lowered / control, cut_vs_size_summary.txt):\n');
    for i = 1:n
        fprintf('  %-17s size x %.4g  ratio %.3f [%.3f, %.3f]  cut %.1f%% [%.1f, %.1f]\n', B.name{i}, B.x(i), ...
            B.y(i), B.ylo(i), B.yhi(i), B.cut(i), B.ci(i, :));
    end
    fprintf('  CONTROL table rows (design 1, n 20) and run E-H blocks (design 1, a 0.70, sd, cut): %s\n', figure_utils.pass_fail(ok));
    fprintf('  CONTROL every cut and CI against its run''s own main.mat (0.05 pp): %s\n', figure_utils.pass_fail(okOwn));
    fprintf('  CONTROL every sd against its run''s main.mat and the expected size (exact): %s\n', figure_utils.pass_fail(okX));
    if ~ok, error('fig6_tracking_load:control', 'panel B data do not reproduce cut_vs_size_summary.txt'); end
    if ~okOwn, error('fig6_tracking_load:own', 'a panel B cut or CI differs from its run''s own result'); end
    if ~okX, error('fig6_tracking_load:x', 'a panel B mutation sd differs from its run or the expected size'); end
end


% ===================== drawing =======================================
function draw(R, B, outBase)
    figure_utils.setup_defaults();
    W = 8.7; H = 10.4;
    c.untouched = figure_utils.INK; c.tested = figure_utils.INK;   % ink throughout: no cooperation in this task
    fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', figure_utils.FONT, 'DefaultTextFontName', figure_utils.FONT);
    guard = onCleanup(@() close(fig)); %#ok<NASGU>
    axA = axes('Position', nr([1.7 5.85 6.65 3.9], W, H));     % panel A's box, in cm
    panel_load_cells(axA, R.C.L, c); figure_utils.head_cm(axA, 'A', 'Load per mutation', 1.51);
    axB = axes('Position', nr([1.7 1.0 6.65 3.3], W, H));
    panel_curve(axB, B, c); figure_utils.head_cm(axB, 'B', 'Target held still, responsiveness lowered', 1.51);
    ws = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector', 'Padding', 'figure');
    warning(ws);
    exportgraphics(fig, [outBase '.png'], 'Resolution', 300, 'Padding', 'figure');
    exportgraphics(fig, [outBase '.tif'], 'Resolution', 1200, 'Padding', 'figure');
    imwrite(imread([outBase '.tif']), [outBase '.tif'], 'Compression', 'lzw', 'Resolution', 1200);
    t = findall(fig, '-property', 'FontSize');
    fprintf('fig6_tracking_load: wrote %s.{pdf,png,tif} (%.1f x %.2f cm); font sizes %s\n', outBase, W, H, ...
        mat2str(unique(round(cell2mat(get(t, 'FontSize')), 2))'));
end


function panel_curve(ax, B, c)
% One network (design 1): seven points at mutation size 0.3, 0.5, 1, 2, 3, 5 and 10 times the standard.
    hold(ax, 'on');
    xl = [0.2 15];
    plot(ax, xl, [1 1], '--', 'Color', [0 0 0], 'LineWidth', 0.75, 'HandleVisibility', 'off');
    n = numel(B.x); hb = gobjects(1, n); hm = hb;
    for i = 1:n
        hb(i) = plot(ax, [B.x(i) B.x(i)], [B.ylo(i) B.yhi(i)], '-', 'Color', c.tested, 'LineWidth', 0.8, 'HandleVisibility', 'off');
        hm(i) = plot(ax, B.x(i), B.y(i), 'o', 'MarkerSize', 6, 'LineWidth', 0.9, 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', c.tested, 'HandleVisibility', 'off');   % panel A's "target held still" mark
    end
    ok = all(abs([hm.YData] - B.y) < 1e-6) && all(abs([hm.XData] - B.x) < 1e-12) ...
        && all(abs(cell2mat(get(hb(:), 'YData')) - [B.ylo(:) B.yhi(:)]) < 1e-6, 'all');
    fprintf('  CONTROL plotted B values = 1 - cut / 100 (1e-6, read back): %s\n', figure_utils.pass_fail(ok));
    if ~ok, error('fig6_tracking_load:control', 'a plotted panel B value differs from 1 - cut / 100'); end
    % legend: the control line only (a marker entry inside the data region would read as a fourth point;
    % 'responsiveness lowered' is in the caption, every point is that arm)
    h = plot(ax, NaN, NaN, '--', 'Color', [0 0 0], 'LineWidth', 0.75);
    set(ax, 'XScale', 'log', 'XLim', xl, 'XTick', [0.3 1 3 10], 'XTickLabel', {'0.3', '1', '3', '10'}, ...
        'YLim', [0 1.05], 'YTick', 0:0.25:1, 'YTickLabel', {'0', '0.25', '0.5', '0.75', '1'});
    xlabel(ax, ['mutation size (' char(215) ' standard, log scale)']); ylabel(ax, 'mutation load / control');
    lg = legend(ax, h, {'control level'}, 'Box', 'off', ...
        'FontSize', figure_utils.SIZE_SMALL, 'Location', 'southeast', 'TextColor', figure_utils.INK);
    lg.ItemTokenSize = [8 7];
    hold(ax, 'off'); figure_utils.style_panel(ax);
    set(ax, 'XMinorTick', 'off');
end


% ---- panel A ----
function panel_load_cells(ax, D, c)
% Panel A: each cell's measured against predicted load per mutation, with 95% intervals on both
% axes, the diagonal and a band of +-10%. Every open cell (target held still) is drawn first at a
% larger size and every filled cell (moving target) on top at a smaller one, so a ring stays
% visible around any overlap.
    hold(ax, 'on');
    xl = [0.006 0.45]; yl = xl;
    xs = logspace(log10(xl(1)), log10(xl(2)), 64);
    fill(ax, [xs fliplr(xs)], [0.9 * xs fliplr(1.1 * xs)], [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, xl, xl, '-', 'Color', figure_utils.DIAGONAL_GRAY, 'LineWidth', 0.75, 'HandleVisibility', 'off');
    tk = [0.01 0.03 0.1 0.3]; tl = {'0.01', '0.03', '0.1', '0.3'};
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xl, 'YLim', yl, 'XTick', tk, 'XTickLabel', tl, ...
        'YTick', tk, 'YTickLabel', tl);
    for k = 1:numel(D.x)
        errorbar(ax, D.x(k), D.L(k), D.L(k) - D.Llo(k), D.Lhi(k) - D.L(k), D.x(k) - D.xlo(k), D.xhi(k) - D.x(k), ...
            'LineStyle', 'none', 'Color', figure_utils.GRAY, 'LineWidth', 0.8, 'CapSize', 2, 'HandleVisibility', 'off');
    end
    for pass = [2 1]                   % 2 = held still (open, under), 1 = moving (filled, on top)
        for k = find(D.ai(:)' == pass)
            if pass == 1
                plot(ax, D.x(k), D.L(k), 'o', 'MarkerSize', 3.6, 'MarkerFaceColor', c.untouched, ...
                    'MarkerEdgeColor', c.untouched, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            else
                plot(ax, D.x(k), D.L(k), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', c.untouched, 'LineWidth', 0.9, 'HandleVisibility', 'off');
            end
        end
    end
    hM = plot(ax, NaN, NaN, 'o', 'MarkerSize', 3.6, 'MarkerFaceColor', c.untouched, 'MarkerEdgeColor', c.untouched);
    hF = plot(ax, NaN, NaN, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c.untouched, 'LineWidth', 0.9);
    lg = legend(ax, [hM hF], {'moving target', 'target held still'}, 'Box', 'off', 'FontSize', figure_utils.SIZE_SMALL, ...
        'Location', 'southeast', 'TextColor', figure_utils.INK);
    lg.ItemTokenSize = [8 7];
    xlabel(ax, 'predicted load per mutation'); ylabel(ax, 'measured load per mutation');
    hold(ax, 'off'); figure_utils.style_panel(ax);
    set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
end


% ---- end of panel A ----


function p = nr(r, W, H), p = [r(1) / W, r(2) / H, r(3) / W, r(4) / H]; end
