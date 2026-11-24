function R = figS1_draw(cfg)
% FIGS1_DRAW  Draws and checks Fig. S1: in 65 populations about to fail, the share of runs that
% collapse within 100 generations as observed, under the sham repair and under the repair, with
% and without the defector added (A); and the cumulative collapse with the defector added (B).
%
% Implements : Fig. S1 (called by figS1_repair_test); reprints as checks the levels and
%              differences of SI S4.3 and Table S10
% Inputs     : cfg struct: pubDir (default the publication folder under the data root), figDir
%              (default output/figures)
%              the per-population shares and the population set, from
%                repair_test_summary('main') (which rewrites repair_matched.txt)
%              publication/repair_matched/repair_matched_build_targets.mat and repair_matched_main.mat
%                (written by run_repair_test; W99s: targets whose cooperation is at least 0.99)
%              publication/repair_test_observed/repair_test_observed.mat   (written by
%                run_repair_test_observed)
% Outputs    : R (levels, differences); [figDir '/figS1_draw'] .pdf (17.8 x 5.4 cm), .png
%              (300 dpi) and .tif (1200 dpi, LZW); prints the values and checks
% Seed       : rng(20261013, 'twister'), then randi(n, n, 10000) over the n populations: the
%              summary's resamples, so the intervals equal the report's (the caller's stream is
%              restored afterwards)
%
% Arms: observed (stored arm `A`; `A0` without the defector), sham repair (`Sp`, `S0p`), repair
% (`Bp`, `B0p`). A bar is the mean over populations of each one's share of runs that collapse
% (cooperation below 0.5 in some generation); whiskers, 95% percentile bootstrap intervals. The run
% stops unless there are 65 populations, the paired differences reproduce the summary's
% (repair_test_summary), and the time course at generation 100 equals the bars; the levels and
% differences are also printed against the values the paper states (PASS or FAIL, report only).

    if nargin < 1, cfg = struct(); end
    thisDir = fileparts(mfilename('fullpath'));
    if ~isfield(cfg, 'pubDir'), cfg.pubDir = grn_data_file('publication'); end
    if ~isfield(cfg, 'figDir'), cfg.figDir = fullfile(fileparts(thisDir), 'output', 'figures'); end
    S = repair_test_summary('main');          % the summary's own per-event shares
    X = S.X;
    rs = rng; rng(20261013, 'twister');            % the summary's draws
    bi = randi(numel(X.pA), numel(X.pA), 10000);
    rng(rs);
    R = levels(X, bi);
    R.D = differences(X, bi);
    check(R);
    check_levels(R);
    T = trajectories(cfg.pubDir, X, R);
    draw_combined(R, T, cfg);
end


function R = levels(X, bi)
% Rows: defector added, no defector added. Columns: observed, sham repair, repair.
    v = {X.pA, X.pSp, X.pBp; X.pA0, X.pS0p, X.pB0p};
    R.n = numel(X.pA); R.mu = nan(2, 3); R.lo = R.mu; R.hi = R.mu;
    for g = 1:2
        for a = 1:3
            [R.mu(g, a), ci] = boot(v{g, a}, bi); R.lo(g, a) = ci(1); R.hi(g, a) = ci(2);
        end
    end
end


function D = differences(X, bi)
% Row 1 sham - repair, row 2 observed - sham; columns defector added, no defector added.
    d = {X.pSp - X.pBp, X.pS0p - X.pB0p; X.pA - X.pSp, X.pA0 - X.pS0p};
    D.mu = nan(2, 2); D.lo = D.mu; D.hi = D.mu;
    for r = 1:2
        for g = 1:2
            [D.mu(r, g), ci] = boot(d{r, g}, bi); D.lo(r, g) = ci(1); D.hi(r, g) = ci(2);
        end
    end
end


function [m, ci] = boot(x, bi)
    m = mean(x); ci = quantile(mean(x(bi), 1), [0.025 0.975]);
end


function check(R)
    D = R.D;
    got = [D.mu(1, 1) D.lo(1, 1) D.hi(1, 1); D.mu(2, 1) D.lo(2, 1) D.hi(2, 1); D.mu(2, 2) D.lo(2, 2) D.hi(2, 2)];
    ref = [0.266 0.161 0.372; 0.089 0.015 0.168; 0.128 0.045 0.215];
    if R.n ~= 65, error('figS1_draw:control', '%d populations, expected 65', R.n); end
    ok = all(abs(round(got, 3) - ref) < 1e-9, 'all');
    fprintf('REPAIR (matched start), n = %d populations; observed / sham / repair\n', R.n);
    nm = {'defector added', 'no defector added'};
    for g = 1:2
        fprintf('  %-13s %s\n', nm{g}, sprintf('%.3f [%.3f, %.3f]  ', [R.mu(g, :); R.lo(g, :); R.hi(g, :)]));
        fprintf('  %-13s sham - repair %+.3f [%+.3f, %+.3f] | observed - sham %+.3f [%+.3f, %+.3f]\n', '', ...
            D.mu(1, g), D.lo(1, g), D.hi(1, g), D.mu(2, g), D.lo(2, g), D.hi(2, g));
    end
    fprintf('  CONTROL S''-B'' +0.266 [+0.161, +0.372], A-S'' +0.089 [+0.015, +0.168], A0-S0'' +0.128 [+0.045, +0.215] (the paper''s): %s\n', figure_utils.pass_fail(ok));
end


function check_levels(R)
% The drawn bars and whisker ends against the six intervals written in below (columns mu, lo,
% hi; rows observed, sham, repair; defector added, then no defector added).
    ref = {[0.840 0.774 0.902; 0.751 0.665 0.832; 0.485 0.378 0.595], ...
           [0.732 0.655 0.806; 0.605 0.506 0.702; 0.263 0.182 0.355]};
    ok = true;
    for g = 1:2
        got = [R.mu(g, :); R.lo(g, :); R.hi(g, :)]';
        ok = ok && all(abs(round(got, 3) - ref{g}) < 1e-9, 'all');
    end
    fprintf('  CONTROL drawn levels and whisker ends against the paper''s, six arms to 3 decimals: %s\n', figure_utils.pass_fail(ok));
end


function T = trajectories(pub, X, R)
% Cumulative collapse (hc below 0.5 at any generation up to g) of A, S', B', defector added:
% mean over populations of each population's share of its runs.
    dir0 = fullfile(pub, 'repair_matched');
    Bd = load(fullfile(dir0, 'repair_matched_build_targets.mat'), 'out'); t = Bd.out.targets;
    t = t([t.inRun] & ~[t.excluded]);
    assert(isequal({t.id}', X.id), 'event order differs from the summary''s');
    M = load(fullfile(dir0, 'repair_matched_main.mat'), 'out'); M = M.out.runs;
    I = load(fullfile(pub, 'repair_test_observed', 'repair_test_observed.mat'), 'out'); I = I.out.runs;
    pools = {I, M, M}; arms = {'A', 'Sp', 'Bp'}; nG = 100;
    T.g = 1:nG; T.mu = nan(3, nG);
    for a = 1:3
        V = nan(numel(t), nG);
        for i = 1:numel(t)
            hc = double(pick(pools{a}, t(i).num, arms{a}).hc);
            if size(hc, 2) ~= nG, error('figS1_draw:trace', 'hc has %d generations, expected %d', size(hc, 2), nG); end
            V(i, :) = mean(cummax(double(hc < 0.5), 2), 1);
        end
        T.mu(a, :) = mean(V, 1);
    end
    dev = max(abs(T.mu(:, end)' - R.mu(1, :)));
    fprintf('  Fig. S1B at g = 1, 25, 100: %s\n', sprintf('[%.3f %.3f %.3f] ', T.mu(:, [1 25 100])'));
    fprintf('  CONTROL cumulative collapse at g = 100 vs endpoint means, max |diff| %.4f: %s\n', dev, figure_utils.pass_fail(dev < 0.005));
    if dev >= 0.005, error('figS1_draw:control', 'cumulative collapse at g = 100 does not match the endpoint means'); end
end


function r = pick(pool, num, arm)
    j = find([pool.event] == num & strcmp({pool.arm}, arm));
    if numel(j) ~= 1, error('figS1_draw:trace', 'event %d arm %s: %d records', num, arm, numel(j)); end
    r = pool(j).r;
end


function [col, mk] = arm_style()
% Observed in ink (square), sham repair mid grey (triangle), repair in the defection color
% (circle): the repair restores retaliation, defecting back, which COLOR_DEFECT means.
    col = [figure_utils.INK; figure_utils.GRAY; figure_utils.COLOR_DEFECT]; mk = {'s', '^', 'o'};
end


function draw_bars(R, cfg)
% Fig. S1A, drawn into the shared figure of draw_combined (cfg.fig, cfg.W, cfg.H): per group the
% three arms' bars with their whiskers, each bar named under itself.
    figure_utils.setup_defaults();
    W = cfg.W; H = cfg.H; fig = cfg.fig;
    [col, ~] = arm_style();
    axA = axes(fig, 'Position', [1.35 0.9 7.0 3.9] ./ [W H W H]); hold(axA, 'on');
    off = [-0.27 0 0.27];
    for g = 1:2
        for a = 1:3
            x = g + off(a);
            bar(axA, x, R.mu(g, a), 0.25, 'FaceColor', col(a, :), 'EdgeColor', 'none', 'LineWidth', 0.75);
            plot(axA, [x x], [R.lo(g, a) R.hi(g, a)], '-', 'Color', figure_utils.INK, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            if a == 1   % over the dark observed bar the whisker is white, ink where it leaves the bar
                plot(axA, [x x], [R.lo(g, a) min(R.hi(g, a), R.mu(g, a))], '-', 'Color', 'w', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        end
    end
    set(axA, 'XLim', [0.5 2.5], 'XTick', 1:2, 'XTickLabel', {});
    style3(axA); set(axA, 'YLim', [0 1.2]);
    axA.YLabel.Units = 'normalized'; axA.YLabel.Position(1) = -0.14;   % clear of the tick labels
    grpNames = {{'defector added', sprintf('%d populations', R.n)}, {'no defector added', sprintf('%d populations', R.n)}};
    % Each bar named under itself, so color is not the only cue; the group names sit above their
    % clusters and there is no color legend
    armNames = {'observed', {'sham', 'repair'}, 'repair'};
    for g = 1:2
        for a = 1:3
            text(axA, g + off(a), -0.03 * diff(axA.YLim), armNames{a}, 'FontSize', figure_utils.SIZE_SMALL, 'Color', figure_utils.INK, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Clipping', 'off');
        end
        text(axA, g, axA.YLim(2), grpNames{g}, 'FontSize', figure_utils.SIZE_SMALL, 'Color', figure_utils.INK, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
    report_overlap(axA);
    hold(axA, 'off');
    set(axA.YLabel, 'FontSize', figure_utils.SIZE_LABEL);   % axis title 8 pt, as panel B
    head3(axA, 'A', 'Collapse within 100 generations', 1.16, 8);
end


function report_overlap(ax)
% Control: no two text labels of the panel overlap, and no label reaches into a bar or whisker
% (extents in data units; an error if any does).
    t = findall(ax, 'Type', 'text', 'Units', 'data');
    e = cell2mat(arrayfun(@(h) h.Extent, t, 'UniformOutput', false));   % [x y w h]
    hit = 0;
    for i = 1:numel(t)
        for j = i + 1:numel(t)
            if e(i, 1) < e(j, 1) + e(j, 3) && e(j, 1) < e(i, 1) + e(i, 3) && e(i, 2) < e(j, 2) + e(j, 4) && e(j, 2) < e(i, 2) + e(i, 4)
                hit = hit + 1; fprintf('  overlap: "%s" / "%s"\n', strjoin(cellstr(t(i).String), ' '), strjoin(cellstr(t(j).String), ' '));
            end
        end
    end
    l = findall(ax, 'Type', 'line'); top = max(arrayfun(@(h) max(h.YData), l));
    low = min(e(e(:, 2) > 0, 2));   % the group names above the clusters
    fprintf('  CONTROL panel A labels: %d overlaps; group names start at %.3f, highest whisker %.3f: %s\n', ...
        hit, low, top, figure_utils.pass_fail(hit == 0 && low > top));
    if hit > 0 || low <= top, error('figS1_draw:control', 'panel A labels overlap'); end
end


function draw_time_course(T, cfg)
% Fig. S1B, drawn into the shared figure of draw_combined (cfg.fig, cfg.W, cfg.H, cfg.x0):
% cumulative collapse over generations, plain lines with one end marker.
    figure_utils.setup_defaults();
    [col, mk] = arm_style();
    fig = cfg.fig; W = cfg.W; H = cfg.H; pos = [cfg.x0 + 1.35, 0.9, 7.0, 3.9];
    ax = axes(fig, 'Position', pos ./ [W H W H]); hold(ax, 'on');
    for a = 1:3
        plot(ax, T.g, T.mu(a, :), '-', 'Color', col(a, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(ax, T.g(end), T.mu(a, end), mk{a}, 'MarkerSize', 5, 'LineWidth', 0.8, 'MarkerFaceColor', col(a, :), ...
            'MarkerEdgeColor', col(a, :), 'Clipping', 'off', 'HandleVisibility', 'off');
        hk(a) = plot(ax, NaN, NaN, '-', 'Color', col(a, :), 'LineWidth', 1.5, 'Marker', mk{a}, ...
            'MarkerSize', 5, 'MarkerFaceColor', col(a, :), 'MarkerEdgeColor', col(a, :)); %#ok<AGROW>
    end
    set(ax, 'XLim', [0 100], 'XTick', 0:25:100);
    xlabel(ax, 'generations after restart');
    lg = legend(ax, hk, {'observed', 'sham repair', 'repair'}, 'Box', 'off', 'FontSize', figure_utils.SIZE_SMALL, ...
        'Location', 'southeast', 'TextColor', figure_utils.INK);
    lg.ItemTokenSize = [14 7]; hold(ax, 'off'); style3(ax);
    head3(ax, 'B', 'Cumulative collapse after the restart', 1.16, 8);
end


function draw_combined(R, T, cfg)
% Fig. S1: bars (A) and time course (B) side by side in one 17.8 x 5.4 cm figure.
    figure_utils.setup_defaults();
    W = 17.8; H = 5.4;
    fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', figure_utils.FONT, 'DefaultTextFontName', figure_utils.FONT);
    guard = onCleanup(@() close(fig)); %#ok<NASGU>
    c = struct('fig', fig, 'W', W, 'H', H);
    draw_bars(R, c);
    c.x0 = 8.9; draw_time_course(T, c);
    export_cm(fig, fullfile(cfg.figDir, 'figS1_draw'), W, H);
    figure_utils.check_pdf_width(fullfile(cfg.figDir, 'figS1_draw.pdf'), W);
end


function export_cm(fig, base, W, H)
    ws = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', 'Padding', 'figure');
    warning(ws);
    exportgraphics(fig, [base '.png'], 'Resolution', 300, 'Padding', 'figure');
    exportgraphics(fig, [base '.tif'], 'Resolution', 1200, 'Padding', 'figure');
    imwrite(imread([base '.tif']), [base '.tif'], 'Compression', 'lzw', 'Resolution', 1200);
    t = findall(fig, '-property', 'FontSize');
    fprintf('wrote %s.{pdf,png,tif} (%.1f x %.1f cm); font sizes %s\n', base, W, H, ...
        mat2str(unique(round(cell2mat(get(t, 'FontSize')), 2))'));
end


function style3(ax)
% One type and line style for both panels, no grid.
    set(ax, 'YLim', [0 1], 'YTick', 0:0.25:1, 'FontName', figure_utils.FONT, 'FontSize', figure_utils.SIZE_SMALL, 'LineWidth', 0.6, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off', 'XColor', figure_utils.INK, ...
        'YColor', figure_utils.INK, 'XTickLabelRotation', 0, 'LabelFontSizeMultiplier', 1);
    grid(ax, 'off');
    ylabel(ax, 'fraction collapsed');
    set([ax.XLabel ax.YLabel], 'FontSize', figure_utils.SIZE_LABEL, 'Color', figure_utils.INK, 'FontName', figure_utils.FONT);
end


function head3(ax, letter, ttl, dxCm, fs)
% Letter 10 pt bold, title 8 pt, one baseline; the letter sits dxCm left of the plot box.
    if nargin < 5, fs = 8; end
    u = ax.Units; ax.Units = 'centimeters'; wcm = ax.Position(3); ax.Units = u;
    if ~isempty(letter), text(ax, -dxCm / wcm, 1.04, letter, 'Units', 'normalized', 'FontSize', figure_utils.SIZE_LETTER, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom'); end
    text(ax, 0, 1.04, ttl, 'Units', 'normalized', 'FontSize', fs, 'Color', figure_utils.INK, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
