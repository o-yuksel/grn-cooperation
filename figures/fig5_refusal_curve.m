function fig5_refusal_curve(pubDir)
% FIG5_REFUSAL_CURVE  Draws Fig. 5: how often the continued saturating-contrast runs collapse
% (cooperation below 0.5) as a growing share of the mutations that end cooperation is refused.
%
% Produces : Fig. 5
% Inputs   : pubDir (optional; default the publication folder under the data root)
%            publication/refuse_less_conditional/report.mat   (per population the first continued
%              generation below 0.5, Inf if none, with and without mutation; run_mutation_off)
%            publication/refused_by_class/rruns_coopend_*.mat   (all refused; written by
%              run_all_refused)
%            publication/partial_refusal/rruns_<arm>_*.mat and partial_refusal_report.txt   (the five
%              partial-refusal arms, first batch; written by run_partial_refusal 'main', 'report')
%            publication/partial_refusal_batch2/rruns_<arm>_*.mat and partial_refusal_batch2_report.txt   (all
%              eight arms, second batch; written by run_partial_refusal 'main_batch2', 'report_batch2')
% Outputs  : output/figures/Fig5.pdf (vector, 8.7 x 6.0 cm), .png (300 dpi) and .tif (1200 dpi,
%            LZW), under the release root; prints the counts and x positions
% Run      : fig5_refusal_curve
% Seed     : none (no random draws)
%
% Each point pools 150 runs: populations 1 to 100, and 1 to 50 again with new seeds. A
% partial-refusal arm (`'refuse050'` to `'refuse098'`, a share q of 0.5 to 0.98 refused) sits at
% x = 100 (1 - S), S the share of the cooperation-ending mutations accepted with mutation as in the
% runs that the arm still accepts (pooled, partial_refusal_batch2_report.txt). Mutation as in the runs sits at 0,
% all refused at 100; mutation off is the grey line. Intervals: 95% Clopper-Pearson (binofit). The
% run stops unless every count and share equals the reports and the values set in load_points.
% Colors: the refusal arms remove cooperation-ending mutations, so they take the cooperation color
% (COLOR_COOP); mutation as in the runs is ink, mutation off mid grey.

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1, pubDir = fullfile(grn_data_file('publication')); end
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outBase = fullfile(outDir, 'Fig5');
    D = load_points(pubDir);
    draw(D, outBase);
    fprintf('fig5_refusal_curve: wrote %s.{pdf,png,tif}\n', outBase);
end


function D = load_points(pubDir)
% Pooled collapse counts (n = 150) and x positions for mutation as in the runs, the five
% partial-refusal arms, all refused and mutation off.
    S = load(fullfile(pubDir, 'refuse_less_conditional', 'report.mat')); T = S.T;
    arms = {'refuse050', 'refuse082', 'refuse090', 'refuse095', 'refuse098'};
    first = [sum(isfinite(T.control.cg)), zeros(1, 5), count_arm(fullfile(pubDir, 'refused_by_class'), 'coopend', 100), sum(isfinite(T.off.cg))];
    rep1 = fileread(fullfile(pubDir, 'partial_refusal', 'partial_refusal_report.txt'));
    for a = 1:5
        first(1 + a) = count_arm(fullfile(pubDir, 'partial_refusal'), arms{a}, 100);
        c = regexp(rep1, ['\n\s+' either(arms{a}) '\s+collapses\s+(\d+)'], 'tokens', 'once');
        if str2double(c{1}) ~= first(1 + a), error('fig5:report', '%s: report COUNTS %s, runs %d', arms{a}, c{1}, first(1 + a)); end
    end
    sdir = fullfile(pubDir, 'partial_refusal_batch2');
    names = [{'on'}, arms, {'all_refused', 'off'}];
    second = cellfun(@(a) count_arm(sdir, a, 50), names);
    rep2 = fileread(fullfile(sdir, 'partial_refusal_batch2_report.txt'));
    for i = 1:8
        c = regexp(rep2, ['\n' either(names{i}) '\s+\S+\s+\|\s+(\d+) '], 'tokens', 'once');
        if str2double(c{1}) ~= second(i), error('fig5:report_batch2', '%s: batch-2 report %s, runs %d', names{i}, c{1}, second(i)); end
    end
    x = cellfun(@(a) pooled_share(rep2, a), arms);
    pooled = first + second;
    if ~isequal(first, [36 33 27 17 21 10 7 7]) || ~isequal(second, [24 24 16 11 13 11 6 4]) ...
            || ~isequal(pooled, [60 57 43 28 34 21 13 11]) || ~isequal(round(x * 1e4) / 1e4, [0.5442 0.2190 0.1212 0.0621 0.0227])
        error('fig5:counts', 'collapse counts or x values differ from the paper''s values');
    end
    D.N = 150; D.x = 100 * (1 - [1 x 0]); D.n = pooled(1:7); D.off = pooled(8);
    fprintf('fig5_refusal_curve: first batch %s, second batch %s\n', mat2str(first), mat2str(second));
    fprintf('fig5_refusal_curve: pooled (of 150) on %d, partial refusal %s, all refused %d, off %d\n', pooled(1), mat2str(pooled(2:6)), pooled(7), pooled(8));
    fprintf('fig5_refusal_curve: x (percent removed) = %s\n', mat2str(D.x, 4));
end


function c = count_arm(sdir, arm, nWant)
% Collapses (any hc < 0.5) among populations 1..nWant of one arm; checks the labels.
    f = dir(fullfile(sdir, ['rruns_' arm '_*.mat'])); R = {}; lab = {};
    for i = 1:numel(f), S = load(fullfile(sdir, f(i).name)); R = [R; S.R(:)]; lab = [lab; S.lab]; end %#ok<AGROW>
    k = cell2mat(lab(:, 1));
    if numel(R) ~= nWant || ~isequal(sort(k(:))', 1:nWant) || ~all(strcmp(stored_value('arm', lab(:, 2)), arm))
        error('fig5:arm', '%s: expected populations 1-%d with that arm label', arm, nWant);
    end
    c = sum(cellfun(@(r) any(r.hc < 0.5), R));
end


function r = either(arm)
% The arm's name in a report: its current name or the spelling of the stored text (stored_value),
% as a regular expression.
    r = ['(?:' arm '|' stored_value('arm', arm, 'stored') ')'];
end


function s = pooled_share(rep, arm)
% The arm's pooled share S of the cooperation-ending mutations accepted with mutation as in the
% runs, from its line in the [SHARE LEFT] block of partial_refusal_batch2_report.txt.
    t = regexp(rep, ['\n\s+' either(arm) '\s+new 50[^\n]*pooled \d+ / \d+ = ([\d.]+)'], 'tokens', 'once');
    if isempty(t), error('fig5:report', 'no pooled share for %s in partial_refusal_batch2_report.txt', arm); end
    s = str2double(t{1});
end


function draw(D, outBase)
    figure_utils.setup_defaults();
    W = 8.7; H = 6.0; FS = 7; INK = figure_utils.INK; BLUE = figure_utils.COLOR_COOP; GRAY = figure_utils.GRAY;
    BARBLUE = 1 - 0.4 * (1 - BLUE);                        % pale tint of BLUE for the interval bars of the refusal arms
    BARINK = 1 - 0.4 * (1 - INK);                          % same tint of INK for the mutation-on bar
    fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', figure_utils.FONT, 'DefaultTextFontName', figure_utils.FONT);
    guard = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Position', [1.3 1.15 7.15 4.65] ./ [W H W H]); hold(ax, 'on');
    xl = [-3 103]; yl = [0 50];
    pct = D.n * 100 / D.N; offPct = D.off * 100 / D.N;
    rot = -atan2d(31.33 * 4.65 / diff(yl), 100 * 7.15 / diff(xl));   % slope of the dashed line on the page (31.33 = 40.0 - 8.67 percent; places the label)
    plot(ax, xl, [offPct offPct], '-', 'Color', GRAY, 'LineWidth', 0.75);
    plot(ax, [0 100], [pct(1) pct(end)], '--', 'Color', INK, 'LineWidth', 0.8);
    [~, ci] = binofit(D.n, D.N, 0.05); ci = ci * 100;
    for i = 1:numel(D.x)
        col = BARBLUE; if i == 1, col = BARINK; end
        plot(ax, D.x([i i]), ci(i, :), '-', 'Color', col, 'LineWidth', 1.6);
    end
    plot(ax, D.x(1), pct(1), 'o', 'MarkerSize', 4.5, 'MarkerFaceColor', INK, 'MarkerEdgeColor', INK, 'LineWidth', 0.75);
    plot(ax, D.x(2:end), pct(2:end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'LineWidth', 0.75);
    text(ax, 56, 20.4, 'in proportion', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', FS, 'Color', INK, 'Rotation', rot);
    text(ax, 2, 4.6, 'no mutation', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', FS, 'Color', INK);
    set(ax, 'XLim', xl, 'YLim', yl, 'XTick', 0:25:100, 'YTick', 0:10:50, 'FontName', figure_utils.FONT, 'FontSize', FS, ...
        'LineWidth', 0.6, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off', 'XColor', INK, 'YColor', INK, ...
        'LabelFontSizeMultiplier', 1, 'Layer', 'top');
    grid(ax, 'off');
    xlabel(ax, 'cooperation-ending mutations removed (%)'); ylabel(ax, 'runs that collapsed (%)');
    set([ax.XLabel ax.YLabel], 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK, 'FontName', figure_utils.FONT);
    export_all(fig, outBase);
end


function export_all(fig, base)
    ws = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', 'Padding', 'figure');
    warning(ws);
    exportgraphics(fig, [base '.png'], 'Resolution', 300, 'Padding', 'figure');
    exportgraphics(fig, [base '.tif'], 'Resolution', 1200, 'Padding', 'figure');
    imwrite(imread([base '.tif']), [base '.tif'], 'Compression', 'lzw', 'Resolution', 1200);
end
