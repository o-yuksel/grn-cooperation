function fig = fig2_draw(cfg)
% FIG2_DRAW  Draws and checks Fig. 2: per signaling architecture, mean cooperation over
% generations 1 to 500 with its 95% band and, above it, the direction of selection on what
% offspring inherit (A); the histogram of generation-500 cooperation with the number of the 100
% populations that succeed (B).
%
% Implements : Fig. 2 (called by fig2_signaling_architecture)
% Inputs     : cfg   figure_config struct (default figure_config()); cfg.outStem (optional),
%                    the output path without extension (default output/figures/fig2_draw)
%              the curves and histograms, from fig2_data
%              publication/selection_over_time/selection_z_<arch>_r<NN>.mat, selection_sign_<arch>_r<NN>.mat
%                and lineage_records_<arch>_r<NN>.mat, replicates 1 to 20   (written by run_selection_z,
%                run_selection_sign and run_lineage_records)
%              publication/selection_over_time/summary_selection_inherited.txt   (written by
%                selection_inherited; read from the same publication folder)
% Outputs    : fig; [outStem '.tif'] (600 dpi, LZW), '.pdf' and '.png' (150 dpi), through
%              figure_utils.export; prints the checks
% Seed       : rng(20261022), default generator, before the bootstrap (the caller's stream is
%              restored afterwards)
%
% S_G = Cov(n, g) in each generation, n an individual's offspring count and g the mean cooperation
% of the individuals sharing its wiring (as in selection_inherited). A strip block reads "dilemma
% escaped" when the 95% bootstrap interval of the 20-run mean (10,000 resamples of runs) lies above
% zero, "dilemma holds" when below, else "no detectable selection". The strip takes the two fills
% figure_utils.DILEMMA_HOLDS_TINT and DILEMMA_ESCAPED_TINT, the filled cells of Fig. 1C, and leaves
% "no detectable selection" white. The run stops unless the layout checks pass, the strip equals
% the summary, and the drawn counts and curves equal fig2_data's.

    if nargin < 1
        cfg = figure_config();
        figure_utils.setup_defaults();
    end
    D = fig2_data(cfg);
    G = load_selection_background(cfg, D);
    check_summary(G, cfg);

    L = layout_in();
    fig = figure_utils.create_figure(L.W, L.H);
    ov = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', 'HitTest', 'off', 'Tag', 'legend_overlay');
    xlim(ov, [0 L.W]);  ylim(ov, [0 L.H]);  hold(ov, 'on');

    rowsA = draw_rows_A(fig, D, G, L);
    check_patches(rowsA, G);
    rowsB = draw_strips_B(fig, D, L);
    draw_annotations(ov, D, G, L);

    check_curves_unchanged(rowsA, D);
    check_alignment(rowsA, rowsB);
    check_counts(ov, D);
    check_clearances(ov, rowsB, D);
    check_polish(ov, rowsA, rowsB, D);
    check_strip_strength(rowsA, rowsB, D);
    check_encoding(rowsA, rowsB, ov, D);
    check_text_and_lines(fig);

    if isfield(cfg, 'outStem') && ~isempty(cfg.outStem)
        base = cfg.outStem;
    else
        outDir = fullfile(fileparts(cfg.analysisDir), 'output', 'figures');
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        base = fullfile(outDir, 'fig2_draw');
    end
    [~, stemName] = fileparts(base);
    spec = struct('nPanels', 2, 'panelLabels', {{'A', 'B'}}, ...
        'yMaxExpected', 100, 'name', stemName);
    figure_utils.export(fig, base, 600, spec);   % TIFF at 600 dpi
    fprintf('fig2_draw exported: %s\n', base);
end

function L = layout_in()
%LAYOUT_IN  All positions in inches on the canvas (origin bottom-left).
    L.W = 7.08;  L.rh = 0.50;  L.gap = 0.16;  L.bot = 0.56;  L.top = 0.52;
    L.H = L.bot + 4 * L.rh + 3 * L.gap + L.top;
    L.xA = 1.00;  L.wA = 2.56;  L.xB = 3.94;  L.wB = 2.56;
    L.fs = 7;                                  % one text size for the figure
    L.histAlpha   = 0.55;                      % panel-B histogram FaceAlpha only
    L.tintHolds   = figure_utils.DILEMMA_HOLDS_TINT;     % strip and legend swatches, opaque; Fig. 1C's "dilemma holds" fill
    L.tintEscaped = figure_utils.DILEMMA_ESCAPED_TINT;   % Fig. 1C's "dilemma escaped" fill
    L.tintNone    = [1 1 1];                   % "no detectable selection": white inside the strip outline
    L.arch  = figure_utils.architecture_palette().colors;
    L.strip = 0.075;  L.sgap = 0.012;
    L.lw = 1.5;                                % curve and row-swatch width (pt)
    L.tickLen = 0.015;                         % tick length (fraction of axis length)
    L.dark = 0.72;
    L.bandCol = L.dark * L.arch;               % band color = curve color
    L.bandAlpha = 0.30;
    L.base = figure_utils.BLACK;
    L.axisInk = figure_utils.INK;              % axes, ticks and the y-axis line
end

function p = pos_in(L, x, y, w, h)
    p = [x / L.W, y / L.H, w / L.W, h / L.H];
end

function yb = row_bottom(L, m)
%ROW_BOTTOM  Bottom (inches) of row m (m = 1 graded separate, bottom row).
    yb = L.bot + (m - 1) * (L.rh + L.gap);
end

function C = state_colors(L)
%STATE_COLORS  Rows: dilemma holds (-1), no detectable selection (0, white), dilemma escaped (+1); constants in layout_in.
    C = [L.tintHolds; L.tintNone; L.tintEscaped];
end

function rows = draw_rows_A(fig, D, G, L)
    rows = struct('ax', cell(D.nModes, 1), 'strip', cell(D.nModes, 1));
    for m = 1:D.nModes
        C = state_colors(L);
        col = L.dark * L.arch(m, :);
        ax = axes(fig, 'Position', pos_in(L, L.xA, row_bottom(L, m), L.wA, L.rh));
        hold(ax, 'on');
        if L.strip > 0   % state strip above the plot area; the plot background stays white
            sax = axes(fig, 'Position', pos_in(L, L.xA, row_bottom(L, m) + L.rh + L.sgap, L.wA, L.strip));
            hold(sax, 'on');  set(sax, 'XScale', 'log', 'Visible', 'off');  xlim(sax, [1 500]);  ylim(sax, [0 100]);
            draw_state_patches(sax, G.state(:, m), C);
            patch(sax, [1 500 500 1], [0 0 100 100], 'w', 'FaceColor', 'none', 'EdgeColor', figure_utils.GRAY, 'LineWidth', 0.3, ...
                'Clipping', 'off', 'HandleVisibility', 'off', 'Tag', 'strip_outline');
            rows(m).strip = sax;
        else
            draw_state_patches(ax, G.state(:, m), C);
        end
        fill(ax, [D.gens, fliplr(D.gens)], [D.coopCiLo(m, :), fliplr(D.coopCiHi(m, :))], ...
            L.bandCol(m, :), 'FaceAlpha', L.bandAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off', 'Tag', 'curve_band');
        plot(ax, D.gens, D.meanCoopOverTime(m, :), '-', 'Color', col, 'LineWidth', L.lw, ...
            'HandleVisibility', 'off', 'Tag', 'curve_line');
        set(ax, 'XScale', 'log', 'Layer', 'top', 'Box', 'off', 'TickDir', 'out', 'TickLength', [L.tickLen L.tickLen], ...
            'FontName', figure_utils.FONT, 'FontSize', L.fs, 'LineWidth', 0.6, ...
            'XColor', L.axisInk, 'YColor', L.axisInk, 'Color', 'none');
        xlim(ax, [1 500]);  ylim(ax, [0 100]);
        yticks(ax, [0 100]);  yticklabels(ax, {});   % labels drawn in draw_annotations
        xticks(ax, [1 10 100 500]);
        ax.XMinorTick = ternary_(m == 1, 'on', 'off');   % full log ticks on the bottom row only
        plot(ax, [1 1], [0 100], '-', 'Color', L.axisInk, 'LineWidth', 0.6, 'Clipping', 'off', ...
            'HandleVisibility', 'off', 'Tag', 'yaxis_line');
        if m == 1
            xticklabels(ax, {'1', '10', '100', '500'});
            xlabel(ax, 'generation', 'FontSize', figure_utils.SIZE_LABEL, 'FontName', figure_utils.FONT, 'FontWeight', 'normal');
        else
            xticklabels(ax, {});
        end
        rows(m).ax = ax;
    end
end

function draw_state_patches(ax, st, C)
%DRAW_STATE_PATCHES  One patch per run of equal state; edges at geometric midpoints, 1 and 500 at the ends.
    gs = (1:numel(st))';
    e = [1; sqrt(gs(1:end-1) .* gs(2:end)); numel(st)];
    i = 1;
    while i <= numel(st)
        j = i;
        while j < numel(st) && st(j + 1) == st(i), j = j + 1; end
        c = C(st(i) + 2, :);
        patch(ax, [e(i) e(j+1) e(j+1) e(i)], [0 0 100 100], c, 'EdgeColor', c, 'LineWidth', 0.25, ...
            'HandleVisibility', 'off', 'Tag', 'bg_patch', 'UserData', [i j st(i)]);
        i = j + 1;
    end
end

function rows = draw_strips_B(fig, D, L)
%DRAW_STRIPS_B  One histogram strip per architecture on the same rows as panel A (common count scale).
    nBins = 20;  edges = linspace(0, 1, nBins + 1);
    cnt = cell(D.nModes, 1);  cmax = 0;
    for m = 1:D.nModes
        cnt{m} = histcounts(D.finalCoopByMode{m}, edges);
        cmax = max(cmax, max(cnt{m}));
    end
    rows = struct('ax', cell(D.nModes, 1));
    for m = 1:D.nModes
        ax = axes(fig, 'Position', pos_in(L, L.xB, row_bottom(L, m), L.wB, L.rh));
        hold(ax, 'on');
        yTopB = 1.15 * cmax;   % band and axes above the tallest bar
        fill(ax, [0.95 1 1 0.95], [0 0 yTopB yTopB], [0.85 0.85 0.85], 'EdgeColor', 'none', ...
            'HandleVisibility', 'off', 'Tag', 'success_shade');
        xp = [edges(1), reshape([edges(1:end-1); edges(2:end)], 1, []), edges(end)];
        yp = [0, reshape([cnt{m}; cnt{m}], 1, []), 0];
        fill(ax, xp, yp, L.arch(m, :), 'FaceColor', L.arch(m, :), 'FaceAlpha', L.histAlpha, 'EdgeColor', L.arch(m, :), 'LineWidth', 0.9, ...
            'HandleVisibility', 'off', 'Tag', 'hist_bars');
        plot(ax, [0.95 0.95], [0 yTopB], ':', 'Color', [0.40 0.40 0.40], 'LineWidth', 0.75, ...
            'HandleVisibility', 'off', 'Tag', 'success_line');
        set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [L.tickLen L.tickLen], 'Layer', 'top', ...
            'FontName', figure_utils.FONT, 'FontSize', L.fs, 'LineWidth', 0.6, 'XColor', L.axisInk, 'YColor', 'none', ...
            'YTick', [], 'Color', 'none');
        xlim(ax, [0 1]);  ylim(ax, [0 yTopB]);
        xticks(ax, [0 0.5 1]);  ax.XMinorTick = 'off';
        if m == 1
            xticklabels(ax, {'0', '0.5', '1'});
            xlabel(ax, 'cooperation at generation 500', 'FontSize', figure_utils.SIZE_LABEL, 'FontName', figure_utils.FONT, 'FontWeight', 'normal');
        else
            xticklabels(ax, {});
        end
        if m == D.nModes
            barLen = 20;  y0 = 0.30 * cmax;  x0 = 0.12;  cap = 0.012;
            plot(ax, [x0 x0], [y0 y0 + barLen], '-', 'Color', L.base, 'LineWidth', 0.75, 'HandleVisibility', 'off');
            plot(ax, x0 + cap * [-1 1], [y0 y0], '-', 'Color', L.base, 'LineWidth', 0.75, 'HandleVisibility', 'off');
            plot(ax, x0 + cap * [-1 1], [y0 y0] + barLen, '-', 'Color', L.base, 'LineWidth', 0.75, 'HandleVisibility', 'off');
            text(ax, x0 + 0.04, y0 + barLen / 2, '20 populations', 'FontSize', L.fs, 'FontName', figure_utils.FONT, ...
                'Color', L.base, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
        rows(m).ax = ax;
    end
end

function draw_annotations(ov, D, G, L)
%DRAW_ANNOTATIONS  Letters, shared y label, row names, success counts and label, legend (overlay axes in inches).
    yTop = row_bottom(L, D.nModes) + L.rh;
    yHead = yTop + 0.20;
    txt = @(x, y, s, varargin) text(ov, x, y, s, 'FontSize', L.fs, 'FontName', figure_utils.FONT, 'Color', L.base, varargin{:});
    text(ov, 0.02, yHead, 'A', 'FontSize', figure_utils.SIZE_LETTER, 'FontWeight', 'bold', 'FontName', figure_utils.FONT, 'Color', L.base, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Tag', 'letterA');
    text(ov, L.xB, yHead, 'B', 'FontSize', figure_utils.SIZE_LETTER, 'FontWeight', 'bold', 'FontName', figure_utils.FONT, 'Color', L.base, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Tag', 'letterB');
    yMid = row_bottom(L, 1) + (4 * L.rh + 3 * L.gap) / 2;
    txt(0.675, yMid, 'cooperation (%)', 'FontSize', figure_utils.SIZE_LABEL, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    for m = 1:D.nModes
        yc = row_bottom(L, m) + L.rh / 2;
        txt(L.xA - 0.08, row_bottom(L, m) + L.rh, '100', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Tag', 'ytick100');
        txt(L.xA - 0.08, row_bottom(L, m), '0', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Tag', 'ytick0');
        parts = strsplit(D.mode_labels{m}, ' ');
        plot(ov, [0.005 0.145], [yc yc], '-', 'Color', L.dark * L.arch(m, :), 'LineWidth', L.lw, 'Tag', 'rowswatch');
        txt(0.17, yc, parts, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Tag', 'rowname');
        txt(L.xB + L.wB + 0.055, yc, sprintf('%d of 100', sum(D.finalCoopByMode{m} >= 0.95)), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Tag', 'count');
    end
    % two lines, so the label stays inside the 17.8 cm canvas
    txt(L.xB + L.wB * 0.95 - 0.02, yTop + 0.05, {'success', ['(' char(8805) ' 0.95)']}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Tag', 'success_label');
    draw_legend(ov, G, L, yHead + 0.05, yHead + 0.215);
end

function draw_legend(ov, G, L, y, yLead)
%DRAW_LEGEND  Grey lead-in line "selection on cooperation:" above one row of three swatches; asserted clear of B's letter.
    text(ov, L.xA, yLead, 'selection on cooperation:', 'FontSize', L.fs, 'FontName', figure_utils.FONT, ...
        'Color', [0.20 0.20 0.20], 'VerticalAlignment', 'middle', 'Tag', 'legendlead');
    C = state_colors(L);
    x = L.xA;  sw = 0.08;  swH = 0.08;
    for k = 1:3
        h = text(ov, 0, 0, G.names{k}, 'FontSize', L.fs, 'FontName', figure_utils.FONT, 'Color', L.base, ...
            'VerticalAlignment', 'middle', 'Tag', 'legendtext');
        w = h.Extent(3);
        rectangle(ov, 'Position', [x y - swH/2 sw swH], 'FaceColor', C(k, :), ...
            'EdgeColor', [0.45 0.45 0.45], 'LineWidth', 0.5);
        h.Position = [x + sw + 0.03, y, 0];
        x = x + sw + 0.03 + w + 0.07;
    end
    xEnd = x - 0.07;
    bL = findall(ov, 'Tag', 'letterB');
    if xEnd > bL.Extent(1) - 0.05
        error('fig2_draw:legend', 'Legend (ends %.2f in) comes within 0.05 in of the B letter (%.2f in)', xEnd, bL.Extent(1));
    end
    fprintf(['LEGEND: swatch row ends %.2f in, B letter starts %.2f in' newline], xEnd, bL.Extent(1));
end

function check_alignment(rowsA, rowsB)
%CHECK_ALIGNMENT  Row k of A and strip k of B: equal bottom and top; bottom axes on one baseline.
    for m = 1:numel(rowsA)
        pa = rowsA(m).ax.Position;  pb = rowsB(m).ax.Position;
        if abs(pa(2) - pb(2)) > 1e-6 || abs((pa(2) + pa(4)) - (pb(2) + pb(4))) > 1e-6
            error('fig2_draw:align', 'Row %d: A and B strip differ in bottom or top', m);
        end
        if m > 1 && rowsA(m).ax.Position(2) <= rowsA(m-1).ax.Position(2) + rowsA(m-1).ax.Position(4)
            error('fig2_draw:align', 'Rows %d and %d touch or overlap', m-1, m);
        end
    end
    if abs(rowsA(1).ax.Position(2) - rowsB(1).ax.Position(2)) > 1e-6, error('fig2_draw:align', 'Bottom axes differ'); end
    if abs(rowsA(1).ax.Position(3) - rowsB(1).ax.Position(3)) > 1e-6, error('fig2_draw:align', 'Plot widths differ'); end
    fprintf('CONTROL ALIGN PASS: rows of A and strips of B share bottom and top (1e-6), one bottom baseline, equal widths, gaps positive\n');
end

function check_encoding(rowsA, rowsB, ov, D)
%CHECK_ENCODING  Per row: the row name; the curve and swatch solid in the row's color; the bars filled
%   and edged in the architecture color, which names the architecture.
    L = layout_in();
    out = {'graded', 'graded', 'saturating', 'saturating'};   % rows 1-4 in fig2_data's architecture order
    sens = {'separate', 'contrast', 'separate', 'contrast'};
    for m = 1:D.nModes
        if ~strcmp(D.mode_labels{m}, [out{m} ' ' sens{m}]), error('fig2_draw:encoding', 'Row %d label %s unexpected', m, D.mode_labels{m}); end
        col = L.dark * L.arch(m, :);
        ln = findall(rowsA(m).ax, 'Tag', 'curve_line');
        if ~strcmp(ln.LineStyle, '-') || ~isequal(ln.Color, col), error('fig2_draw:encoding', 'Curve style/colour wrong in row %d', m); end
        bars = findall(rowsB(m).ax, 'Tag', 'hist_bars');
        if ~isequal(bars.FaceColor, L.arch(m, :)) || ~isequal(bars.EdgeColor, L.arch(m, :)), error('fig2_draw:encoding', 'Bars wrong in row %d', m); end
    end
    sw = findall(ov, 'Tag', 'rowswatch');
    ys = arrayfun(@(x) x.YData(1), sw);  [~, o] = sort(ys);  sw = sw(o);   % bottom row first
    for m = 1:D.nModes
        if ~strcmp(sw(m).LineStyle, '-') || ~isequal(sw(m).Color, L.dark * L.arch(m, :)) || sw(m).LineWidth ~= L.lw
            error('fig2_draw:encoding', 'Row swatch %d wrong', m);
        end
    end
    fprintf('CONTROL ENCODING PASS: per row, curve and swatch color and solid line, and filled bars in the architecture color\n');
end

function r = ternary_(c, a, b), if c, r = a; else, r = b; end, end

function check_clearances(ov, rowsB, D)
%CHECK_CLEARANCES  Tick-label clearance to strips, strip gaps, no overlapping annotation texts, scale bar clear of bars.
    L = layout_in();  tol = 1e-9;
    % (1) y tick labels vs the strips (overlay inches): bounding-box distance >= 0.03 in
    lab = [findall(ov, 'Tag', 'ytick100'); findall(ov, 'Tag', 'ytick0')];
    minClear = Inf;
    for m = 1:D.nModes
        sb = [L.xA, row_bottom(L, m) + L.rh + L.sgap, L.wA, L.strip];   % strip box [x y w h]
        for k = 1:numel(lab)
            e = lab(k).Extent;
            dx = max([sb(1) - (e(1) + e(3)), e(1) - (sb(1) + sb(3)), 0]);
            dy = max([sb(2) - (e(2) + e(4)), e(2) - (sb(2) + sb(4)), 0]);
            minClear = min(minClear, hypot(dx, dy));
        end
    end
    if minClear < 0.03, error('fig2_draw:clear', 'Tick label within %.3f in of a strip', minClear); end
    % (2) strip gaps: own gap < 1/4 of the gap to the previous row's x axis
    own = L.sgap;  prev = L.gap - L.sgap - L.strip;
    if own > prev / 4 + tol, error('fig2_draw:gaps', 'Strip gap %.3f vs %.3f is not below 1:4', own, prev); end
    % (3) no overlaps among overlay texts (0.01 in margin)
    t = findall(ov, 'Type', 'text');  t = t(arrayfun(@(x) ~isempty(strtrim(strjoin(string(x.String), ''))), t));
    for a = 1:numel(t)
        for b = a+1:numel(t)
            ea = t(a).Extent;  eb = t(b).Extent;
            if ea(1) < eb(1) + eb(3) + 0.01 && eb(1) < ea(1) + ea(3) + 0.01 && ea(2) < eb(2) + eb(4) + 0.01 && eb(2) < ea(2) + ea(4) + 0.01
                error('fig2_draw:overlap', 'Texts overlap: "%s" and "%s"', strjoin(string(t(a).String), ' '), strjoin(string(t(b).String), ' '));
            end
        end
    end
    % (4) the scale bar and its label (top strip of B) clear of the histogram bars, and "success" clear of the top strip
    axT = rowsB(D.nModes).ax;
    sc = findall(axT, 'Type', 'text');  sc = sc(arrayfun(@(x) strcmp(x.String, '20 populations'), sc));
    cnt = histcounts(D.finalCoopByMode{D.nModes}, linspace(0, 1, 21));
    xr = [0.12 - 0.012, sc.Extent(1) + sc.Extent(3)];
    binE = linspace(0, 1, 21);  hit = binE(1:end-1) < xr(2) & binE(2:end) > xr(1);
    barTop = max(cnt(hit));  yLow = min(sc.Extent(2), 0.30 * max(axT.YLim));
    if barTop >= yLow - 1, error('fig2_draw:scalebar', 'Scale bar/label reaches histogram bars (%d vs %.1f)', barTop, yLow); end
    fprintf(['CONTROL CLEARANCES PASS: tick labels clear of strips by %.3f in; strip gaps %.3f in own vs %.3f in to the previous row (1:%.1f); %d annotation texts mutually clear; scale bar clear of bars (max bar %d counts under the bar base %.1f)' newline], ...
        minClear, own, prev, prev / own, numel(t), barTop, yLow);
end

function check_polish(ov, rowsA, rowsB, D)
%CHECK_POLISH  Axis line per row, tick-label clearance to the axis ticks, strip outline, tick types, text colors.
    L = layout_in();
    for m = 1:D.nModes
        ln = findall(rowsA(m).ax, 'Tag', 'yaxis_line');
        if numel(ln) ~= 1 || ln.LineWidth ~= 0.6 || ~isequal(ln.YData, [0 100]), error('fig2_draw:polish', 'Row %d y-axis line missing or wrong', m); end
        ol = findall(rowsA(m).strip, 'Tag', 'strip_outline');
        if numel(ol) ~= 1 || ol.LineWidth ~= 0.3 || ~isequal(ol.EdgeColor, figure_utils.GRAY), error('fig2_draw:polish', 'Strip %d outline missing or wrong', m); end
        wantMinor = ternary_(m == 1, 'on', 'off');
        if ~strcmp(rowsA(m).ax.XMinorTick, wantMinor) || ~strcmp(rowsB(m).ax.XMinorTick, 'off'), error('fig2_draw:polish', 'Minor ticks wrong in row %d', m); end
        if ~isequal(rowsA(m).ax.XTick, [1 10 100 500]) || ~isequal(rowsB(m).ax.XTick, [0 0.5 1]), error('fig2_draw:polish', 'Major ticks wrong in row %d', m); end
    end
    tickEnd = L.xA - L.tickLen * max(L.wA, L.rh);                  % where the outward tick marks end (inches)
    lab = [findall(ov, 'Tag', 'ytick100'); findall(ov, 'Tag', 'ytick0')];
    clr = min(arrayfun(@(x) tickEnd - (x.Extent(1) + x.Extent(3)), lab));
    if clr < 0.03, error('fig2_draw:polish', 'y tick labels only %.3f in from the tick marks', clr); end
    t = findall(ov, 'Type', 'text');
    for k = 1:numel(t)
        if isempty(strtrim(strjoin(string(t(k).String), ''))), continue; end
        if max(t(k).Color) > 0.20 + 1e-9, error('fig2_draw:polish', 'Text "%s" lighter than 0.20', strjoin(string(t(k).String), ' ')); end
    end
    fprintf(['CONTROL POLISH PASS: y-axis line (0.6 pt) on all 4 rows; strip outlines 0.3 pt; tick labels %.3f in clear of the tick marks; minor ticks only on the bottom row, upper rows and B strips major only; all annotation text at or below 0.20' newline], clr);
end

function check_strip_strength(rowsA, rowsB, D)
%CHECK_STRIP_STRENGTH  Histogram FaceAlpha = L.histAlpha (0.55); strip colors = figure_utils.DILEMMA_HOLDS_TINT / DILEMMA_ESCAPED_TINT, opaque.
    L = layout_in();
    for m = 1:D.nModes
        bars = findall(rowsB(m).ax, 'Tag', 'hist_bars');
        if bars.FaceAlpha ~= L.histAlpha, error('fig2_draw:strength', 'Histogram FaceAlpha %.2f in row %d, expected %.2f', bars.FaceAlpha, m, L.histAlpha); end
        p = findall(rowsA(m).strip, 'Tag', 'bg_patch');
        for k = 1:numel(p)
            if p(k).FaceAlpha ~= 1, error('fig2_draw:strength', 'Strip patch is transparent'); end
        end
    end
    wantH = figure_utils.DILEMMA_HOLDS_TINT;
    wantE = figure_utils.DILEMMA_ESCAPED_TINT;
    C = state_colors(L);
    if max(abs(C(1, :) - wantH)) > 1e-12 || max(abs(C(3, :) - wantE)) > 1e-12 || ~isequal(C(2, :), L.tintNone)
        error('fig2_draw:strength', 'State colors are not figure_utils.DILEMMA_HOLDS_TINT / DILEMMA_ESCAPED_TINT');
    end
    fprintf(['CONTROL STRIP STRENGTH PASS: histogram FaceAlpha %.2f in all 4 rows; strip blocks opaque (holds %s, escaped %s), none %s' newline], ...
        L.histAlpha, mat2str(round(C(1, :), 3)), mat2str(round(C(3, :), 3)), mat2str(C(2, :)));
end

function check_counts(ov, D)
%CHECK_COUNTS  Success counts read back from the drawn text.
    h = findall(ov, 'Tag', 'count');
    ys = arrayfun(@(x) x.Position(2), h);  [~, o] = sort(ys);   % bottom row first = architecture 1
    got = arrayfun(@(x) sscanf(x.String, '%d of 100'), h(o))';
    want = arrayfun(@(m) sum(D.finalCoopByMode{m} >= 0.95), 1:D.nModes);
    if ~isequal(got, want), error('fig2_draw:counts', 'Success counts %s differ from %s', mat2str(got), mat2str(want)); end
    fprintf('CONTROL COUNTS PASS: success counts %s read back from the drawn text\n', mat2str(got));
end

function check_text_and_lines(fig)
%CHECK_TEXT_AND_LINES  Printed text 6.5-12 pt (canvas = print width); data lines >= 0.75 pt; axes 0.5-0.75 pt.
    t = findall(fig, 'Type', 'text');
    n = 0;
    for k = 1:numel(t)
        if isempty(strtrim(strjoin(string(t(k).String), ''))), continue; end
        if t(k).FontSize < 6.5 || t(k).FontSize > 12, error('fig2_draw:text', 'Text size %.1f pt: %s', t(k).FontSize, strjoin(string(t(k).String), ' ')); end
        n = n + 1;
    end
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        if ax(k).FontSize < 6.5 || ax(k).FontSize > 12, error('fig2_draw:text', 'Axes tick font %.1f pt', ax(k).FontSize); end
        if ~strcmp(ax(k).Tag, 'legend_overlay') && (ax(k).LineWidth < 0.5 || ax(k).LineWidth > 0.75)
            error('fig2_draw:lines', 'Axes line width %.2f', ax(k).LineWidth);
        end
    end
    for tag = {'curve_line', 'hist_bars', 'success_line'}
        h = findall(fig, 'Tag', tag{1});
        for k = 1:numel(h)
            if isprop(h(k), 'LineWidth') && h(k).LineWidth < 0.75, error('fig2_draw:lines', '%s line width %.2f', tag{1}, h(k).LineWidth); end
        end
    end
    fprintf('CONTROL TEXT/LINES PASS: %d text objects and %d axes within 6.5-12 pt; data lines >= 0.75 pt; axes 0.5-0.75 pt\n', n, numel(ax));
end

function check_patches(rows, G)
%CHECK_PATCHES  Drawn states equal computed states; holds and escaped patches identical across rows and equal to the state colors.
    L = layout_in();  C = state_colors(L);
    seen = zeros(0, 4);
    for m = 1:numel(rows)
        axp = rows(m).ax;  if L.strip > 0, axp = rows(m).strip; end
        p = findall(axp, 'Tag', 'bg_patch');
        drawn = nan(500, 1);
        for k = 1:numel(p)
            u = p(k).UserData;
            if ~isequal(p(k).FaceColor, C(u(3) + 2, :)), error('fig2_draw:patches', 'Patch colour mismatch (row %d)', m); end
            seen = [seen; u(3), p(k).FaceColor]; %#ok<AGROW>
            if any(~isnan(drawn(u(1):u(2)))), error('fig2_draw:patches', 'Overlapping patches'); end
            drawn(u(1):u(2)) = u(3);
        end
        if ~isequal(drawn, G.state(:, m)), error('fig2_draw:patches', 'Drawn states differ from computed states (%s)', G.modes{m}); end
    end
    for st = [-1 1]
        if size(unique(seen(seen(:, 1) == st, 2:4), 'rows'), 1) ~= 1, error('fig2_draw:patches', 'State %d patches differ across rows', st); end
    end
    fprintf('CONTROL PATCHES PASS: drawn states equal computed states at all 500 generations x 4 rows; holds and escaped patches identical across rows and equal to the figure_utils tints' + string(char(10)));
end

function check_curves_unchanged(rows, D)
%CHECK_CURVES_UNCHANGED  The curves equal the loader's.
    for m = 1:numel(rows)
        ln = findall(rows(m).ax, 'Tag', 'curve_line');
        if ~isequal(ln.YData, D.meanCoopOverTime(m, :))
            error('fig2_draw:curves', 'Curve data of architecture %d differ from the loader', m);
        end
    end
    fprintf('CONTROL CURVES PASS: four curves identical to the loader\n');
end

function G = load_selection_background(cfg, D)
%LOAD_SELECTION_BACKGROUND  Ensemble inherited selection S_G, bootstrap CI and state per architecture and generation.
    dataDir = cfg.dataDir;
    if isfolder(fullfile(dataDir, 'publication')), dataDir = fullfile(dataDir, 'publication'); end   % the data root, or already its publication folder
    ed = fullfile(dataDir, 'selection_over_time');
    modes = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'};
    bcfg = struct('bootN', 10000, 'bootSeed', 20261022);
    A = cell(D.nModes, 1);
    for m = 1:D.nModes
        A{m} = zeros(500, 20);
        for r = 1:20
            tag = sprintf('%s_r%02d', modes{m}, r);
            Zf = load(fullfile(ed, ['selection_z_' tag '.mat']), 'Z', 'par');
            Sf = load(fullfile(ed, ['selection_sign_' tag '.mat']), 'S');
            Rf = load(fullfile(ed, ['lineage_records_' tag '.mat']), 'par', 'idx');
            if ~isequal(Zf.par, Rf.par), error('fig2_draw:join', '%s: parent records differ', tag); end
            S = zeros(500, 1);  SG = S;
            for g = 1:500
                [S(g), SG(g)] = parts_of(Zf.par(:, g), Zf.Z(:, g), Rf.idx(:, g));
            end
            if max(abs(S - Sf.S)) > 1e-12, error('fig2_draw:join', '%s: S differs from stored S', tag); end
            A{m}(:, r) = SG;
        end
    end
    ci = boot_all(A, bcfg);
    G.modes = modes;  G.mS = zeros(500, D.nModes);  G.ci = zeros(500, 2, D.nModes);  G.state = zeros(500, D.nModes);
    for m = 1:D.nModes
        G.mS(:, m) = mean(A{m}, 2);
        G.ci(:, :, m) = ci{m};
        st = zeros(500, 1);  st(ci{m}(:, 1) > 0) = 1;  st(ci{m}(:, 2) < 0) = -1;
        G.state(:, m) = st;
    end
    G.names = {'dilemma holds', 'no detectable selection', 'dilemma escaped'};
end

function ci = boot_all(Sc, cfg)
% Bootstrap over runs, as in selection_inherited: per generation, the 2.5 and 97.5 percentiles of
% the mean over resampled runs. The caller's random stream is restored afterwards.
    rs = rng; restore = onCleanup(@() rng(rs));
    rng(cfg.bootSeed);
    ci = cell(numel(Sc), 1); chunk = 1000;
    for a = 1:numel(Sc)
        S = Sc{a}; nr = size(S, 2); G = size(S, 1);
        bm = zeros(G, cfg.bootN);
        for b0 = 1:chunk:cfg.bootN
            nb = min(chunk, cfg.bootN - b0 + 1);
            idx = randi(nr, nr, nb);
            X = reshape(S(:, idx(:)), G, nr, nb);
            bm(:, b0:b0 + nb - 1) = reshape(mean(X, 2), G, nb);
        end
        ci{a} = prctile(bm, [2.5 97.5], 2);
    end
end

function check_summary(G, cfg)
%CHECK_SUMMARY  Checks against summary_selection_inherited.txt: (1) the first generation favoring
%   cooperation in each architecture; (2) the S_G means and intervals of its table, 48 rows, to
%   0.5e-4; (3) the shares of generations 16 to 500 in each state. The summary is read from the
%   same publication folder as the runs (cfg.dataDir).
    dataDir = cfg.dataDir;
    if isfolder(fullfile(dataDir, 'publication')), dataDir = fullfile(dataDir, 'publication'); end   % the data root, or already its publication folder
    f = fullfile(dataDir, 'selection_over_time', 'summary_selection_inherited.txt');
    txt = fileread(f);
    first = zeros(1, 4);
    for m = 1:4
        first(m) = find(G.state(:, m) == 1, 1);
        tk = regexp(txt, [G.modes{m} '\s+S_G (\d+) \|'], 'tokens', 'once');
        if isempty(tk) || str2double(tk{1}) ~= first(m), error('fig2_draw:first', 'summary first-generation line differs for %s', G.modes{m}); end
    end
    fprintf('CONTROL 1 PASS: first "cooperating pays" generation under S_G = %s, equal to the summary\n', mat2str(first));
    s0 = strfind(txt, 'Ensemble mean [95% bootstrap CI]');
    lines = strsplit(txt(s0(1):end), newline);
    cur = 0;  n = 0;  worst = 0;
    for i = 2:numel(lines)
        nm = regexp(lines{i}, '^\s+(\w+) \(\d+ runs\)', 'tokens', 'once');
        if ~isempty(nm), cur = find(strcmp(G.modes, nm{1}), 1); continue; end
        tk = regexp(lines{i}, '^\s+(\d+)\s+[+-]\d\.\d{4} \[[^\]]*\]\s+([+-]\d\.\d{4}) \[([+-]\d\.\d{4}),([+-]\d\.\d{4})\]', 'tokens', 'once');
        if isempty(tk) || cur == 0, continue; end
        g = str2double(tk{1});  v = str2double(tk(2:4));
        d = max(abs([G.mS(g, cur), G.ci(g, :, cur)] - v));
        worst = max(worst, d);
        if d > 0.5e-4 + 1e-9, error('fig2_draw:report', '%s g %d: report differs by %.2e', G.modes{cur}, g, d); end
        n = n + 1;
    end
    if n ~= 48, error('fig2_draw:report', 'Checked %d report rows, expected 48', n); end
    fprintf('CONTROL 2 PASS: S_G mean and CI equal the summary table (%d rows, 4 architectures x 12 generations), max diff %.1e\n', n, worst);
    nS = 0;
    for m = 1:4
        s2 = G.state(16:500, m);
        got = [mean(s2 == -1), mean(s2 == 0), mean(s2 == 1)];
        tk = regexp(txt, [G.modes{m} '\s+S_G defecting pays ([\d.]+), no clear advantage ([\d.]+), cooperating pays ([\d.]+)'], 'tokens', 'once');
        if isempty(tk) || max(abs(got - str2double(tk))) > 0.0005 + 1e-9
            error('fig2_draw:shares', 'State shares (generations 16-500) differ from the summary for %s', G.modes{m});
        end
        nS = nS + 1;
    end
    fprintf('CONTROL 3 PASS: state shares over generations 16-500 equal the summary S_G shares (%d architectures)\n', nS);
end

