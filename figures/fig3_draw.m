function fig3_draw(cfg)
% FIG3_DRAW  Draws Fig. 3: success against the generation-500 share of members, one point per
% setting: the four signaling architectures (circles), five spatial settings (squares) and groups
% split by payoff or at random (diamonds), the last two families under saturating contrast.
%
% Implements : Fig. 3 (called by fig3_success_vs_share_of_members)
% Inputs     : cfg   struct; cfg.outStem, the output path without extension
%              the signaling and spatial points, from fig3_data('data')
%              publication/css_predicts_spatial.mat   (condition_names; written by
%                gen500_share_of_members_spatial)
%              publication/group_two_level/group_fission_analysis_N1000.mat   (struct A:
%                the random-split control, `fission_sham`, and group selection, `fission`;
%                written by group_selection_arms)
% Outputs    : [outStem '.pdf'] (vector, 8.7 x 7.05 cm), '.png' (300 dpi) and '.tif' (600 dpi,
%              LZW, through figure_utils.export_tiff); prints r, each point with its interval, and
%              the fonts used
% Seed       : none (no random draws)
%
% Intervals are exact (Clopper-Pearson) 95% binomial intervals on success, with 100 runs per
% circle and square and each arm's stored count per diamond. r, the correlation over all points,
% is computed before two points near a share of 0.20 are drawn 0.006 to either side (ring k = 2
% left, group selection right) so both stay visible; the printed list gives the unshifted shares.

    W = 8.7; H = 7.05;   % figure size, cm
    figure_utils.setup_defaults();
    OD = fig3_data('data', figure_config());
    OD = add_parameters(OD);
    OD = load_group_arms(OD);

    fig = figure('Units','centimeters', 'Position',[1 1 W H], 'Color','w');
    set(fig, 'DefaultAxesFontName',figure_utils.FONT, 'DefaultTextFontName',figure_utils.FONT);
    guard = onCleanup(@() close(fig)); %#ok<NASGU>

    ax = axes(fig, 'Position', [1.25/W, 1.0/H, 7.15/W, 5.65/H]);   % axes box 7.15 x 5.65 cm
    draw_panel(ax, OD);
    report_fonts(fig);

    stem = cfg.outStem;
    w = warning('off','MATLAB:print:ContentTypeImageSuggested');
    exportgraphics(fig, [stem '.pdf'], 'ContentType','vector', 'Padding','figure');
    warning(w);
    exportgraphics(fig, [stem '.png'], 'Resolution', 300, 'Padding','figure');
    figure_utils.export_tiff(fig, [stem '.tif'], 600);
    fprintf('\nSaved %s.pdf + .png + .tif (%.1f x %.1f cm)\n', stem, W, H);
end


function OD = load_group_arms(OD)
%LOAD_GROUP_ARMS  the two group arms at N = 1000 (20 groups of 50): random
%   splits first, group selection second, as draw_panel and group_labels
%   expect; r over all points is computed here, before any drawing shift.
    cfg = figure_config();
    pubDir = fullfile(cfg.dataDir, 'publication');
    f = fullfile(pubDir, 'group_two_level', 'group_fission_analysis_N1000.mat');
    G = load(f); A = G.A;
    OD.grp.reach_occupancy   = [A.fission_sham.occ.occupancy; A.fission.occ.occupancy];
    OD.grp.success_empirical = 100 * [A.fission_sham.rate; A.fission.rate];
    OD.grp_n = [A.fission_sham.n; A.fission.n];
    OD.grp_names = {'random splits', 'group selection'};
    OD.x_all = [OD.sens.reach_occupancy(:); OD.spat.reach_occupancy(:); OD.grp.reach_occupancy(:)];
    OD.y_all = [OD.sens.success_empirical(:); OD.spat.success_empirical(:); OD.grp.success_empirical(:)];
    OD.r_pool = corr(OD.x_all, OD.y_all);
    fprintf('fig3_draw: %d points, pooled r = %+.4f\n', numel(OD.x_all), OD.r_pool);
end


function OD = add_parameters(OD)
%ADD_PARAMETERS  which topology each square is, read from the file the
%   spatial points come from (the diamonds' names are set in load_group_arms).
    cfg = figure_config();
    pubDir = fullfile(cfg.dataDir, 'publication');
    sp = load(fullfile(pubDir, 'css_predicts_spatial.mat'), 'condition_names');
    OD.spat_names = sp.condition_names(:)';
end


% ===================== spec constants ================================
function f = fonts()
    f.axis = 8; f.tick = 7; f.label = 7;
end

function c = ink()
    c.axis = figure_utils.INK;
    c.note = figure_utils.NOTE;
    c.key  = figure_utils.INK;
end

function m = marker_sizes()
%MARKER_SIZES  circle, square and diamond matched by enclosed area on the
%   page. At a given MarkerSize MATLAB draws a square smaller than a circle
%   and a diamond taller than wide (measured on a 600 dpi render: 9.1 pixels
%   per MarkerSize unit for a circle, 7.0 for a square, 9.9 wide by 13.1 tall
%   for a diamond), so equal area with a 4.5 pt circle needs MarkerSize 5.18
%   for the square and 4.51 for the diamond. One table for the data points,
%   the legend and the shape key.
    d = 4.5;
    m.o = d;
    m.s = 5.18;
    m.d = 4.51;
end


% ===================== the panel =====================================
function draw_panel(ax, OD)
    F = fonts(); C = ink(); hold(ax,'on'); set(ax,'Color','w');
    % Two of the three points near a share of 0.20 are drawn shifted so all
    % stay visible (ring k = 2 by -0.006, group selection by +0.006; drawing
    % only). load_group_arms computes r from the true values, and print_mapping
    % prints the true shares (OD.x_true).
    OD.x_true = OD.x_all;
    OD.spat.reach_occupancy(3) = OD.spat.reach_occupancy(3) - 0.006;
    OD.grp.reach_occupancy(2) = OD.grp.reach_occupancy(2) + 0.006;
    OD.x_all = [OD.sens.reach_occupancy(:); OD.spat.reach_occupancy(:); OD.grp.reach_occupancy(:)];
    cols = OD.palette.colors;
    satcon = cols(4,:);
    nrep = [100 * ones(numel(OD.x_all) - numel(OD.grp_n), 1); OD.grp_n];

    xl = [0, max(OD.x_all)*1.18]; yl = [0 100];
    set(ax, 'XLim', xl, 'YLim', yl, 'YTick', 0:25:100);

    [lo, hi] = clopper_pearson(OD.y_all, nrep);
    errorbar(ax, OD.x_all, OD.y_all, OD.y_all - lo, hi - OD.y_all, ...
        'LineStyle','none', 'Color',figure_utils.GRAY, 'LineWidth',0.8, ...
        'CapSize',2, 'HandleVisibility','off');

    for m = 1:4
        mark(ax, OD.sens.reach_occupancy(m), OD.sens.success_empirical(m), 'o', cols(m,:));
    end
    for c = 1:numel(OD.spat.reach_occupancy)
        mark(ax, OD.spat.reach_occupancy(c), OD.spat.success_empirical(c), 's', satcon);
    end
    for g = 1:numel(OD.grp.reach_occupancy)
        mark(ax, OD.grp.reach_occupancy(g), OD.grp.success_empirical(g), 'd', satcon);
    end

    spatial_labels(ax, OD, F, C);
    group_labels(ax, OD, F, C);
    architecture_legend(ax, cols, satcon, F, C);
    % r and n in the top-right corner (places the text; not a result).
    text(ax, xl(2) - 0.0263*diff(xl), 96.7, sprintf('r = %+.2f, n = %d', OD.r_pool, numel(OD.x_all)), ...
        'FontSize', F.label, 'Color', C.axis, ...
        'HorizontalAlignment','right', 'VerticalAlignment','top');
    shape_key(ax, xl, F, C);

    % Plain-words label: the fraction of the population that are members of
    % the conditional genotype network; no symbol or subscript to size.
    xlabel(ax, 'share of members at generation 500', ...
        'FontSize', F.axis, 'Color', C.axis);
    ylabel(ax, 'success rate (%)', 'FontSize', F.axis, 'Color', C.axis);
    style(ax);
    print_mapping(OD, lo, hi);
end


function mark(ax, x, y, shape, col, face)
%MARK  a matched-area glyph (marker_sizes) with a 0.8 pt edge in col, filled
%   with col unless face is given ('none' draws the outline only, as in the
%   shape key).
    if nargin < 6, face = col; end
    sizes = marker_sizes();
    switch shape
        case 'o', sz = sizes.o;
        case 's', sz = sizes.s;
        case 'd', sz = sizes.d;
        otherwise, error('mark: unknown shape %s', shape);
    end
    plot(ax, x, y, shape, 'MarkerSize', sz, 'MarkerFaceColor', face, ...
        'MarkerEdgeColor', col, 'LineWidth', 0.8);
end


function spatial_labels(ax, OD, F, C)
%SPATIAL_LABELS  one label per square: the neighborhood k, or "well mixed";
%   each label sits dx (data units) beside its square.
    o = OD.spat.reach_occupancy; s = OD.spat.success_empirical;
    dx = 0.0145;
    put(ax, o(1), 34.95, 'well mixed', 'center', F, C);    % wellmixed
    put(ax, o(2) + dx, s(2), '1',  'left',  F, C);         % ring k = 1
    put(ax, o(3) - dx, s(3), '2',  'right', F, C);         % ring k = 2, left of its square
    put(ax, o(4) - dx, s(4), '5',  'right', F, C);         % ring k = 5
    put(ax, o(5) + dx, s(5), '10', 'left',  F, C);         % ring k = 10
end


function group_labels(ax, OD, F, C)
%GROUP_LABELS  one direct label per diamond.
    o = OD.grp.reach_occupancy; s = OD.grp.success_empirical;
    dx = 0.0145;   % the squares' label gap
    % The diamond already marks groups; the labels name how the splitting
    % group is chosen.
    put(ax, o(1) + dx, s(1), 'split at random', 'left', F, C);
    put(ax, o(2) + dx, s(2), 'split by payoff', 'left', F, C);
end


function put(ax, x, y, str, ha, F, C)
    text(ax, x, y, str, 'FontSize', F.label, 'Color', C.key, ...
        'HorizontalAlignment', ha, 'VerticalAlignment','middle');
end


function architecture_legend(ax, cols, satcon, F, C)
%ARCHITECTURE_LEGEND  hand drawn, because one row carries three glyphs: the
%   saturating-contrast color, which the squares and diamonds also wear,
%   shows circle, square and diamond, then the name. All four names share one
%   x, clear of that row's three glyphs; the three upper rows leave the space
%   blank. Positions are data units (places the legend; not a result).
    names = figure_utils.architecture_palette().labels;
    gx = 0.0604;                      % circle column
    sx = 0.0781;                      % square column
    dx = 0.0982;                      % diamond column (diamond is the widest
                                       % glyph, so it needs the widest gap in
                                       % front of it to clear the square by
                                       % the same margin the other pairs get)
    tx = 0.1133;                      % common name column, after all three
    yr = [94.5 87.3 80.1 72.9];
    for m = 1:4
        mark(ax, gx, yr(m), 'o', cols(m,:));
        if m == 4
            % Only the saturating-contrast row wears the square and the
            % diamond: that color is the one shared with the spatial and
            % group-selection points, so this is the one row that needs them.
            mark(ax, sx, yr(m), 's', cols(m,:));
            mark(ax, dx, yr(m), 'd', cols(m,:));
        end
        text(ax, tx, yr(m), names{m}, 'FontSize', F.label, 'Color', C.key, ...
            'VerticalAlignment','middle');
    end
    if ~isequal(cols(4,:), satcon), error('palette mismatch'); end
end


function shape_key(ax, xl, F, C)
%SHAPE_KEY  what each shape is and at what population size; no header, to
%   match the legend above. gx sits in the empty space between the low- and
%   high-share points (places the key; not a result). The three glyphs are
%   outlines only, edge dark grey (C.note) at the markers' 0.8 pt, so the key
%   reads as shape, apart from the filled data points.
    shapes = {'o','s','d'};
    % Shape and population size only; what the numerals beside the squares
    % are, and the group sizes of the diamonds, are left to the caption.
    labels = {'signaling, pop. 1000', 'spatial, pop. 50', 'groups, pop. 1000'};
    gx = 0.275; tx = gx + 0.0276; yr = [20.7 13.0 5.3];
    for i = 1:3
        mark(ax, gx, yr(i), shapes{i}, C.note, 'none');
        text(ax, tx, yr(i), labels{i}, 'FontSize', F.label, 'Color', C.note, ...
            'VerticalAlignment','middle');
    end
    if xl(2) < 0.5, warning('key may overrun: xlim %.3f', xl(2)); end
end


function [lo, hi] = clopper_pearson(pct, n)
%CLOPPER_PEARSON  exact binomial 95% interval, in percent. The point estimate
%   is untouched; this only puts an interval around it.
    n = n(:) .* ones(numel(pct), 1);
    k = round(pct(:) .* n / 100);
    lo = zeros(size(k)); hi = zeros(size(k));
    for i = 1:numel(k)
        if k(i) > 0, lo(i) = betainv(0.025, k(i), n(i) - k(i) + 1) * 100; end
        if k(i) < n(i), hi(i) = betainv(0.975, k(i) + 1, n(i) - k(i)) * 100; else, hi(i) = 100; end
    end
end


function print_mapping(OD, lo, hi)
%PRINT_MAPPING  every point, its parameter and its interval, for the caption.
    names = [OD.palette.labels, OD.spat_names(:)', OD.grp_names(:)'];
    fprintf('\nFig 3 point mapping (share of members, success, 95%% CI):\n');
    for i = 1:numel(OD.x_true)
        fprintf('  %-22s %.4f  %2.0f%%  [%4.1f, %4.1f]\n', names{i}, ...
            OD.x_true(i), OD.y_all(i), lo(i), hi(i));
    end
end


function style(ax)
%STYLE  0.6 pt axes, ticks out, no box, no grid; 7 pt ticks and 8 pt axis labels.
    F = fonts(); C = ink();
    set(ax, 'FontName',figure_utils.FONT, 'FontSize', F.tick, 'LineWidth', 0.6, ...
        'TickDir','out', 'TickLength',[0.015 0.015], 'Box','off', ...
        'XColor', C.axis, 'YColor', C.axis, ...
        'XTickLabelRotation', 0, 'YTickLabelRotation', 0);
    % Setting the axes font size rescales the axis labels through
    % LabelFontSizeMultiplier, which is 1.1 by default and turned 8 pt labels
    % into 7.7. Pin the multiplier and restate the size.
    set(ax, 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
    set([ax.XLabel ax.YLabel], 'FontSize', F.axis, 'Color', C.axis, 'FontName',figure_utils.FONT);
end


function report_fonts(fig)
%REPORT_FONTS  every distinct font family, size and weight in the figure, before export.
    t = findall(fig, 'Type','text');
    a = findall(fig, 'Type','axes');
    sz = [arrayfun(@(h) h.FontSize, t); arrayfun(@(h) h.FontSize, a)];
    fam = unique([arrayfun(@(h) string(h.FontName), t); ...
                  arrayfun(@(h) string(h.FontName), a)]);
    wt = unique(arrayfun(@(h) string(h.FontWeight), t));
    fprintf('\nFONTS: family %s | sizes %s | weights %s\n', ...
        strjoin(cellstr(fam), ', '), mat2str(unique(sz)'), strjoin(cellstr(wt), ', '));
end
