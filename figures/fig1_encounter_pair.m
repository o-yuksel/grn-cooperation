function fig1_encounter_pair(dataFile)
% FIG1_ENCOUNTER_PAIR  Draws Fig. 1: two networks coupled through their sensors (A); network 1's
% readout over development when network 2 is a copy of it (conditional) and when it is an
% unconditional defector (B); and a two-by-two table of what the two networks of one encounter
% gain by defecting on each other, by kind of pair (C).
%
% Produces : Fig. 1
% Inputs   : dataFile   (optional) a saved cache of fig1_partner_data (struct out), read only;
%                       default: fig1_partner_data(), which reads or builds its cache
%                       output/figures/fig1_partner_data.mat from the member that
%                       fig1_member_and_mutant picks in publication/saturating_contrast_100reps.mat
%                       (written by run_signaling_architectures)
% Outputs  : output/figures/Fig1.pdf (vector), .png (300 dpi) and .tif (600 dpi, RGB, LZW), under
%            the release root; prints the layout, the fonts and node fills, and the gaps between
%            texts, curves, dots and grid lines (pt)
% Run      : fig1_encounter_pair
% Seed     : none here (the partner search in fig1_partner_data seeds itself)
%
% The figure is 17.8 x 7.2 cm, Arial, 7 pt text and ticks, 8 pt axis titles. Panel B plots
% development steps 1 to 300 (log axis), the horizon the runs develop for: solid, network 2 a copy
% of network 1; dashed, network 2 an unconditional defector; each curve ends in a dot colored by
% network 1's action (COLOR_COOP cooperates, COLOR_DEFECT defects). Panel C's cell texts follow
% from the run's payoffs T, R, P, S: the incentive to defect on a network is what the other earns
% by defecting on it minus what it earns by cooperating with it, P - R = -2 on a conditional
% network, which defects back, and P - S = +1 or T - R = +2 on an unconditional one (a defector or
% a cooperator); the run stops unless the payoffs give -2, +1 and +2. A network gains by defecting
% where its partner is unconditional. The both-conditional cell ("dilemma escaped") and the
% both-unconditional cell ("dilemma holds") take figure_utils.DILEMMA_ESCAPED_TINT and
% DILEMMA_HOLDS_TINT, the fills of Fig. 2A's strip; the two mixed cells are white. The table shows
% single pairs, not averages over encounters, and does not draw a network that cooperates with a
% defector and defects on a cooperator (T - S = +5).

    if nargin < 1, dataFile = ''; end
    root = fileparts(fileparts(mfilename('fullpath')));
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outBase = fullfile(outDir, 'Fig1');

    if isempty(dataFile)
        dat = fig1_partner_data();
    else
        S = load(dataFile, 'out'); dat = S.out;
        fprintf('fig1_encounter_pair: data read from %s\n', dataFile);
    end

    fig = draw_figure(dat);

    w = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
    warning(w);
    exportgraphics(fig, [outBase '.png'], 'Resolution', 300);
    figure_utils.export_tiff(fig, [outBase '.tif'], 600);
    close(fig);
    fprintf('fig1_encounter_pair: wrote %s.{pdf,png,tif}  (%.2f x %.2f in canvas)\n', ...
        outBase, 7.00, 7.2/2.54);
end

% =========================================================================
function fig = draw_figure(dat)
%DRAW_FIGURE  The three panels and the top strip on one 17.8 x 7.2 cm canvas.
    cfg = struct();
    cfg.figWidth  = 7.00;            % 17.8 cm
    cfg.figHeight = 7.2/2.54;        % 7.2 cm

    % --- type and strokes ---------------------------------------------------
    cfg.fsPanel  = 10;      % panel letters, the only bold text
    cfg.fsAxis   = 8;       % axis titles
    cfg.fsTick   = 7;       % tick labels
    cfg.fsText   = 7;       % every in-panel text
    cfg.colAxis  = figure_utils.INK;
    cfg.colNote  = figure_utils.NOTE;
    cfg.colFaint = [0.50 0.50 0.50];
    cfg.colRef   = [0.50 0.50 0.50];    % reference lines
    cfg.colLink  = [0.75 0.75 0.75];    % C's lines from B's end dots
    cfg.lwAxis   = 0.6;
    cfg.lwData   = 1.2;
    cfg.lwRef    = 0.6;
    cfg.lwLink   = 0.5;
    cfg.lwSchem  = 0.8;     % schematic strokes (arrows, gene rims)
    cfg.mkSize   = 4.5;
    cfg.mkEdge   = 0.8;
    cfg.gapPt    = 2.0;     % white between an arrow or line end and the outer
                            % edge of the stroke it meets (pt)
    cfg.fanPt    = 3.0;     % the coupling's fan point sits this far out (pt)
                            % so the two arrows into a sensor pair keep a stem
    cfg.stickGapPt = 2.4;   % the strip's stick stops this far below the text
                            % box of "development" (pt)
    cfg.headTip  = 1.0;    % an arrowhead's stroked point overshoots its end
                            % point by about this much (pt); added at every tip

    % --- color ---------------------------------------------------------------
    cfg.coop     = figure_utils.COLOR_COOP;     % #0072B2 cooperation
    cfg.defect   = figure_utils.COLOR_DEFECT;   % #D55E00 defection
    cfg.ink      = figure_utils.INK;            % panel B's two curves, C's dilemma lines
    cfg.fsNodeFill = 7;                         % white letters on the filled nodes
    % White on #D55E00 has 3.9:1 contrast against 5.2:1 on #0072B2, so the
    % two D letters are one point larger.
    cfg.fsNodeD    = 8;
    % Panel C's two filled cells = Fig. 2A's strip fills.
    cfg.tintHolds   = figure_utils.DILEMMA_HOLDS_TINT;     % #E0BFB4
    cfg.tintEscaped = figure_utils.DILEMMA_ESCAPED_TINT;   % #C1CEE9

    % --- layout (figure fractions; pt on the 504 pt canvas where noted) -----
    W = cfg.figWidth * 72;
    cfg.yWords   = 0.89731;  % the life-cycle line
    cfg.yDrop    = 0.83969;  % the bracket under panels A and B
    cfg.letterY  = 0.78460;  % panel letters, one baseline

    % B's box; A's two networks are fitted to its height (network 1's header
    % on B's top edge, network 2's nodes on B's baseline). The horizontal
    % positions are then set from panel C's measured labels (below).
    bLeft = 0.435 - 0.018;
    cfg.posR     = [bLeft, 0.188, (343 - bLeft*W)/W, 0.5465];
    topEdge = cfg.posR(2) + cfg.posR(4);
    botEdge = cfg.posR(2);
    ah = 0.371;
    cfg.posA = [0.078, topEdge - 0.9405*ah, 0.2587, ah];
    cfg.posP = [0.078, botEdge - 0.3571*ah, 0.2587, ah];

    % Panel B legend and key words (pt from B's lower-left corner).
    cfg.legX      = 5;      % legend's left edge
    cfg.legBottom = 5;      % bottom of the legend's last text line
    cfg.dashPt    = [5 2.5];   % panel B's dashed curve: dash, gap (pt along the curve)
    cfg.legSample = 12.5;   % line-sample length: two dashes and one gap
    cfg.legGap    = 3;      % sample to text
    % Legend rows, top to bottom: {line style ('' = heading, no sample), text lines}.
    cfg.legRows   = {{'', {'network 2 is'}}, {'-', {'conditional'}}, ...
                     {'--', {'an unconditional defector'}}};
    cfg.headGapPt = 3;      % ink gap from B's heading to B's top edge (pt)
    cfg.bYTicks   = [-1 0 1];
    cfg.bYlim     = [-1.85 2.05];
    cfg.keyLinePt = 8.6;    % baseline-to-baseline step of the stacked key words (pt)
    cfg.keyGapPt  = 1.5;    % ink gap between the dotted line and the words next to it (pt)
    cfg.sideRight = 7;      % right edge of "action threshold", pt in from B's right end
    cfg.inkAsc    = 0.733;  % ink height of a lower-case word with ascenders ("defects"), em,
                            % as measured on the 300 dpi render

    cfg.table = table_layout(dat.config);

    figure_utils.setup_defaults();
    fig = figure_utils.create_figure(cfg.figWidth, cfg.figHeight);
    set(fig, 'GraphicsSmoothing', 'on');
    cfg.ptW = cfg.figWidth * 72;     % points per normalized figure unit
    cfg.ptH = cfg.figHeight * 72;

    % Horizontal layout. Panel C's labels are measured on a stand-in table
    % first: the rotated row labels and "network 2" to its right (Rr) and the
    % depth of its column labels and "network 1" below it. A's two networks and B move
    % L.shiftPt left, B then L.bShiftPt right; the table's left edge stands
    % L.abGap pt right of B's end dots.
    Lt = cfg.table; rD = cfg.mkSize/2 + cfg.mkEdge/2;
    cfg.posC = [400/W, botEdge, 80/W, 80/cfg.ptH];
    [axM, Sst, depth, Rr] = draw_table(fig, cfg, true); delete(axM);
    sh = Lt.shiftPt/W;
    cfg.posA(1) = cfg.posA(1) - sh; cfg.posP(1) = cfg.posP(1) - sh; cfg.posR(1) = cfg.posR(1) - sh;
    cfg.posR(1) = cfg.posR(1) + Lt.bShiftPt/W;
    bEndPt = (cfg.posR(1) + cfg.posR(3))*W;
    gL = bEndPt + rD + Lt.abGap;
    fprintf('LAYOUT: B axes %.1f..%.1f pt (%.1f wide); C labels: left %.1f, right %.1f, depth %.1f pt\n', ...
        cfg.posR(1)*W, bEndPt, cfg.posR(3)*W, Sst, Rr, depth);

    infoA = draw_agent(fig, cfg.posA, 'network 1', false, cfg);
    infoP = draw_agent(fig, cfg.posP, 'network 2', true,  cfg);
    draw_coupling(fig, infoA, infoP, cfg);
    cfg.head = infoA.head;
    fprintf('ARROWHEAD (all arrows): length %.2f pt, width %.2f pt, edge %.2f pt\n', ...
        cfg.head(1), 2*cfg.head(2), cfg.head(3));

    axR  = draw_readout(fig, dat, cfg);
    % Canvas position (pt) of B's two end dots, for C's dots and connector lines.
    Hz = dat.horizonEvolve; yl = cfg.bYlim; aB = get(axR, 'Position');
    cfg.dotYAbs = (aB(2) + aB(4)*([dat.r_cond(Hz) dat.r_uncond(Hz)] - yl(1))/diff(yl))*cfg.ptH;
    cfg.bDotXAbs = (aB(1) + aB(3))*cfg.ptW;

    % B's heading, centered over B, its ink cfg.headGapPt above B's top edge;
    % C's heading takes the same baseline.
    hB = text(axR, 0, 0, 'network 1 is conditional', 'FontSize', cfg.fsText, 'FontName', figure_utils.FONT, ...
        'Color', cfg.colNote, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline', ...
        'Units', 'points', 'Clipping', 'off', 'Interpreter', 'none');
    ibB = place_ink(hB, 1, 1);
    axPosB = get(axR, 'Position');
    hB.Position(1:2) = [axPosB(3)*cfg.ptW/2 - (ibB(1)+ibB(3))/2, axPosB(4)*cfg.ptH + cfg.headGapPt - ibB(2)];
    cfg.headBaseAbs = axPosB(2)*cfg.ptH + hB.Position(2);   % pt from the canvas bottom
    set(hB, 'Units', 'data');

    % The table: square, from gL to the right edge that leaves room for its
    % right-hand labels, its bottom edge the label depth above the top of B's
    % x title (so C's "network 1" stands level with B's x title).
    hb = axR.XLabel; set(hb, 'Units', 'points'); E = hb.Extent; ins = ink_insets(hb);
    bTop = cfg.posR(2)*cfg.ptH + E(2) + E(4) + ins(4); set(hb, 'Units', 'data');
    gBot = bTop + depth;
    gR = Lt.gridRight - Rr - 0.4;   % 0.4 pt: the rotated labels' boxes reach past their ink
    Wg = gR - gL; Hg = Wg;
    topPt = topEdge*cfg.ptH;
    if gBot + Hg > topPt, Hg = topPt - gBot; Wg = Hg; end
    cfg.posC = [gL/W, gBot/cfg.ptH, Wg/W, Hg/cfg.ptH];
    fprintf('TABLE PLACE: B x-title ink top %.2f pt; grid %.1f..%.1f x %.1f..%.1f pt\n', ...
        bTop, gL, gL + Wg, gBot, gBot + Hg);
    axC = draw_table(fig, cfg);
    draw_strip(fig, cfg);

    axLbl = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', ...
        'FontName', figure_utils.FONT, 'FontSize', cfg.fsText);
    xlim(axLbl, [0 1]); ylim(axLbl, [0 1]);
    panel_label(axLbl, cfg.posA(1) - 0.032, cfg.letterY, 'A', cfg);
    panel_label(axLbl, cfg.posR(1) - 0.082, cfg.letterY, 'B', cfg);   % clears B's tick labels and y title
    panel_label(axLbl, cfg.posC(1) - (Sst + 4)/W, cfg.letterY, 'C', cfg);
    uistack(axR, 'top');
    uistack(axC, 'top');

    report_fonts(fig, cfg);
end

% =========================================================================
function report_fonts(fig, cfg)
%REPORT_FONTS  Type inventory and the node fills.
    h = findall(fig, '-property', 'FontSize');
    sizes = []; fams = {}; bolds = 0;
    for k = 1:numel(h)
        if isprop(h(k), 'String') && isempty(get(h(k), 'String')) && ...
                ~isa(h(k), 'matlab.graphics.axis.Axes')
            continue;
        end
        sizes(end+1) = get(h(k), 'FontSize'); %#ok<AGROW>
        fams{end+1}  = get(h(k), 'FontName'); %#ok<AGROW>
        if isprop(h(k), 'FontWeight') && strcmp(get(h(k), 'FontWeight'), 'bold')
            bolds = bolds + 1;
        end
    end
    fprintf('FONT REPORT: families {%s}; sizes [%s] pt; bold objects: %d\n', ...
        strjoin(unique(fams), ', '), ...
        strjoin(arrayfun(@(v) sprintf('%.3g', v), unique(sizes), 'Uni', 0), ' '), bolds);
    fprintf('NODE FILLS: C = %s (cfg.coop), D = %s (cfg.defect); every gene circle rimmed %s at %.1f pt\n', ...
        hex(cfg.coop), hex(cfg.defect), hex(cfg.colAxis), cfg.lwSchem);
end

function s = hex(c)
    s = sprintf('#%02X%02X%02X', round(255*c(1)), round(255*c(2)), round(255*c(3)));
end

% =========================================================================
function panel_label(ax, x, y, s, cfg)
    text(ax, x, y, s, 'FontSize', cfg.fsPanel, 'FontWeight', 'bold', ...
        'FontName', figure_utils.FONT, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'Color', [0 0 0]);
end

% =========================================================================
function draw_strip(fig, cfg)
%DRAW_STRIP  The six stages of one generation in order, joined by arrows;
%   "development" is the only dark word, and a stick drops from just below it
%   to a bracket spanning panels A and B.

    ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', ...
        'FontName', figure_utils.FONT, 'FontSize', cfg.fsText);
    xlim(ax, [0 1]); ylim(ax, [0 1]); hold(ax, 'on');

    words  = {'random pairing', 'development', 'payoff', 'selection', ...
              'reproduction', 'mutation'};
    isDev  = [false true false false false false];
    pad    = 0.008;          % white either side of an arrow glyph

    % Lay the line out from measured extents, then center it.
    h = gobjects(1, 2*numel(words)-1);
    x = 0; k = 0;
    for w = 1:numel(words)
        if w > 1
            k = k + 1;
            % Unicode arrow set in Arial (TeX's \rightarrow would bring in a
            % second font family).
            h(k) = text(ax, x + pad, cfg.yWords, char(8594), 'FontSize', cfg.fsText, ...
                'FontName', figure_utils.FONT, 'Color', cfg.colFaint, 'Interpreter', 'none', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
            x = x + 2*pad + h(k).Extent(3);
        end
        k = k + 1;
        if isDev(w), col = cfg.colAxis; else, col = cfg.colFaint; end
        h(k) = text(ax, x, cfg.yWords, words{w}, 'FontSize', cfg.fsText, ...
            'FontName', figure_utils.FONT, 'Color', col, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        x = x + h(k).Extent(3);
    end
    total = x;
    x0 = 0.5 - total/2;
    devX = NaN;
    for k = 1:numel(h)
        h(k).Position(1) = h(k).Position(1) + x0;
        if strcmp(h(k).String, 'development')
            devX = h(k).Position(1) + h(k).Extent(3)/2;
            devBot = h(k).Extent(2);
        end
    end

    % The stick from just below "development" down to the A-B bracket.
    pL = cfg.posA(1);
    pR = cfg.posR(1) + cfg.posR(3);
    dropTick = 0.020;   % long enough that the span reads as a bracket
    stickTop = devBot - cfg.stickGapPt/cfg.ptH;
    plot(ax, [devX devX], [stickTop cfg.yDrop], '-', 'Color', cfg.colFaint, ...
        'LineWidth', cfg.lwRef);
    plot(ax, [pL pR], [cfg.yDrop cfg.yDrop], '-', 'Color', cfg.colFaint, 'LineWidth', cfg.lwRef);
    plot(ax, [pL pL], cfg.yDrop + [0 -dropTick], '-', 'Color', cfg.colFaint, 'LineWidth', cfg.lwRef);
    plot(ax, [pR pR], cfg.yDrop + [0 -dropTick], '-', 'Color', cfg.colFaint, 'LineWidth', cfg.lwRef);
end

% =========================================================================
function info = draw_agent(fig, pos, netLabel, mirrored, cfg)
%DRAW_AGENT  One 4-gene network, sensors -> wiring -> outputs. Every gene
%   circle carries the same dark rim; the output genes are filled with the
%   cooperation and defection colors, their letters white. Returns node
%   centers and radii in FIGURE coordinates and the arrowhead size in pt.

    ax = axes(fig, 'Position', pos, 'FontName', figure_utils.FONT, 'FontSize', cfg.fsText);
    hold(ax, 'on'); axis(ax, 'off');
    xlim(ax, [-0.24 1.24]); ylim(ax, [-0.34 1.34]);

    if ~mirrored
        xSens = 0.0; xProc = 0.5; xOut = 1.0;
    else
        xSens = 1.0; xProc = 0.5; xOut = 0.0;
    end
    yTop = 0.74; yBot = 0.26;

    figPos = get(fig, 'Position');
    xAsp = (pos(4)*figPos(4)/diff(ylim(ax))) / (pos(3)*figPos(3)/diff(xlim(ax)));
    r = 0.155; ry = r; rx = r * xAsp;
    theta = linspace(0, 2*pi, 80);

    ex = 0.19; ey = 0.30;
    fill(ax, xProc + ex*cos(theta), 0.5 + ey*sin(theta), [0.96 0.96 0.96], ...
        'EdgeColor', [0.55 0.55 0.55], 'LineWidth', cfg.lwRef, 'LineStyle', '--');
    text(ax, xProc, 0.5, 'wiring', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', cfg.fsText, 'FontName', figure_utils.FONT, ...
        'Color', cfg.colNote);

    node(ax, xSens, yTop, rx, ry, [1 1 1],    'S_{\fontsize{7.5}C}', cfg.colAxis, cfg.fsText,     cfg);
    node(ax, xSens, yBot, rx, ry, [1 1 1],    'S_{\fontsize{7.5}D}', cfg.colAxis, cfg.fsText,     cfg);
    node(ax, xOut,  yTop, rx, ry, cfg.coop,   'C',   [1 1 1],     cfg.fsNodeFill, cfg);
    node(ax, xOut,  yBot, rx, ry, cfg.defect, 'D',   [1 1 1],     cfg.fsNodeD,    cfg);

    % Each arrow runs from a node center towards an aim point on the wiring
    % block (or back) and stops cfg.gapPt short of the outer edge of the
    % stroke it meets at both ends: the node rim (cfg.lwSchem) and the
    % dashed block outline (cfg.lwRef); at a tip the arrowhead's overshoot
    % (cfg.headTip) is added. Distances in pt.
    sgn = 1; if mirrored, sgn = -1; end
    ptX = pos(3)*cfg.ptW/diff(xlim(ax));     % pt per data unit, x
    ptY = pos(4)*cfg.ptH/diff(ylim(ax));     % pt per data unit, y
    gNodeTail = cfg.gapPt + cfg.lwSchem/2;
    gNodeTip  = cfg.gapPt + cfg.lwSchem/2 + cfg.headTip;
    gBlkTail  = cfg.gapPt + cfg.lwRef/2;
    gBlkTip   = cfg.gapPt + cfg.lwRef/2 + cfg.headTip;
    E = [xProc 0.5 ex ey];
    hd = zeros(0, 3);   % per arrow: head length, half-width (pt), edge width
    for yy = [yTop 0.62; yBot 0.38]'
        % sensor -> wiring
        C = [xSens yy(1)]; T = [xProc - sgn*ex, yy(2)];
        [px, py] = node_edge(C(1), C(2), T(1), T(2), rx, ry);
        P = step_pt([px py], T, gNodeTail, ptX, ptY);
        Q = step_pt(ray_ellipse(C, T, E), C, gBlkTip, ptX, ptY);
        hd(end+1, :) = flow(ax, P(1), P(2), Q(1), Q(2), cfg, ptX, ptY); %#ok<AGROW>
        % wiring -> output
        O = [xOut yy(1)]; S = [xProc + sgn*ex, yy(2)];
        [ox, oy] = node_edge(O(1), O(2), S(1), S(2), rx, ry);
        P = step_pt(ray_ellipse(O, S, E), O, gBlkTail, ptX, ptY);
        Q = step_pt([ox oy], S, gNodeTip, ptX, ptY);
        hd(end+1, :) = flow(ax, P(1), P(2), Q(1), Q(2), cfg, ptX, ptY); %#ok<AGROW>
    end

    text(ax, 0.5, 1.24, netLabel, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', cfg.fsText, ...
        'FontName', figure_utils.FONT, 'Color', cfg.colNote);

    xl = xlim(ax); yl = ylim(ax);
    toFig = @(xd, yd) [pos(1) + (xd - xl(1))/diff(xl)*pos(3), ...
                       pos(2) + (yd - yl(1))/diff(yl)*pos(4)];
    info.sensC = toFig(xSens, yTop);
    info.sensD = toFig(xSens, yBot);
    info.outC  = toFig(xOut,  yTop);
    info.outD  = toFig(xOut,  yBot);
    info.rx    = rx/diff(xl)*pos(3);
    info.ry    = ry/diff(yl)*pos(4);
    info.outSign = sign(xOut - 0.5);
    info.head = mean(hd, 1);   % the arrowhead every other arrow is drawn to
end

function node(ax, cx, cy, rx, ry, fillCol, name, letterCol, letterFs, cfg)
    theta = linspace(0, 2*pi, 80);
    fill(ax, cx + rx*cos(theta), cy + ry*sin(theta), fillCol, ...
        'EdgeColor', cfg.colAxis, 'LineWidth', cfg.lwSchem);
    text(ax, cx, cy, name, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', letterFs, ...
        'FontName', figure_utils.FONT, 'Color', letterCol);
end

function [px, py] = node_edge(cx, cy, tx, ty, rx, ry)
    dx = tx - cx; dy = ty - cy;
    s = 1 / sqrt((dx/rx)^2 + (dy/ry)^2);
    px = cx + s*dx; py = cy + s*dy;
end

function P = step_pt(P, towards, d, ptX, ptY)
%STEP_PT  Move P by d pt towards the point `towards` (data units scaled
%   ptX, ptY pt per unit in x, y).
    v = [(towards(1) - P(1))*ptX, (towards(2) - P(2))*ptY];
    v = v / norm(v);
    P = [P(1) + d*v(1)/ptX, P(2) + d*v(2)/ptY];
end

function X = ray_ellipse(C, T, E)
%RAY_ELLIPSE  First point where the ray from C (outside the ellipse) through
%   T meets the ellipse E = [cx cy ax ay].
    d = T - C; c0 = C - E(1:2);
    a = (d(1)/E(3))^2 + (d(2)/E(4))^2;
    b = 2*(c0(1)*d(1)/E(3)^2 + c0(2)*d(2)/E(4)^2);
    c = (c0(1)/E(3))^2 + (c0(2)/E(4))^2 - 1;
    t = (-b - sqrt(b^2 - 4*a*c)) / (2*a);
    X = C + t*d;
end

function hd = flow(ax, x1, y1, x2, y2, cfg, ptX, ptY)
%FLOW  Panel A's arrow inside a network: a cfg.lwSchem line with a filled
%   triangular head. Returns the head's length and half-width in pt and its
%   edge width, which every other arrow in the figure copies.
    col = cfg.colNote;
    plot(ax, [x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', cfg.lwSchem);
    dx = x2 - x1; dy = y2 - y1; d = sqrt(dx^2 + dy^2);
    ux = dx/d; uy = dy/d; aLen = 0.055; aW = 0.032;
    hp = fill(ax, [x2, x2 - aLen*ux + aW*uy, x2 - aLen*ux - aW*uy], ...
             [y2, y2 - aLen*uy - aW*ux, y2 - aLen*uy + aW*ux], col, 'EdgeColor', col);
    hd = [hypot(aLen*ux*ptX, aLen*uy*ptY), hypot(aW*uy*ptX, aW*ux*ptY), hp.LineWidth];
end

function arrow_pt(ax, P, Q, col, cfg, head, ptX, ptY)
%ARROW_PT  The figure's arrow anywhere: a cfg.lwSchem line from P to Q (data
%   units of ax, ptX x ptY pt per unit) with head `head` at Q, all in col.
    plot(ax, [P(1) Q(1)], [P(2) Q(2)], '-', 'Color', col, 'LineWidth', cfg.lwSchem);
    head_at(ax, Q, [(Q(1)-P(1))*ptX, (Q(2)-P(2))*ptY], col, head, ptX, ptY);
end

function head_at(ax, Q, v, col, head, ptX, ptY)
%HEAD_AT  Triangle head with its point at Q, pointing along v (a pt vector):
%   head(1) pt long, 2*head(2) pt wide, edge head(3) pt.
    v = v / norm(v); n = [-v(2) v(1)];
    b = -head(1)*v;                                   % tip -> base center, pt
    xs = Q(1) + [0, (b(1) + head(2)*n(1))/ptX, (b(1) - head(2)*n(1))/ptX];
    ys = Q(2) + [0, (b(2) + head(2)*n(2))/ptY, (b(2) - head(2)*n(2))/ptY];
    fill(ax, xs, ys, col, 'EdgeColor', col, 'LineWidth', head(3));
end

% =========================================================================
function draw_coupling(fig, infoA, infoP, cfg)
%DRAW_COUPLING  Each network's two output genes bundle, sweep wide, and fan
%   onto BOTH of the other's sensor genes: under contrast sensing each sensor
%   receives the signed difference of the other's C and D, so neither
%   sensor has a private channel (simulate_pair). Drawn on one invisible
%   figure-wide axes (data = normalized figure units).
    axO = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', 'Color', 'none', ...
        'FontName', figure_utils.FONT, 'FontSize', cfg.fsText);   % no text; font and size set so the type inventory counts them
    xlim(axO, [0 1]); ylim(axO, [0 1]); hold(axO, 'on');
    one_way(axO, infoA, infoP, cfg);
    one_way(axO, infoP, infoA, cfg);
end

function one_way(axO, src, dst, cfg)
    col = cfg.colNote;
    s   = src.outSign;
    gap = src.rx + 0.016;
    J = [mean([src.outC(1) src.outD(1)]) + s*gap, mean([src.outC(2) src.outD(2)])];
    K = [mean([dst.sensC(1) dst.sensD(1)]) + s*gap, mean([dst.sensC(2) dst.sensD(2)])];

    % The bundle starts cfg.gapPt off the output rims and the two arrows
    % stop cfg.gapPt short of the sensor rims (ptW x ptH pt per unit).
    g = cfg.gapPt + cfg.lwSchem/2;
    for nd = {src.outC, src.outD}
        p = step_pt(rim(nd{1}, src.rx, src.ry, J), J, g, cfg.ptW, cfg.ptH);
        plot(axO, [p(1) J(1)], [p(2) J(2)], '-', 'Color', col, 'LineWidth', cfg.lwSchem);
    end

    ctrl = [max([J(1) K(1)]*s)*s + s*0.055, mean([J(2) K(2)])];
    K(1) = K(1) + s*cfg.fanPt/cfg.ptW;      % fan point, cfg.fanPt further out
    tt = linspace(0, 1, 60)';
    bx = (1-tt).^2*J(1) + 2*(1-tt).*tt*ctrl(1) + tt.^2*K(1);
    by = (1-tt).^2*J(2) + 2*(1-tt).*tt*ctrl(2) + tt.^2*K(2);
    plot(axO, bx, by, '-', 'Color', col, 'LineWidth', cfg.lwSchem);

    for nd = {dst.sensC, dst.sensD}
        p = step_pt(rim(nd{1}, dst.rx, dst.ry, K), K, g + cfg.headTip, cfg.ptW, cfg.ptH);
        arrow_pt(axO, K, p, col, cfg, src.head, cfg.ptW, cfg.ptH);
    end
end

function p = rim(c, rx, ry, towards)
    dx = towards(1) - c(1); dy = towards(2) - c(2);
    s = 1 / sqrt((dx/rx)^2 + (dy/ry)^2);
    p = [c(1) + s*dx, c(2) + s*dy];
end

% =========================================================================

% =========================================================================
function ax = draw_readout(fig, dat, cfg)
%DRAW_READOUT  Network 1's readout x_C - x_D against developmental time, over
%   the horizon the population evolved under. Line style says what network 2
%   is (solid: conditional; dashed: an unconditional defector), named in a
%   boxless legend in the lower left under the heading "network 2 is"; both
%   curves are ink. Each ends in a dot colored by network 1's action (blue
%   cooperates, red-orange defects); height says the same: the dotted line is
%   labelled "action threshold" at its right end, with "cooperates" (blue)
%   one text line above it and "defects" (red-orange) one text line below
%   it, the key to the two dot colors.

    H = dat.horizonEvolve;
    t = dat.t(1:H);
    rC = dat.r_cond(1:H); rU = dat.r_uncond(1:H);

    ax = axes(fig, 'Position', cfg.posR);
    hold(ax, 'on');

    hZero = yline(ax, 0, ':');
    hZero.Color = cfg.colRef;
    hZero.LineWidth = cfg.lwRef;
    hZero.FontName = figure_utils.FONT;
    hZero.FontSize = cfg.fsText;

    set(ax, 'XScale', 'log');
    xlim(ax, [1 H]);
    set(ax, 'XTick', [1 3 10 30 100 300]);
    set(ax, 'XTickLabelRotation', 0);
    ylim(ax, cfg.bYlim);
    set(ax, 'YTick', cfg.bYTicks);

    % Curves and legend are laid out in pt from the axes' lower-left corner.
    pos = get(ax, 'Position');
    axW = pos(3)*cfg.ptW; axH = pos(4)*cfg.ptH;
    xl = log10(xlim(ax)); yl = ylim(ax);
    toData = @(xp, yp) [10.^(xl(1) + xp/axW*diff(xl)), yl(1) + yp/axH*diff(yl)];
    toPt   = @(xd, yd) [axW*(log10(xd) - xl(1))/diff(xl), axH*(yd - yl(1))/diff(yl)];

    plot(ax, t, rC, '-', 'Color', cfg.ink, 'LineWidth', cfg.lwData);
    dashed_pt(ax, toPt(t, rU), toData, cfg.dashPt, cfg.ink, cfg.lwData);

    plot(ax, H, rC(H), 'o', 'MarkerSize', cfg.mkSize, ...
        'MarkerFaceColor', cfg.coop,   'MarkerEdgeColor', 'w', 'LineWidth', cfg.mkEdge);
    plot(ax, H, rU(H), 'o', 'MarkerSize', cfg.mkSize, ...
        'MarkerFaceColor', cfg.defect, 'MarkerEdgeColor', 'w', 'LineWidth', cfg.mkEdge);
    xlabel(ax, 'developmental time (steps)', 'FontSize', cfg.fsAxis);
    hY = ylabel(ax, 'network 1''s readout {\itx}_{\fontsize{7.5}C} - {\itx}_{\fontsize{7.5}D}', 'FontSize', cfg.fsAxis);
    figure_utils.style_axis(ax);
    set(ax, 'FontSize', cfg.fsTick, 'LineWidth', cfg.lwAxis, ...
        'XColor', cfg.colAxis, 'YColor', cfg.colAxis);
    ax.XLabel.FontSize = cfg.fsAxis;
    ax.YLabel.FontSize = cfg.fsAxis;
    report_ylabel_room(ax, hY, cfg, 'B');

    tp = {'FontSize', cfg.fsText, 'FontName', figure_utils.FONT, 'Units', 'points'};

    % The key sits at the right end, clear of the end dots: "action threshold"
    % just above the dotted line, "cooperates" one text line above it,
    % "defects" one text line below the line (ink cfg.keyGapPt from it).
    y0 = toPt(1, 0); y0 = y0(2);
    xR = axW - cfg.sideRight;                        % right edge of the words, pt
    asc = cfg.inkAsc*cfg.fsText;
    lineBot = y0 - cfg.lwRef/2; lineTop = y0 + cfg.lwRef/2;
    g = cfg.keyGapPt;
    hWords = [ ...
        text(ax, xR, lineTop + g, 'action threshold', tp{:}, 'Color', cfg.colNote, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline'), ...
        text(ax, xR, lineTop + g + cfg.keyLinePt, 'cooperates', tp{:}, 'Color', cfg.coop, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline'), ...
        text(ax, xR, lineBot - g - asc, 'defects', tp{:}, 'Color', cfg.defect, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline')];

    % The legend, built bottom-up from one text line's height; each sample
    % sits on the first line of its row, and the heading row has no sample
    % and starts at the legend's left edge.
    rows = cfg.legRows;
    probe = text(ax, 0, 0, 'Xg', tp{:}); lh = probe.Extent(4); delete(probe);
    nLines = sum(cellfun(@(r) numel(r{2}), rows));
    yTopRow = cfg.legBottom + nLines*lh;            % top of the legend block, pt
    hLeg = gobjects(0); boxes = zeros(0, 4);
    yCur = yTopRow;
    for r = 1:numel(rows)
        st = rows{r}{1}; lines = rows{r}{2};
        yMid = yCur - lh/2;
        textX = cfg.legX + cfg.legSample + cfg.legGap;
        if isempty(st)
            textX = cfg.legX;
        elseif strcmp(st, '--')
            dashed_pt(ax, [cfg.legX yMid; cfg.legX + cfg.legSample yMid], toData, ...
                cfg.dashPt, cfg.ink, cfg.lwData);
        else
            P = toData(cfg.legX, yMid); Q = toData(cfg.legX + cfg.legSample, yMid);
            plot(ax, [P(1) Q(1)], [P(2) Q(2)], '-', 'Color', cfg.ink, 'LineWidth', cfg.lwData);
        end
        if ~isempty(st)
            boxes(end+1, :) = [cfg.legX, yMid - cfg.lwData/2, cfg.legX + cfg.legSample, yMid + cfg.lwData/2]; %#ok<AGROW>
        end
        for k = 1:numel(lines)
            h = text(ax, textX, yCur - (k - 0.5)*lh, lines{k}, tp{:}, ...
                'Color', cfg.colAxis, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
            hLeg(end+1) = h; %#ok<AGROW>
        end
        yCur = yCur - numel(lines)*lh;
    end
    report_text_clearance([hWords hLeg], boxes, toPt, t, {rC, rU}, ...
        [H rC(H); H rU(H)], cfg);
    set([hWords hLeg], 'Units', 'data');
end

function dashed_pt(ax, Ppt, toData, pattern, col, lw)
%DASHED_PT  Draw the polyline Ppt (pt, axes frame) dashed by arc length:
%   pattern = [dash gap] in pt, starting with a dash, so the dashes keep one
%   printed length however steep the curve. One line object, NaN-separated.
    seg = hypot(diff(Ppt(:,1)), diff(Ppt(:,2)));
    s = [0; cumsum(seg)];
    period = sum(pattern);
    X = []; Y = [];
    for a = 0:period:s(end)
        b = min(a + pattern(1), s(end));
        ss = [a; s(s > a & s < b); b];
        P = [interp1(s, Ppt(:,1), ss), interp1(s, Ppt(:,2), ss)];
        D = toData(P(:,1), P(:,2));
        X = [X; D(:,1); NaN]; Y = [Y; D(:,2); NaN]; %#ok<AGROW>
    end
    plot(ax, X, Y, '-', 'Color', col, 'LineWidth', lw);
end

function report_text_clearance(hText, lineBoxes, toPt, t, curves, dots, cfg)
%REPORT_TEXT_CLEARANCE  Print, in pt, the smallest gap from each panel-B text
%   box and legend line sample to the two curves (stroke edge), the two end
%   dots (outer edge) and the axes' left and bottom lines. Text boxes carry
%   padding, so a small negative value can still mean clear ink.
    tf = exp(linspace(log(t(1)), log(t(end)), 6000))';
    C = zeros(0, 2);
    for k = 1:numel(curves), C = [C; toPt(tf, interp1(t, curves{k}, tf))]; end %#ok<AGROW>
    Dp = toPt(dots(:,1), dots(:,2));
    rDot = cfg.mkSize/2 + cfg.mkEdge/2;
    bd = @(P, b) min(hypot(max([b(1) - P(:,1), zeros(size(P,1),1), P(:,1) - b(3)], [], 2), ...
                           max([b(2) - P(:,2), zeros(size(P,1),1), P(:,2) - b(4)], [], 2)));
    items = {};
    for k = 1:numel(hText)
        E = get(hText(k), 'Extent'); items(end+1, :) = {hText(k).String, [E(1) E(2) E(1)+E(3) E(2)+E(4)]}; %#ok<AGROW>
    end
    for k = 1:size(lineBoxes, 1), items(end+1, :) = {sprintf('legend sample %d', k), lineBoxes(k, :)}; end %#ok<AGROW>
    for k = 1:size(items, 1)
        b = items{k, 2};
        fprintf('B TEXT "%s": to curves %.2f pt, to end dots %.2f pt, to axes %.2f pt\n', ...
            items{k, 1}, bd(C, b) - cfg.lwData/2, bd(Dp, b) - rDot, min(b(1), b(2)));
    end
end

function report_ylabel_room(ax, hY, cfg, name)
%REPORT_YLABEL_ROOM  Print how far the y-axis title's text box stops short of
%   the axis ends (pt; negative = runs past).
    set(hY, 'Units', 'points'); E = get(hY, 'Extent'); set(hY, 'Units', 'data');
    pos = get(ax, 'Position'); axH = pos(4)*cfg.ptH;
    fprintf('YLABEL %s: box spans %.2f..%.2f pt along an axis of %.2f pt; room at bottom %.2f pt, at top %.2f pt\n', ...
        name, E(2), E(2) + E(4), axH, E(2), axH - (E(2) + E(4)));
end


function g = box_gap(a, b)
    dx = max([b(1) - a(3), a(1) - b(3), 0]);
    dy = max([b(2) - a(4), a(2) - b(4), 0]);
    if dx == 0 && dy == 0
        g = -min(min(a(3), b(3)) - max(a(1), b(1)), min(a(4), b(4)) - max(a(2), b(2)));
    else
        g = hypot(dx, dy);
    end
end

function ib = place_ink(h, ptX, ptY)
%PLACE_INK  Ink box [x1 y1 x2 y2] (pt) of text h placed at (0,0): its
%   Extent plus the ink insets measured on a 1200 dpi render of the string.
    h.Position(1:2) = [0 0];
    E = h.Extent;
    Ep = [E(1)*ptX, E(2)*ptY, (E(1)+E(3))*ptX, (E(2)+E(4))*ptY];
    ins = ink_insets(h);
    ib = Ep + ins;
end

function ins = ink_insets(h)
%INK_INSETS  [left bottom right top] offsets (pt) from a text's Extent box
%   to its ink, from a render of the same string alone.
    persistent cache
    if isempty(cache), cache = containers.Map(); end
    key = sprintf('%s|%s|%s|%g|%g', strjoin(cellstr(h.String), '\n'), h.HorizontalAlignment, ...
        h.VerticalAlignment, h.FontSize, h.Rotation);
    if isKey(cache, key), ins = cache(key); return; end
    W = 300; Hh = 150; res = 1200; s = res/72;
    tf = figure('Visible', 'off', 'Units', 'points', 'Position', [0 0 W Hh], 'Color', 'w', ...
        'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
    ta = axes(tf, 'Position', [0 0 1 1], 'Units', 'normalized'); axis(ta, 'off');
    xlim(ta, [0 W]); ylim(ta, [0 Hh]);
    t = text(ta, W/2, Hh/2, h.String, 'FontSize', h.FontSize, 'FontName', h.FontName, ...
        'HorizontalAlignment', h.HorizontalAlignment, 'VerticalAlignment', h.VerticalAlignment, ...
        'Interpreter', h.Interpreter, 'Rotation', h.Rotation, 'Color', [0 0 0]);
    E = t.Extent;
    img = print(tf, '-RGBImage', sprintf('-r%d', res));
    close(tf);
    sx = size(img, 2)/W; sy = size(img, 1)/Hh;
    m = any(img < 200, 3);
    [rr, cc] = find(m);
    ink = [(min(cc)-1)/sx, Hh - max(rr)/sy, max(cc)/sx, Hh - (min(rr)-1)/sy];
    ins = ink - [E(1), E(2), E(1)+E(3), E(2)+E(4)];
    cache(key) = ins;
    if abs(sx - s) > 0.05*s, warning('fig1_encounter_pair:scale', 'render scale %.2f vs %.2f', sx, s); end
end

% =========================================================================
function L = table_layout(runCfg)
%TABLE_LAYOUT  Kinds, cell texts, dots and spacing of panel C.
%   The incentive to defect on a network X = (the other network's payoff
%   defecting on X) - (its payoff cooperating with X): P - R on a
%   conditional network, which defects back; P - S on an unconditional
%   defector and T - R on an unconditional cooperator. The cell texts say a
%   network gains by defecting exactly where its partner is unconditional,
%   which needs the first negative and the other two positive; the run's
%   payoffs must give -2, +1 and +2 (error otherwise).
    T = runCfg.payoffDC; R = runCfg.payoffCC; P = runCfg.payoffDD; S = runCfg.payoffCD;
    incCond = P - R; incDef = P - S; incCoop = T - R;
    fprintf('TABLE payoffs T, R, P, S = %g, %g, %g, %g; incentive to defect on: a conditional network %+g, an unconditional one %+g (defects) or %+g (cooperates)\n', ...
        T, R, P, S, incCond, incDef, incCoop);
    if incCond ~= -2 || incDef ~= 1 || incCoop ~= 2
        error('fig1_encounter_pair:payoffs', 'payoffs give %g, %g, %g, not -2, +1, +2', incCond, incDef, incCoop);
    end
    L.inc = [incCond incDef incCoop];
    % Cells are indexed [column row]: column 1 = network 1 conditional,
    % column 2 = unconditional; row 1 (bottom) = network 2 unconditional,
    % row 2 (top) = network 2 conditional.
    L.kinds    = {'conditional', 'unconditional'};   % columns, left to right
    L.rowKinds = {'unconditional', 'conditional'};   % rows, bottom to top
    L.xTitle = 'network 1'; L.yTitle = 'network 2';
    L.gap = 3;          % ink gap from the grid to its labels (pt)
    L.titleSep = 5;     % ink gap from the labels to "network 1" / "network 2" (pt)
    L.cellInset = 4;    % ink gap from a cell's text block to its grid line (pt)
    dt = @(cellCR, face) struct('cell', cellCR, 'face', face);
    % Blue: both conditional (B's solid meeting); red-orange: network 1
    % conditional, network 2 unconditional (B's dashed meeting). Same order
    % as B's end dots (r_cond, r_uncond).
    L.dots = [dt([1 2], figure_utils.COLOR_COOP), dt([1 1], figure_utils.COLOR_DEFECT)];
    ct = @(cellCR, lines, nReason) struct('cell', cellCR, 'lines', {lines}, 'nReason', nReason);
    L.cellText = [ct([2 1], {'each gains', 'by defecting', 'dilemma', 'holds'}, 2), ...
                  ct([1 2], {'neither gains', 'by defecting', 'dilemma', 'escaped'}, 2)];
    L.whiteText = struct('cell', {[2 2], [1 1]}, ...
                         'lines', {{'network 2 gains', 'by defecting'}, {'network 1 gains', 'by defecting'}});
    L.escCell = [1 2]; L.holdsCell = [2 1];
    L.shiftPt   = 18;     % A and B move this far left of their home positions (pt)
    L.bShiftPt  = 16.2;   % B then moves this far right, evening the A-B and B-C gaps (pt)
    L.abGap     = 27.7;   % B's end dots (outer edge) to the table's left edge (pt)
    L.gridRight = 500;    % right edge of C's right-hand labels (pt on the canvas)
end

% =========================================================================
function [ax, S, depth, Rr] = draw_table(fig, cfg, measureOnly)
%DRAW_TABLE  Panel C as a two-by-two table: columns = what network 1 is,
%   rows = what network 2 is; categorical labels, no numeric axis. The
%   axes' data units are pt (origin at the grid's lower-left corner).
%   Returns S, the ink reach of the labels left of the grid (pt, >= 0), the
%   depth of the labels below it and Rr, their reach right of it. With
%   measureOnly, draws only the grid and its labels and returns after
%   measuring them.
    if nargin < 3, measureOnly = false; end
    L = cfg.table;
    pos = cfg.posC;
    Wg = pos(3)*cfg.ptW; Hg = pos(4)*cfg.ptH;
    ax = axes(fig, 'Position', pos, 'Visible', 'off', 'FontName', figure_utils.FONT, 'FontSize', cfg.fsText);
    hold(ax, 'on'); xlim(ax, [0 Wg]); ylim(ax, [0 Hg]);
    cw = Wg/2; ch = Hg/2;
    cb = @(c) [(c(1)-1)*cw, c(1)*cw, c(1)*cw, (c(1)-1)*cw; (c(2)-1)*ch, (c(2)-1)*ch, c(2)*ch, c(2)*ch];
    xyE = cb(L.escCell); xyH = cb(L.holdsCell);
    hFill = [fill(ax, xyE(1,:), xyE(2,:), cfg.tintEscaped, 'EdgeColor', 'none'), ...
             fill(ax, xyH(1,:), xyH(2,:), cfg.tintHolds, 'EdgeColor', 'none')];
    for k = 0:2
        plot(ax, [k*cw k*cw], [0 Hg], '-', 'Color', cfg.colAxis, 'LineWidth', cfg.lwAxis);
        plot(ax, [0 Wg], [k*ch k*ch], '-', 'Color', cfg.colAxis, 'LineWidth', cfg.lwAxis);
    end

    mk = @(str, col, fsz, ha, va, rot) text(ax, 0, 0, str, 'FontSize', fsz, 'FontName', figure_utils.FONT, ...
        'Color', col, 'HorizontalAlignment', ha, 'VerticalAlignment', va, 'Rotation', rot, ...
        'Clipping', 'off', 'Interpreter', 'none');
    hT = gobjects(0); inkB = zeros(0, 4); names = {}; cellOf = zeros(0, 2);
    function add(h, ib, d, nm, cc)
        set(h, 'Position', [d(1), d(2), 0]);
        hT(end+1) = h; inkB(end+1, :) = ib + [d d]; names{end+1} = nm; cellOf(end+1, :) = cc;
    end
    G = L.gap;
    probe = mk({'Xg', 'Xg'}, cfg.colNote, cfg.fsText, 'center', 'baseline', 0);
    e2 = probe.Extent(4); set(probe, 'String', 'Xg'); pitch = e2 - probe.Extent(4); delete(probe);

    % ---- below the grid: network 1's kinds on one baseline, then its title ----
    base1 = NaN;
    for c = 1:2
        h = mk(L.kinds{c}, cfg.colAxis, cfg.fsText, 'center', 'baseline', 0);
        ib = place_ink(h, 1, 1);
        if isnan(base1), base1 = -G - cfg.lwAxis/2 - ib(4); end
        add(h, ib, [(c - 0.5)*cw - (ib(1)+ib(3))/2, base1], ['col head ' L.kinds{c}], [0 0]);
    end
    lowest = min(inkB(:, 2));
    % ---- right of the grid: network 2's kinds, rotated, then its title ----
    for r = 1:2
        h = mk(L.rowKinds{r}, cfg.colAxis, cfg.fsText, 'center', 'bottom', 90);
        ib = place_ink(h, 1, 1);
        add(h, ib, [Wg + G + cfg.lwAxis/2 - ib(1), (r - 0.5)*ch - (ib(2)+ib(4))/2], ['row head ' L.rowKinds{r}], [0 0]);
    end
    hx = mk(L.xTitle, cfg.colAxis, cfg.fsAxis, 'center', 'top', 0);
    ib = place_ink(hx, 1, 1);
    add(hx, ib, [Wg/2 - (ib(1)+ib(3))/2, lowest - L.titleSep - ib(4)], 'x title', [0 0]);
    depth = -inkB(end, 4);
    hy = mk(L.yTitle, cfg.colAxis, cfg.fsAxis, 'center', 'bottom', 90);
    ib = place_ink(hy, 1, 1);
    rowsRight = max(inkB(cellfun(@(n) strncmp(n, 'row head', 8), names), 3));
    add(hy, ib, [rowsRight + L.titleSep - ib(1), Hg/2 - (ib(2)+ib(4))/2], 'y title', [0 0]);
    S = max(-min(inkB(:, 1)), 0);
    Rr = max(0, max(inkB(:, 3)) - Wg);
    if measureOnly, return; end

    % ---- heading, on the baseline of B's heading ----
    h = mk('in one encounter', cfg.colNote, cfg.fsText, 'center', 'baseline', 0);
    ib = place_ink(h, 1, 1);
    add(h, ib, [Wg/2 - (ib(1)+ib(3))/2, cfg.headBaseAbs - pos(2)*cfg.ptH], 'heading', [0 0]);

    % ---- dots, level with B's end dots ----
    rDot = cfg.mkSize/2 + cfg.mkEdge/2;
    dotsXY = zeros(numel(L.dots), 2);
    for k = 1:numel(L.dots)
        dotsXY(k, :) = [(L.dots(k).cell(1) - 0.5)*cw, cfg.dotYAbs(k) - pos(2)*cfg.ptH];
        plot(ax, dotsXY(k, 1), dotsXY(k, 2), 'o', 'MarkerSize', cfg.mkSize, ...
            'MarkerFaceColor', L.dots(k).face, 'MarkerEdgeColor', 'w', 'LineWidth', cfg.mkEdge);
    end

    % ---- filled cells: two grey lines (who gains) over two ink lines (the
    % dilemma). A top-row block stands on its cell's bottom edge, clear of
    % the dot above; a bottom-row block hangs from its cell's top edge,
    % clear of the dot below. ----
    rowLine1 = [NaN NaN];
    for k = 1:numel(L.cellText)
        T = L.cellText(k);
        n = numel(T.lines);
        cols = [repmat({cfg.colNote}, 1, T.nReason), repmat({cfg.ink}, 1, n - T.nReason)];
        hs = gobjects(1, n); ibs = zeros(n, 4);
        for j = 1:n
            hs(j) = mk(T.lines{j}, cols{j}, cfg.fsText, 'center', 'baseline', 0);
            ibs(j, :) = place_ink(hs(j), 1, 1);
        end
        if T.cell(2) == 1
            base = T.cell(2)*ch - cfg.lwAxis/2 - L.cellInset - ibs(1, 4) - (n - 1)*pitch;
        else
            base = (T.cell(2) - 1)*ch + cfg.lwAxis/2 + L.cellInset - ibs(n, 2);
        end
        rowLine1(T.cell(2)) = base + (n - 1)*pitch;   % baseline of the block's first line
        for j = 1:n
            add(hs(j), ibs(j, :), [(T.cell(1) - 0.5)*cw - (ibs(j,1)+ibs(j,3))/2, base + (n - j)*pitch], ...
                sprintf('cell %d%d line %d', T.cell, j), T.cell);
        end
    end

    % ---- white cells: two grey lines, first line on the same row's first line ----
    for k = 1:numel(L.whiteText)
        T = L.whiteText(k); n = numel(T.lines);
        for j = 1:n
            hw = mk(T.lines{j}, cfg.colNote, cfg.fsText, 'center', 'baseline', 0);
            ibw = place_ink(hw, 1, 1);
            add(hw, ibw, [(T.cell(1) - 0.5)*cw - (ibw(1)+ibw(3))/2, rowLine1(T.cell(2)) - (j - 1)*pitch], ...
                sprintf('white cell %d%d line %d', T.cell, j), T.cell);
        end
    end

    % ---- faint level lines from B's end dots to the table's dots, broken
    % round any text they would cross; drawn above the fills, below the rest ----
    hL = gobjects(0);
    for k = 1:numel(L.dots)
        yk = dotsXY(k, 2);
        xa = cfg.bDotXAbs + rDot - pos(1)*cfg.ptW;   % outer edge of B's dot
        xb = dotsXY(k, 1) - rDot;                    % outer edge of the table's dot
        segs = [xa xb];
        for i = 1:size(inkB, 1)
            if inkB(i,2) - 1 < yk && yk < inkB(i,4) + 1 && inkB(i,3) + 2 > xa && inkB(i,1) - 2 < xb
                lo = inkB(i,1) - 2.5; hi = inkB(i,3) + 2.5; new = zeros(0, 2);
                for q = 1:size(segs, 1)
                    if hi <= segs(q,1) || lo >= segs(q,2), new(end+1, :) = segs(q, :); continue; end %#ok<AGROW>
                    if lo > segs(q,1), new(end+1, :) = [segs(q,1) lo]; end %#ok<AGROW>
                    if hi < segs(q,2), new(end+1, :) = [hi segs(q,2)]; end %#ok<AGROW>
                end
                segs = new;
            end
        end
        for q = 1:size(segs, 1)
            hL(end+1) = plot(ax, segs(q, :), [yk yk], '-', 'Color', cfg.colLink, ...
                'LineWidth', cfg.lwLink, 'Clipping', 'off'); %#ok<AGROW>
        end
        fprintf('CONNECTOR %d: level y = %.2f pt above the table bottom (canvas %.2f pt), %d piece(s)\n', ...
            k, yk, cfg.dotYAbs(k), size(segs, 1));
    end
    kids = ax.Children; kids = kids(~ismember(kids, [hL hFill]));
    ax.Children = [kids; hL(:); hFill(:)];

    % ---- ink gaps (pt) ----
    gmin = inf;
    for i = 1:size(inkB, 1)
        a = inkB(i, :);
        gT = inf;
        for j = setdiff(1:size(inkB, 1), i)
            if any(cellOf(i,:)) && isequal(cellOf(i,:), cellOf(j,:)), continue; end   % lines of one block
            gT = min(gT, box_gap(a, inkB(j, :)));
        end
        gD = inf; for j = 1:size(dotsXY, 1), gD = min(gD, box_gap(a, [dotsXY(j,:) dotsXY(j,:)]) - rDot); end
        gB = inf;                                   % grid lines
        for k = 0:2
            gB = min([gB, box_gap(a, [k*cw - cfg.lwAxis/2, 0, k*cw + cfg.lwAxis/2, Hg]), ...
                          box_gap(a, [0, k*ch - cfg.lwAxis/2, Wg, k*ch + cfg.lwAxis/2])]);
        end
        g = min([gT gD gB]);
        gmin = min(gmin, g);
        fprintf('INK "%s": texts %.2f, dots %.2f, grid lines %.2f pt\n', names{i}, gT, gD, gB);
    end
    for k = 1:size(dotsXY, 1)
        gb = inf;
        for j = 0:2
            gb = min([gb, abs(dotsXY(k,1) - j*cw), abs(dotsXY(k,2) - j*ch)]);
        end
        fprintf('DOT %d: center to nearest grid line %.2f pt (outer edge %.2f pt clear)\n', k, gb, gb - rDot - cfg.lwAxis/2);
    end
    fprintf('INK MIN panel C: %.2f pt\n', gmin);
    fprintf('PANEL C: grid %.1f x %.1f pt (cells %.1f); labels right %.1f pt; depth %.1f pt\n', Wg, Hg, cw, Rr, depth);
end
