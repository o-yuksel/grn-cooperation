function R = fig6_figS2_data(cfg)
% FIG6_FIGS2_DATA  The data and checks behind Fig. 6, and the drawing of Fig. S2: the cut in the
% tracking task's mutation load when responsiveness or ||x*|| is lowered, the responsiveness of
% populations with the target moving over time, and the per-mutation load cells.
%
% Implements : Fig. 6's data (cfg.style 'data', for fig6_tracking_load) and Fig. S2 (cfg.style
%              'si_moving', for figS2_moving_target)
% Inputs     : cfg: style ('data' or 'si_moving'), pubDir (default the publication folder under the
%              data root), figDir (default output/figures); for 'data' also bigGroup (also read
%              tracking_mutation_regime/size_10x, mutations ten times larger) and extraHeldStill ({folder, a; ...},
%              further runs with the target held still)
%              publication/tracking_cut_after_150/design1 and design2/main.mat and
%                main.txt   (written by tracking_cut_after_150)
%              publication/tracking_moving_target/moving_target.mat and moving_target.txt   (the target
%                moving; written by run_tracking_moving_target)
%              for 'data': publication/tracking_runs/tracking_runs_N1000_load.txt
%                (written by tracking_load_vs_identity) and the .mat file fig6_load_cells('data')
%                reads
% Outputs    : R, each run's results (fields `held_still`, `moving`, and for 'data' `big`, `extra` and
%              the load cells `C`); style 'si_moving': [figDir '/fig6_figS2_data_moving'] .pdf,
%              .png (300 dpi) and .tif (1200 dpi, LZW); each run's changes and the checks, printed
% Seed       : none here; the bootstrap of fig6_load_cells uses its own fixed RandStream (twister)
%
% Each run's change is the paired difference d in mutation load against its control, with its
% interval, from the stored analysis table, also as a percent of the control's median load. The
% run stops unless every value equals its result text (4 decimals; the responsiveness medians of
% Fig. S2 to 0.006) and the eight load cells equal the per-cell means of the _load.txt file.

    if nargin < 1 || ~isfield(cfg, 'style') || ~any(strcmp(cfg.style, {'data', 'si_moving'}))
        error('fig6_figS2_data:style', 'cfg.style must be ''data'' or ''si_moving''');
    end
    thisDir = fileparts(mfilename('fullpath'));
    if ~isfield(cfg, 'pubDir'), cfg.pubDir = grn_data_file('publication'); end
    if ~isfield(cfg, 'figDir'), cfg.figDir = fullfile(fileparts(thisDir), 'output', 'figures'); end
    R.held_still = [held_still_design(cfg, 'tracking_cut_after_150/design1', 0.70), held_still_design(cfg, 'tracking_cut_after_150/design2', 0.20)];
    R.moving = moving_target(cfg);
    if strcmp(cfg.style, 'si_moving')
        draw_si_moving(R.moving, cfg);        % Fig. S2, the time course alone (figS2_moving_target)
        return
    end
    R.big = [];
    if isfield(cfg, 'bigGroup') && cfg.bigGroup   % target held still, a = 0.70, mutations x10 (tracking_mutation_regime/size_10x): Fig. 6B's 10x point
        R.big = held_still_design(cfg, 'tracking_mutation_regime/size_10x', 0.70);
    end
    R.C = fig6_load_cells('data', struct('dataDir', fullfile(cfg.pubDir, 'tracking_runs')));
    check_fig6a(R.C, cfg);
    % cfg.extraHeldStill = {folder, a; ...} reads further runs with the target held still through
    % held_still_design, each checked against its result text
    for j = 1:size(getfield_or(cfg, 'extraHeldStill', cell(0, 2)), 1)
        R.extra(j) = held_still_design(cfg, cfg.extraHeldStill{j, 1}, cfg.extraHeldStill{j, 2});
    end
end


% ===================== data + controls ===============================
function F = held_still_design(cfg, dirName, a)
    S = load(fullfile(cfg.pubDir, dirName, 'main.mat')); T = S.out.analysis.table;
    txt = fileread(fullfile(cfg.pubDir, dirName, 'main.txt'));
    line = regexp(txt, 'MUTATION LOAD:[^\n]*', 'match', 'once');
    F.a = a; F.name = dirName; F.control = median(S.out.analysis.Y(:, 1, 1), 'omitnan');
    tokC = regexp(line, 'CONTROL ([\d.]+) \|', 'tokens', 'once');
    arms = {'CONTROL2', 'ALPHA-DOWN', 'XNORM-DOWN'}; idx = [2 3 4];
    F.d = nan(1, 3); F.ci = nan(3, 2); F.n = nan(1, 3);
    ok = abs(F.control - str2double(tokC{1})) < 5e-5;
    for k = 1:3
        tok = regexp(line, [arms{k} ' [\d.]+, d ([+-][\d.]+) \[([+-][\d.]+), ([+-][\d.]+)\] n (\d+)'], 'tokens', 'once');
        F.d(k) = T(1, idx(k)).d; F.ci(k, :) = T(1, idx(k)).ci(:)'; F.n(k) = T(1, idx(k)).n;
        ok = ok && abs(F.d(k) - str2double(tok{1})) < 5.1e-5 && abs(F.ci(k, 1) - str2double(tok{2})) < 5.1e-5 ...
            && abs(F.ci(k, 2) - str2double(tok{3})) < 5.1e-5 && F.n(k) == str2double(tok{4});
    end
    F.pct = 100 * F.d / F.control; F.pctCI = 100 * F.ci / F.control;
    fprintf('HELD STILL %s (a = %.2f): control load %.4f; change vs control, percent: CONTROL2 %+.1f, ALPHA-DOWN %+.1f [%+.1f, %+.1f], XNORM-DOWN %+.1f [%+.1f, %+.1f]\n', ...
        dirName, a, F.control, F.pct(1), F.pct(2), F.pctCI(2, :), F.pct(3), F.pctCI(3, :));
    fprintf('  absolute d: ALPHA-DOWN %+.4f [%+.4f, %+.4f], XNORM-DOWN %+.4f [%+.4f, %+.4f]; CONTROL vs main.txt (4 dp): %s\n', ...
        F.d(2), F.ci(2, :), F.d(3), F.ci(3, :), figure_utils.pass_fail(ok));
    if ~ok, error('fig6_figS2_data:control', '%s does not reproduce its result text', dirName); end
end


function M = moving_target(cfg)
    S = load(fullfile(cfg.pubDir, 'tracking_moving_target', 'moving_target.mat')); a = S.out.analysis;
    txt = fileread(fullfile(cfg.pubDir, 'tracking_moving_target', 'moving_target.txt'));
    M.k = [0:10:90, 95:5:150];
    M.ctl = a.push.CONTROL.ra; M.low = a.push.ALPHA_DOWN.ra;
    if size(M.ctl, 2) ~= numel(M.k), error('fig6_figS2_data:control', 'time-course columns %d, expected %d', size(M.ctl, 2), numel(M.k)); end
    M.medC = median(M.ctl, 1, 'omitnan'); M.medL = median(M.low, 1, 'omitnan');
    M.q1L = prctile(M.low, 25, 1); M.q3L = prctile(M.low, 75, 1);
    M.q1C = prctile(M.ctl, 25, 1); M.q3C = prctile(M.ctl, 75, 1);
    sel = find(ismember(M.k, 0:30:150));
    ok = true;
    for nm = {'CONTROL', 'ALPHA-DOWN'}
        tok = regexp(txt, [nm{1} '\s+gen\+\[[^\]]*\]: alpha \[([^\]]*)\]'], 'tokens', 'once');
        v = sscanf(tok{1}, '%f')';
        if strcmp(nm{1}, 'CONTROL'), got = M.medC(sel); else, got = M.medL(sel); end
        ok = ok && numel(v) == numel(got) && all(abs(v - got) < 0.006);
    end
    load1 = regexp(txt, 'LOAD:[^\n]*', 'match', 'once');
    ben = regexp(txt, 'benefit \(best genotype window w\):[^\n]*', 'match', 'once');
    M.load = grab(load1, 'ALPHA-DOWN'); M.benefit = grab(ben, 'ALPHA-DOWN');
    tk = regexp(load1, 'CONTROL ([\d.]+) \|', 'tokens', 'once'); M.loadControl = str2double(tk{1});
    tk = regexp(ben, 'CONTROL ([\d.]+) \|', 'tokens', 'once'); M.benControl = str2double(tk{1});
    tk = regexp(ben, 'ALPHA-DOWN ([\d.]+),', 'tokens', 'once'); M.benLow = str2double(tk{1});
    Tl = a.table(1, 3); Tb = a.table(3, 3);   % LOAD row 1, benefit row 3 (moving_target.txt's order)
    ok = ok && abs(Tl.d - M.load(1)) < 5.1e-5 && abs(Tb.d - M.benefit(1)) < 5.1e-5;
    M.loadMat = Tl.d; M.benMat = Tb.d;
    % All three arms' LOAD differences (text, 4 dp) in the order CONTROL2, ALPHA-DOWN,
    % XNORM-DOWN, each checked against the stored table's d (columns 2, 3, 4).
    M.loadArms = [grab(load1, 'CONTROL2'); grab(load1, 'ALPHA-DOWN'); grab(load1, 'XNORM-DOWN')];
    for j = 1:3
        ok = ok && abs(a.table(1, j + 1).d - M.loadArms(j, 1)) < 5.1e-5;
        M.loadArmsTab(j, :) = [a.table(1, j + 1).d, a.table(1, j + 1).ci(:)'];
    end
    if isfield(a, 'Y')        % the control level is a median over replicates, as for the held-still runs
        M.loadControlMedian = median(a.Y(:, 1, 1), 'omitnan');
        ok = ok && abs(M.loadControlMedian - M.loadControl) < 5e-5;
    end
    fprintf('MOVING alpha-down: median alpha ratio at gen+150 %.3f (control %.3f); load d %+.4f [%+.4f, %+.4f]; benefit d %+.4f [%+.4f, %+.4f] (%.4f -> %.4f); controls vs moving_target.txt: %s\n', ...
        M.medL(end), M.medC(end), M.load, M.benefit, M.benControl, M.benLow, figure_utils.pass_fail(ok));
    if ~ok, error('fig6_figS2_data:control', 'moving-target time course or table does not reproduce moving_target.txt'); end
end


function v = grab(line, arm)
    tok = regexp(line, [arm ' [\d.]+, d ([+-][\d.]+) \[([+-][\d.]+), ([+-][\d.]+)\]'], 'tokens', 'once');
    v = str2double(tok);
end


function check_fig6a(TD, cfg)
    txt = fileread(fullfile(cfg.pubDir, 'tracking_runs', 'tracking_runs_N1000_load.txt'));
    seg = txt(1:strfind(txt, '---- L0:') - 1);
    tok = regexp(seg, '(moving|fixed)\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+\[', 'tokens');
    ok = numel(tok) == 8;
    for i = 1:numel(tok)
        arm = strcmp(tok{i}{1}, 'fixed') + 1; a = str2double(tok{i}{2}); m = str2double(tok{i}{3});
        k = find(TD.L.ai == arm & abs(TD.L.aList(TD.L.f) - a) < 1e-9);
        ok = ok && numel(k) == 1 && abs(TD.L.L(k) - m) < 5e-5 * max(1, 10 * m);
    end
    fprintf('Fig. 6A: eight cells against the per-cell means of tracking_runs_N1000_load.txt: %s\n', figure_utils.pass_fail(ok));
    if ~ok, error('fig6_figS2_data:control', 'Fig. 6A cells do not match tracking_runs_N1000_load.txt'); end
end


% ===================== Fig. S2 =======================================
function draw_si_moving(M, cfg)
% Fig. S2: the moving-target time course alone; control in ink (square end marker),
% responsiveness lowered in mid grey (circle end marker); no cooperation in this task, so no blue
% or red-orange.
    figure_utils.setup_defaults();
    W = 8.7; H = 6.2;
    c.untouched = figure_utils.INK; c.tested = figure_utils.GRAY;
    fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', figure_utils.FONT, 'DefaultTextFontName', figure_utils.FONT);
    guard = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Position', nr([1.35 1.05 7.0 4.9], W, H)); hold(ax, 'on');
    fill(ax, [M.k fliplr(M.k)], [M.q1C fliplr(M.q3C)], c.untouched, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax, [M.k fliplr(M.k)], [M.q1L fliplr(M.q3L)], c.tested, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % Each line ends in a filled marker (control square, responsiveness lowered circle), as
    % Fig. S1B draws its end markers; the legend glyphs carry line and marker.
    lines_ = {M.medC, M.medL}; cols_ = {c.untouched, c.tested}; mk_ = {'s', 'o'};
    for a = 1:2
        plot(ax, M.k, lines_{a}, '-', 'Color', cols_{a}, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax, M.k(end), lines_{a}(end), mk_{a}, 'MarkerSize', 5, 'LineWidth', 0.8, 'MarkerFaceColor', cols_{a}, ...
            'MarkerEdgeColor', cols_{a}, 'Clipping', 'off', 'HandleVisibility', 'off');
        hk(a) = plot(ax, NaN, NaN, '-', 'Color', cols_{a}, 'LineWidth', 1.2, 'Marker', mk_{a}, ...
            'MarkerSize', 5, 'MarkerFaceColor', cols_{a}, 'MarkerEdgeColor', cols_{a}); %#ok<AGROW>
    end
    hC = hk(1); hL = hk(2);
    set(ax, 'XLim', [0 150], 'YLim', [0 1.25], 'XTick', 0:50:150, 'YTick', 0:0.25:1.25, ...
        'YTickLabel', {'0', '0.25', '0.5', '0.75', '1', '1.25'});
    xlabel(ax, 'generations after split'); ylabel(ax, 'responsiveness / control');
    lg = legend(ax, [hC hL], {'control', 'responsiveness lowered'}, 'Box', 'off', 'FontSize', figure_utils.SIZE_SMALL, ...
        'Location', 'southeast', 'TextColor', figure_utils.INK);
    lg.ItemTokenSize = [14 7];
    hold(ax, 'off'); figure_utils.style_panel(ax);
    export_page(fig, fullfile(cfg.figDir, 'fig6_figS2_data_moving'), W, H);
end


function export_page(fig, base, W, H)
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


function p = nr(r, W, H), p = [r(1) / W, r(2) / H, r(3) / W, r(4) / H]; end


function v = getfield_or(s, f, d), if isfield(s, f), v = s.(f); else, v = d; end, end
