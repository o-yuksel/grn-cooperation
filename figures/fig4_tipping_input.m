function fig4_tipping_input(pubDir)
% FIG4_TIPPING_INPUT  Draws Fig. 4: (A) the share of mutations that end a member's cooperation or
% retaliation, by the input each adds over its tipping input; (B) the same shares after selection
% for safety, for members that must also retaliate and for those that need only cooperate; (C)
% the share that end cooperation with the clone as the response to the partner is scaled.
%
% Produces : Fig. 4; SI S3.3 (0.6992 and 0.7452 of losses caught by the tipping input, as a check)
% Inputs   : pubDir (optional; default the publication folder under the data root)
%            publication/tipping_input/tipping_input_7000.mat   (written by
%              tipping_input_flag('run', 7000))
%            publication/safety_search/long_flags.mat   (written by safety_search_losses('long'))
%            publication/safety_search/fig4b_retaliation_rest.txt   (written by
%              fig4b_retaliation_from_rest)
%            publication/conditionality_score/response_dial/response_dial_rest.mat and .txt
%              (written by response_dial_from_rest)
% Outputs  : output/figures/Fig4.pdf (vector), .png (300 dpi) and .tif (1200 dpi, LZW), under the
%            release root; prints the checks, panel A's table by range, panel B's summary per
%            group, panel C's means and intervals, and panel C's numbers for the caption
% Run      : fig4_tipping_input
% Seed     : rng(20261017 + 30, 'twister') for panel A's 2,000 bootstraps, rng(20261024,
%            'twister') for panel B's horizontal jitter; each restores the caller's stream
%
% A mutation changes the weight from source gene j to gene i by dW and adds the input e = dW * x_j,
% x_j the source gene's expression against the clone (cooperation) or against the worst-case
% mutant (retaliation); q = |e| / tipping input, the smallest input of that sign at gene i that
% causes the loss. Each mutation is weighted by its member's number of individuals.

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1, pubDir = grn_data_file('publication'); end
    D = load_added_input(fullfile(pubDir, 'tipping_input', 'tipping_input_7000.mat'));
    check_control(D);
    lab2 = label_check(D);
    T = cell(1, 2);
    for s = 1:2, T{s} = bin_counts(D.q{s}, D.y{s}); end
    R = boot_bins(T, D.w, D.arch, 2000, 20261017 + 30);
    print_table(R, lab2);
    B = load_safety(fullfile(pubDir, 'safety_search', 'long_flags.mat'));
    print_safety(B);
    C = load_rest(fullfile(pubDir, 'conditionality_score', 'response_dial'));
    fprintf('\n[PANEL C] from-rest rule, n = %d: L %s, lower %s, upper %s; retaliating percent %s\n', C.n, mat2str(C.L, 4), ...
        mat2str(C.ci(:, 1)', 4), mat2str(C.ci(:, 2)', 4), mat2str(C.pct));
    fprintf('%s\n', C.cap);
    outDir = fullfile(root, 'output', 'figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    draw_figure(R, lab2, B, C, fullfile(outDir, 'Fig4'), pubDir);
end


function D = load_added_input(f)
    L = load(f); S = L.S; R = L.R;
    mat = @(fn) cell2mat(arrayfun(@(s) s.(fn)(:)', S(:), 'UniformOutput', false));
    Rm = @(fn) cell2mat(cellfun(@(o) o.(fn)(:)', R(:), 'UniformOutput', false));
    slot = mat('slot'); dW = mat('dW'); lost = mat('lost') ~= 0; rH = mat('rH');
    gene = mod(slot - 1, 4) + 1;                      % row index i of ind2sub([4 4], slot)
    thr = {Rm('thrC'), Rm('thrR')}; xj = {Rm('xjC'), Rm('xjR')};
    n = numel(S); K = size(slot, 2); D.q = cell(1, 2); D.y = {lost, rH > 0}; D.xj = xj;
    for s = 1:2
        e = dW .* xj{s}; col = gene + 4 * (e < 0);    % thr is 4x2 column-major: (i, sign)
        th = thr{s}(sub2ind([n 8], repmat((1:n)', 1, K), col));
        q = abs(e) ./ th; q(e == 0 | ~isfinite(e) | isinf(th) | ~isfinite(q)) = 0;
        D.q{s} = q;
    end
    D.w = [S.count]'; D.arch = [S.arch]'; D.n = n; D.K = K;
    fprintf('loaded %d members, %d mutations each\n', n, K);
end


function check_control(D)
    want = [0.6992 0.7452; 0.4400 0.6426]; got = zeros(2);
    W = repmat(D.w, 1, D.K);
    for s = 1:2
        P = D.q{s} >= 1; Y = D.y{s};
        got(1, s) = sum(W(P & Y)) / sum(W(Y)); got(2, s) = sum(W(P & Y)) / sum(W(P));
    end
    ok = all(abs(round(got, 4) - want) < 5e-5, 'all');
    fprintf('[CONTROL] q >= 1 (cooperation / retaliation): recall %.4f / %.4f (want 0.6992 / 0.7452), precision %.4f / %.4f (want 0.4400 / 0.6426): %s\n', ...
        got(1, 1), got(1, 2), got(2, 1), got(2, 2), pass_fail(ok));
    if ~ok, error('CONTROL FAILED: q definition does not reproduce the tipping-input flag'); end
end


function s = pass_fail(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end


function lab = label_check(D)
% share (clone-weighted) of the 0 < q < 1/100 group whose source-gene expression before the mutation is < 0.01
    W = repmat(D.w, 1, D.K); sh = zeros(1, 2); grp = zeros(1, 2);
    for s = 1:2
        m = D.q{s} > 0 & D.q{s} < 1e-2;
        sh(s) = sum(W(m & D.xj{s} < 0.01)) / sum(W(m)); grp(s) = sum(W(m)) / sum(W(:));
    end
    fprintf('[LABEL CHECK] group 0 < q < 1/100: share of mutations %.4f / %.4f; of these, source-gene expression < 0.01: %.4f / %.4f (cooperation / retaliation)\n', ...
        grp(1), grp(2), sh(1), sh(2));
    if all(sh >= 0.95), lab = {'source', 'gene off', 'before'}; else, lab = {'added input <', '1/100'}; end
    fprintf('  -> tick label: "%s"\n', strjoin(lab, ' '));
end


function b = bin_index(q)
% 0 = no added input (q = 0); 1 = 0 < q < 1/100; 2..13 = log10 q in [-2 + 0.25 (b-2), -2 + 0.25 (b-1)); 14 = q >= 10.
    b = zeros(size(q));
    b(q > 0 & q < 1e-2) = 1;
    mid = q >= 1e-2 & q < 10;
    b(mid) = min(13, floor((log10(q(mid)) + 2) / 0.25) + 2);
    b(q >= 10) = 14;
end


function T = bin_counts(q, y)
% per member: mutations in each column (1 + bin_index) and those among them with the loss
    b = bin_index(q) + 1; n = size(q, 1); T.cnt = zeros(n, 15); T.los = zeros(n, 15);
    for k = 1:15
        m = b == k; T.cnt(:, k) = sum(m, 2); T.los(:, k) = sum(m & y, 2);
    end
end


function R = boot_bins(T, w, arch, B, seed)
    rs = rng; rng(seed, 'twister'); n = numel(w); nb = 15; K = sum(T{1}.cnt(1, :));
    R.frac = zeros(nb, 2); R.share = zeros(nb, 2); bs = nan(B, nb, 2);
    for s = 1:2
        R.frac(:, s) = (w' * T{s}.los)' ./ (w' * T{s}.cnt)'; R.share(:, s) = (w' * T{s}.cnt)' / (sum(w) * K);
    end
    for b = 1:B
        c = zeros(n, 1);
        for m = unique(arch)'
            ix = find(arch == m); c = c + accumarray(ix(randi(numel(ix), numel(ix), 1)), 1, [n 1]);
        end
        cw = c .* w;
        for s = 1:2, bs(b, :, s) = (cw' * T{s}.los) ./ (cw' * T{s}.cnt); end
    end
    rng(rs);
    R.ci = zeros(nb, 2, 2);
    for s = 1:2, R.ci(:, s, :) = reshape(prctile(bs(:, :, s), [2.5 97.5], 1)', nb, 1, 2); end
end


function print_table(R, lab2)
    fprintf('\n%-16s | cooperation: frac [CI]  share  | retaliation: frac [CI]  share\n', 'log10 q bin');
    for k = 1:15
        if k == 1, lab = 'no added input';
        elseif k == 2, lab = ['q < 1/100 (' strjoin(lab2, ' ') ')'];
        elseif k == 15, lab = 'q >= 10';
        else, lo = -2 + 0.25 * (k - 3); lab = sprintf('[%+.2f, %+.2f)', lo, lo + 0.25); end
        fprintf('%-16s |', lab);
        for s = 1:2
            fprintf(' %.4f [%.4f, %.4f]  %.4f |', R.frac(k, s), R.ci(k, s, 1), R.ci(k, s, 2), R.share(k, s));
        end
        fprintf('\n');
    end
    fprintf('share of q >= 1 (columns 11-15): %.4f / %.4f\n', sum(R.share(11:15, 1)), sum(R.share(11:15, 2)));
end


function draw_figure(R, lab2, B, C, outBase, pubDir)
    BLUE = figure_utils.COLOR_COOP; VER = figure_utils.COLOR_DEFECT; INK = figure_utils.INK; NOTE = figure_utils.NOTE;
    % Sizes in inches: panel A's band 2.85 in above the lower row of B and C (hB) and the key (pad)
    dlt = 0.12; pad = 0.10;   % dlt: bottom slack removed; pad: top band for the key
    W = 3.42; hB = 2.30 + 0.48 - dlt; H = 2.85 + hB + pad;   % +0.48 in: the strips under panels B and C
     xl = [-4.6 2.2]; xNo = -4.3; xOff = -3.2; xTop = 1.7;
    fig = figure('Units', 'inches', 'Position', [1 1 W H], 'Color', 'w');
    ax = axes(fig, 'Position', [0.19 (1.25 + hB) / H 0.62 1.45 / H]);
    st = axes(fig, 'Position', [0.19 (0.70 + hB) / H 0.62 0.45 / H]);
    x = [xNo, xOff, -2 + 0.25 * ((1:12) - 0.5), xTop];  % one x per column of R
    cols = {BLUE, VER};
    mk = {'o', 's'}; msz = [2.8 2.5];            % circle / square: the square is drawn ~10% smaller to read equal
    hold(ax, 'on');
    plot(ax, [0 0], [0 1], '-', 'Color', NOTE, 'LineWidth', 0.7);
    for s = 1:2
        v = R.frac(:, s); lo = R.ci(:, s, 1); hi = R.ci(:, s, 2);
        for seg = {3:10, 11:14}                          % 1/100 <= q < 1 and 1 <= q < 10; nothing else is joined
            plot(ax, x(seg{1}), v(seg{1}), '-', 'Color', cols{s}, 'LineWidth', 0.8);
        end
        for k = 1:15, plot(ax, [x(k) x(k)], [lo(k) hi(k)], '-', 'Color', cols{s}, 'LineWidth', 0.8); end
        plot(ax, x, v, mk{s}, 'MarkerFaceColor', cols{s}, 'MarkerEdgeColor', cols{s}, 'MarkerSize', msz(s), 'LineWidth', 0.8, 'LineStyle', 'none');
    end
    draw_key(ax, cols, mk, msz, INK);
    style_axes(ax, xl, [0 1], 0:0.5:1, INK);
    set(ax, 'XTick', [], 'XColor', 'none');
    ylabel(ax, 'fraction of mutations', 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK);
    draw_strip(st, R, x, xl, [xNo xOff], lab2, cols, INK);
    xlabel(st, ['added input ' char(247) ' tipping input for that loss'], 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK, ...
        'Units', 'normalized', 'Position', [0.5 -0.95 0], 'VerticalAlignment', 'top');
    ymax = draw_safety(fig, B, [0.19 (1.18 - dlt) / H 0.36 1.50 / H], INK);
    draw_safety_strip(fig, [0.19 (0.70 - dlt) / H 0.36 0.38 / H], INK, pubDir);   % same height, gap and style as C's strip
    draw_dial(fig, C, [0.66 (1.18 - dlt) / H 0.27 1.50 / H], [0.66 (0.70 - dlt) / H 0.27 0.38 / H], ymax, INK);
    for a = {ax, 'A', -0.30; findobj(fig, 'Tag', 'panelB'), 'B', -0.30 * 0.62 / 0.36; findobj(fig, 'Tag', 'panelC'), 'C', -0.40}'
        text(a{1}, a{3}, 1.04, a{2}, 'Units', 'normalized', 'FontSize', figure_utils.SIZE_LETTER, 'FontWeight', 'bold', ...
            'FontName', figure_utils.FONT, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Clipping', 'off');
    end
    [~, iss] = validate_figure(fig, struct('name', 'Fig4'));
    for i = 1:numel(iss), fprintf('validate_figure: %s\n', iss{i}); end
    ws = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector', 'Padding', 8);
    warning(ws);
    exportgraphics(fig, [outBase '.png'], 'Resolution', 300, 'Padding', 8);
    exportgraphics(fig, [outBase '.tif'], 'Resolution', 1200, 'Padding', 8);
    imwrite(imread([outBase '.tif']), [outBase '.tif'], 'Compression', 'lzw', 'Resolution', 1200);
    fprintf('fig4_tipping_input: wrote %s.{pdf,png,tif}\n', outBase);
end


function style_axes(ax, xl, yl, yt, INK)
    set(ax, 'XLim', xl, 'YLim', yl, 'YTick', yt, 'FontSize', figure_utils.SIZE_SMALL, 'FontName', figure_utils.FONT, 'LineWidth', 0.6, ...
        'TickDir', 'out', 'Box', 'off', 'XColor', INK, 'YColor', INK);
end


function draw_strip(st, R, x, xl, xCat, lab2, cols, INK)
    bw = 0.12; hold(st, 'on');
    plot(st, [0 0], [0 0.5], '-', 'Color', figure_utils.NOTE, 'LineWidth', 0.7);
    for s = 1:2
        for k = 1:15
            fill(st, x(k) + (s - 1.5) * bw + bw / 2 * [-1 1 1 -1] * 0.9, [0 0 R.share(k, s) R.share(k, s)], cols{s}, 'EdgeColor', 'none');
        end
    end
    style_axes(st, xl, [0 0.5], [0 0.25 0.5], INK);
    set(st, 'XTick', [xCat -2 -1 0 1 x(15)], 'XTickLabel', {'', '', '1/100', '1/10', '1', '10', [char(8805) '10']});
    set(st, 'XTickLabelRotation', 0);
    labs = {{'no', 'input'}, lab2};
    for k = 1:2
        text(st, xCat(k), -0.06, labs{k}, 'FontSize', figure_utils.SIZE_SMALL, 'FontName', figure_utils.FONT, 'Color', INK, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Clipping', 'off');
    end
    ylabel(st, {'share of', 'mutations'}, 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK);
end


function B = load_safety(f)
% Per network, the share of its mutations with each loss (stored flags lostC, lostR): C0, R0 at the
% start (grp 1); CS, RS after selection for safety in the group that must stay conditional (grp 2:
% cooperate with the clone and retaliate); CC in the group that need only cooperate (grp 3); E0, ES
% either loss. Panel B draws CS, RS and CC.
    L = load(f); N = L.N; if ~iscell(N), N = num2cell(N); end
    grp = [L.J.grp]; col = @(v) v(:);
    share = @(g, fn) col(cellfun(@(n) mean(n.(fn)), N(grp == g)));
    either = @(g) col(cellfun(@(n) mean(n.lostR | n.lostC), N(grp == g)));
    B.C0 = share(1, 'lostC'); B.R0 = share(1, 'lostR'); B.CS = share(2, 'lostC'); B.RS = share(2, 'lostR');
    B.CC = share(3, 'lostC'); B.E0 = either(1); B.ES = either(2);
end


function print_safety(B)
    fprintf('\n[PANEL B] per-network fraction of mutations causing the loss\n');
    fprintf('%-4s %4s %9s %9s %9s %6s\n', 'var', 'n', 'median', 'min', 'max', 'zeros');
    for f = {'C0', 'R0', 'CS', 'RS', 'CC', 'E0', 'ES'}
        v = B.(f{1}); fprintf('%-4s %4d %9.4f %9.4f %9.4f %6d\n', f{1}, numel(v), median(v), min(v), max(v), sum(v == 0));
    end
    fprintf('(E0/ES = either loss, mean of lostC|lostR; C0,R0,E0 start, CS,RS,ES after the search that must also retaliate; CC after the search that need only cooperate)\n');
end


function ymax = draw_safety(fig, B, pos, INK)
    BLUE = figure_utils.COLOR_COOP; VER = figure_utils.COLOR_DEFECT;
    ax = axes(fig, 'Position', pos, 'Tag', 'panelB'); hold(ax, 'on');
    rs = rng; rng(20261024, 'twister');
    clouds = {B.CS, B.RS, B.CC}; xs = [0.85 1.15 2]; cc = {BLUE, VER, BLUE}; mk = {'o', 's', 'o'}; sz = [14 11 14];
    for i = 1:3
        y = clouds{i}; x = xs(i) + (rand(size(y)) - 0.5) * 0.18;
        scatter(ax, x, y, sz(i), cc{i}, mk{i}, 'filled', 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
        m = median(y); plot(ax, xs(i) + [-0.13 0.13], [m m], '-', 'Color', 'k', 'LineWidth', 1.4);
    end
    rng(rs);
    ymax = 0.25; top = max([B.CS; B.RS; B.CC]);
    if top > ymax, ymax = ceil(top / 0.05) * 0.05; fprintf('[PANEL B] y limit raised to %.2f (max value %.4f)\n', ymax, top); end
    style_axes(ax, [0.5 2.5], [0 ymax], 0:0.05:ymax, INK);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'', ''});      % group labels and x label sit under the strip below
    ylabel(ax, 'fraction of mutations', 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK);
end


function bw = bar_width(fig, pos, xrange)
% Bar width in data units for a printed width of BAR_IN inches, shared by the strips of B and C.
    BAR_IN = 0.123;
    bw = BAR_IN * xrange / (pos(3) * fig.Position(3));
end


function draw_safety_strip(fig, pos, INK, pubDir)
% Strip under panel B: the share of each group's networks that retaliate from rest at the evolved
% response (panel C's rule at s = 1), read from the report of fig4b_retaliation_from_rest. The
% run stops unless the report's check passed, each group has 30 networks and every retaliating
% network is also conditional (cooperates with the clone and retaliates).
    NOTE = figure_utils.NOTE;
    f = fullfile(pubDir, 'safety_search', 'fig4b_retaliation_rest.txt');
    t = fileread(f);
    if contains(t, 'Not computed') || ~contains(t, '-> PASS'), error('fig4b_retaliation_rest report: CHECK not PASS in %s', f); end
    nm = {'must stay conditional', 'need only cooperate'}; fr = nan(1, 2);
    for k = 1:2
        tk = regexp(t, ['  ' nm{k} ': n (\d+); cooperate (\d+); retaliate (\d+); conditional (\d+);'], 'tokens', 'once');
        if isempty(tk), error('count line for "%s" not found in %s', nm{k}, f); end
        v = str2double(tk);
        if v(1) ~= 30, error('panel B strip: %s has n = %d, not 30', nm{k}, v(1)); end
        if v(3) ~= v(4), error('panel B strip: %s retaliate %d differs from conditional %d', nm{k}, v(3), v(4)); end
        fr(k) = v(4) / v(1);
    end
    pct = round(100 * fr);
    fprintf('[PANEL B] conditional from rest: %s percent\n', mat2str(pct));
    st = axes(fig, 'Position', pos, 'Tag', 'panelBstrip'); hold(st, 'on'); bw = bar_width(fig, pos, 2);   % same printed width as C's bars (bar_width)
    for k = 1:2
        fill(st, k + bw / 2 * [-1 1 1 -1], [0 0 fr(k) fr(k)], NOTE, 'EdgeColor', 'none');
        text(st, k, fr(k) + 0.06, sprintf('%d', pct(k)), 'Color', INK, 'FontSize', figure_utils.SIZE_SMALL, 'FontName', figure_utils.FONT, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');
    end
    style_axes(st, [0.5 2.5], [0 1.3], [0 1], INK);
    set(st, 'XTick', [1 2], 'XTickLabel', {'', ''});      % two-line labels drawn as text (narrow axes, panel C beside)
    labs = {{'must also', 'retaliate'}, {'need only', 'cooperate'}};
    for k = 1:2
        text(st, k, -0.1, labs{k}, 'FontSize', figure_utils.SIZE_SMALL, 'FontName', figure_utils.FONT, 'Color', INK, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Clipping', 'off');
    end
    set(st, 'YTick', []);                                 % bars and percentages carry the values
    ylabel(st, {'networks', 'that', 'retaliate', '(%)'}, 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK);   % serves C's strip too (same scale); 4 lines keep it inside the strip's height
    xlabel(st, 'after selection for safety', 'FontName', figure_utils.FONT, 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK, ...
        'Units', 'normalized', 'Position', [0.5 -1.25 0], 'VerticalAlignment', 'bottom');
end


function C = load_rest(rd)
% Panel C data: L, per network the share of its mutations that end cooperation with the clone at
% each response setting s = 0, 0.5, 1, 1.5, 2 (from rest); the mean over networks with a 95%
% cluster-t interval over replicates. The run stops unless there are 267 networks and the means
% equal CHECK 2 of response_dial_rest.txt to 4 significant digits.
    Z = load(fullfile(rd, 'response_dial_rest.mat')); S = Z.S; L = Z.L; rep = Z.rep(:);
    t = fileread(fullfile(rd, 'response_dial_rest.txt'));
    tv = regexp(t, 'CHECK 2[^\n]*?(\[[\d. ]+\]) \(want', 'tokens', 'once');
    if isempty(tv), error('CHECK 2 vector not found in response_dial_rest.txt'); end
    ok = size(L, 1) == 267 && isequal(mat2str(mean(L, 1), 4), tv{1});
    fprintf('[CONTROL] panel C: n %d (want 267), mean L %s (want %s): %s\n', size(L, 1), mat2str(mean(L, 1), 4), tv{1}, pass_fail(ok));
    if ~ok, error('CONTROL FAILED: panel C n or means differ from response_dial_rest.txt CHECK 2'); end
    C.n = size(L, 1); C.s = [0 0.5 1 1.5 2]; C.L = mean(L, 1); C.ci = zeros(5, 2);
    for k = 1:5, [~, C.ci(k, :)] = cluster_t(L(:, k), rep); end
    C.fr = mean(S.ret, 1); C.pct = round(100 * C.fr);
    C.cap = caption_k(Z, t);
end


function c = caption_k(Z, t)
% The networks conditional at both s = 1 and s = 2 (K), printed for the caption, not drawn; the run
% stops unless their count and mean L at both settings equal the report's n_K and "on K" line.
    S = Z.S; L = Z.L; tn = regexp(t, 'n_K = (\d+)', 'tokens', 'once');
    tm = regexp(t, 'on K: mean L\(1\) ([\d.]+), mean L\(2\) ([\d.]+)', 'tokens', 'once');
    if isempty(tn) || isempty(tm), error('n_K or "on K" line not found in response_dial_rest.txt'); end
    K = S.cond(:, 3) & S.cond(:, 5); m = [mean(L(K, 3)) mean(L(K, 5))]; want = [str2double(tm{1}) str2double(tm{2})];
    ok = sum(K) == str2double(tn{1}) && isequal(sprintf('%.5f %.5f', m), sprintf('%.5f %.5f', want));
    fprintf('[CONTROL] panel C caption K: n_K %d (want %s), means %.5f / %.5f (want %s / %s): %s\n', sum(K), tn{1}, m(1), m(2), tm{1}, tm{2}, pass_fail(ok));
    if ~ok, error('CONTROL FAILED: K subset differs from response_dial_rest.txt'); end
    ci = zeros(2, 2); [~, ci(1, :)] = cluster_t(L(K, 3), Z.rep(K)); [~, ci(2, :)] = cluster_t(L(K, 5), Z.rep(K));
    c = sprintf('[PANEL C, caption] K n %d: %.5f -> %.5f; intervals [%.5f, %.5f] and [%.5f, %.5f]', sum(K), m(1), m(2), ci(1, :), ci(2, :));
end


function draw_dial(fig, C, pos, posStrip, ymax, INK)
% Panel C: the share of mutations that end cooperation with the clone, by response setting (from
% rest), mean over the networks, with its interval.
    BLUE = figure_utils.COLOR_COOP; LIGHT = 0.45 * BLUE + 0.55 * [1 1 1]; xs = C.s;
    ax = axes(fig, 'Position', pos, 'Tag', 'panelC'); hold(ax, 'on');
    plot(ax, xs, C.L, '-', 'Color', LIGHT, 'LineWidth', 0.8);
    for j = 1:numel(xs), plot(ax, [xs(j) xs(j)], C.ci(j, :), '-', 'Color', BLUE, 'LineWidth', 0.8); end
    plot(ax, xs, C.L, 'o', 'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'MarkerSize', 3.6, 'LineWidth', 0.8, 'LineStyle', 'none');
    style_axes(ax, [-0.25 2.25], [0 ymax], 0:0.05:ymax, INK);
    set(ax, 'XTick', xs, 'XTickLabel', {});          % x axis line and ticks as B's main axes; labels sit on the strip below
    draw_retaliating(fig, posStrip, xs, C.fr, C.pct, INK);
end


function draw_retaliating(fig, pos, xs, fr, pct, INK)
% Strip under panel C: share of the 267 networks that retaliate at each setting (computed from S.ret in load_rest).
    NOTE = figure_utils.NOTE;
    st = axes(fig, 'Position', pos, 'Tag', 'panelCstrip'); hold(st, 'on'); bw = bar_width(fig, pos, 2.5);
    for j = 1:numel(xs)
        fill(st, xs(j) + bw / 2 * [-1 1 1 -1], [0 0 fr(j) fr(j)], NOTE, 'EdgeColor', 'none');
        text(st, xs(j), fr(j) + 0.06, sprintf('%d', pct(j)), 'Color', INK, 'FontSize', figure_utils.SIZE_SMALL, 'FontName', figure_utils.FONT, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');
    end
    style_axes(st, [-0.25 2.25], [0 1.3], [0 1], INK);
    set(st, 'XTick', xs, 'XTickLabel', {'0', char(189), '1', '1.5', '2'}, 'XTickLabelRotation', 0, 'FontSize', figure_utils.SIZE_SMALL);
    set(st, 'YTick', []);                            % bars and percentages carry the values
    xlabel(st, ['response (' char(215) ' evolved)'], 'FontName', figure_utils.FONT, 'FontSize', figure_utils.SIZE_LABEL, 'Color', INK, ...
        'Units', 'normalized', 'Position', [0.42 -1.25 0], 'VerticalAlignment', 'bottom');
end


function draw_key(ax, cols, mk, msz, INK)
% One key for the whole figure, above panel A: colored marker, dark text.
    nm = {'cooperation lost', 'retaliation lost'}; xm = [-4.4 -1.3]; y = 1.10;
    for k = 1:2
        plot(ax, xm(k), y, mk{k}, 'MarkerFaceColor', cols{k}, 'MarkerEdgeColor', cols{k}, 'MarkerSize', msz(k), ...
            'LineWidth', 0.8, 'LineStyle', 'none', 'Clipping', 'off');
        text(ax, xm(k) + 0.2, y, nm{k}, 'Color', INK, 'FontSize', figure_utils.SIZE_SMALL, 'FontName', figure_utils.FONT, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Clipping', 'off');
    end
end
