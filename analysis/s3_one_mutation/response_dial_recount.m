function response_dial_recount()
% RESPONSE_DIAL_RECOUNT  Do losses rise along the response dial when exactly the same mutations are
% counted at every setting? From the stored outcomes of response_dial (no dynamics are run),
% recounts the share of each network's kept mutations that end cooperation (toggles of the four
% scaled connections left out, since their size scales with s), with four contrasts between
% settings, the checks on them, and alpha and the response slope along the dial.
%
% Produces : SI S3.5; Table S8 (prints the mutations left out, 400 of 8,010 from rest; the medians
%            of alpha and the response slope along the dial; that kept mutations change the wiring
%            alike at every setting)
% Inputs   : publication/conditionality_score/response_dial/response_dial_inputs.mat and chunks/
%              (written by response_dial)
% Outputs  : publication/conditionality_score/response_dial/response_dial_recount.txt
% Run      : response_dial_recount   (under 1 min, stored report)
% Seed     : K.seed = 20261028, Mersenne twister: every sign-flip test K.seed + 54000 (the same
%            10,000 sign vectors); check 3 K.seed + 70000 + i (permutation) and 72000 + i (signs),
%            i = 1 to 200
%
% L(n, s): the share of network n's kept mutants that lose cooperation at setting s (from rest:
% rM; developmental rules at s <= 1: cMut or eMut, as in dial_loss_shares). Contrasts: STEP
% L(1) - L(0), SLOPE L(2) - L(0.5), LOWER L(1) - L(0.5), UPPER L(2) - L(1).

    t0 = tic;
    [K, Hd] = response_dial('handles');
    L = load(fullfile(K.rd, 'response_dial_inputs.mat'), 'I');
    D = Hd.rest_data(K, L.I, fullfile(K.rd, 'chunks'));     % the from-rest set (338) and its subset in reaction fifths 2 to 5 (267; field prim)
    f = fullfile(K.rd, 'response_dial_recount.txt');

    %% Report
    report = @() recount_report(K, Hd, D, f, t0); %#ok<NASGU> called inside evalc
    rep = evalc('report()');
    fprintf('%s', rep);

    %% Save
    fid = fopen(f, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function recount_report(K, Hd, D, f, t0)
% The checks, the recount and the descriptive tables; prints the report's lines.
    fprintf('===== RESPONSE DIAL RECOUNT %s UTC =====\n', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('chunks %s; set %d, primary %d\n', ...
        D.dir, sum(D.set), sum(D.prim));
    E = Hd.exclusion(D);
    C = struct('c1', check1(K, D), 'c2', check2(Hd, D, E));
    P = recount_base(D, E);
    C.c3 = check3(K, P);
    R = contrasts4(K, P.L, P.rp);
    recount_lines(R, C);
    fprintf('\n==================== also reported ====================\n');
    removed(K, D, E, P, R);
    unsettled_too(K, D, E);
    dev_rules(K, Hd, D, E);
    alpha_table(K, D);
    decomp(Hd, D, E);
    fprintf('\ncompute time %.1f min; report %s\n', toc(t0) / 60, f);
end


function L = Lx(lost, keep)
% L_x: share of kept mutants lost, per network and setting (NaN where none is kept).
    L = squeeze(sum(lost & keep, 2)) ./ sum(keep, 2);
end


function P = recount_base(D, E)
    p = D.prim; keep = E.keep(p, :); nk = sum(keep, 2); use = nk > 0; rp = D.rep(p);
    L = Lx(D.lost(p, :, :), keep);
    P = struct('L', L(use, :), 'rp', rp(use), 'nk', nk(use), 'ndrop', sum(~use));
end


% =========================================================================
% CONTRASTS, TEST AND LABELS
% =========================================================================
function R = contrasts4(K, L, rp)
% Four contrasts between settings (columns 1 to 5 are s = 0, 0.5, 1, 1.5, 2), each with a
% cluster-t interval, a sign-flip p and the Holm-adjusted p.
    C = {'STEP',  'L(1) - L(0)',   3, 1; 'SLOPE', 'L(2) - L(0.5)', 5, 2; ...
         'LOWER', 'L(1) - L(0.5)', 3, 2; 'UPPER', 'L(2) - L(1)',   5, 3};
    R = struct('name', C(:, 1), 'def', C(:, 2), 'e', 0, 'ci', [0 0], 'p', 1);
    for j = 1:4
        d = L(:, C{j, 3}) - L(:, C{j, 4}); [R(j).e, R(j).ci] = cluster_t(d, rp);
        R(j).p = signflip(K, d, rp);
    end
    h = holm([R.p]); for j = 1:4, R(j).holm = h(j); end
end


function h = holm(p)
% Holm: adjusted p_(i) = max over j <= i of min(1, (m - j + 1) p_(j)), m = numel(p).
    m = numel(p); [ps, o] = sort(p); a = min(1, (m - (1:m) + 1) .* ps); a = cummax(a); h = zeros(size(p)); h(o) = a;
end


function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end


% =========================================================================
% CHECKS 1-3
% =========================================================================
function ok = check1(K, D)
% Check 1, a reproduction check: with no mutant left out, L by setting and the STEP and SLOPE
% contrasts print as in the stored response_dial_analysis.txt (the fixed strings below); reported.
    p = D.prim; L = Lx(D.lost(p, :, :), true(sum(p), 30)); rp = D.rep(p);
    sL = mat2str(mean(L, 1), 4); want = '[0.01211 0.01898 0.04257 0.07728 0.08602]';
    s = cell(2, 1); W = {'0.0305 [0.0244, 0.0365] p 0.0001', '0.0670 [0.0580, 0.0761] p 0.0001'}; ks = [3 1; 5 2];
    for j = 1:2
        d = L(:, ks(j, 1)) - L(:, ks(j, 2)); [e, ci] = cluster_t(d, rp);
        s{j} = sprintf('%.4f [%.4f, %.4f] p %.4f', e, ci, signflip(K, d, rp));
    end
    ok = strcmp(sL, want) && strcmp(s{1}, W{1}) && strcmp(s{2}, W{2});
    fprintf('\nCHECK 1 (reproduction, no mutant excluded): L by setting %s (want %s); STEP %s (want %s); SLOPE %s (want %s) -> %s\n', ...
        sL, want, s{1}, W{1}, s{2}, W{2}, pf(ok));
end


function ok = check2(Hd, D, E)
% Check 2: the change each mutant makes to the wiring, dW_s = Mm .* Bm - M .* B_s, is the same at
% every setting for kept mutants and s times its value at s = 1 for left-out ones (to 1e-12).
% Printed beside, not part of the check: the recomputed constant input against the stored one.
    Z = D.Z; id = find(D.prim); S = Z.S; dK = 0; dE = 0; dB = 0; nKp = 0; nEx = 0; dial = Hd.dial; apply_mut = Hd.apply_mut;
    for i = id'
        M = Z.M(:, :, i); B = Z.B(:, :, i); W = M .* B; x = Z.xstar(i, :)';
        for t = 1:30
            dW = zeros(4, 4, 5);
            for k = 1:5
                [~, bias, Bs] = dial(M, B, W, x, S(k)); [Mm, Bm] = apply_mut(M, Bs, Z.q(i, t), Z.tog(i, t), Z.d(i, t));
                dW(:, :, k) = Mm .* Bm - M .* Bs; dB = max(dB, max(abs(bias(1:2) - Z.bias(i, :, k)')));
            end
            if E.ex(i, t)
                nEx = nEx + 1; dE = max(dE, max(abs(dW - dW(:, :, 3) .* reshape(S, 1, 1, 5)), [], 'all'));
            else
                nKp = nKp + 1; dK = max(dK, max(abs(dW - dW(:, :, 3)), [], 'all'));
            end
        end
    end
    ok = dK <= 1e-12 && dE <= 1e-12;
    fprintf('\nCHECK 2 (same mutation at every setting; %d kept, %d excluded mutants, primary set): max |dW_s - dW_1| kept %.3g; max |dW_s - s dW_1| excluded %.3g -> %s\n', ...
        nKp, nEx, dK, dE, pf(ok));
    fprintf('  beside, not gating: recomputed bias against the stored bias, max |difference| %.3g\n', dB);
end


function ok = check3(K, P)
% Check 3, null calibration: settings permuted within each network, 200 seeds, 2,000 sign vectors
% per seed for all four contrasts; passes if each contrast has p < 0.05 in at most 7% of seeds.
    nS = 200; x = false(nS, 4); ks = [3 1; 5 2; 3 2; 5 3]; t = tic; [~, ~, loc] = unique(P.rp); G = max(loc);
    for i = 1:nS
        X = perm(P.L, K.seed + 70000 + i); s = flips(K.seed + 72000 + i, G);
        for j = 1:4, x(i, j) = flip_p(X(:, ks(j, 1)) - X(:, ks(j, 2)), loc, s) < 0.05; end
    end
    r = mean(x, 1); ok = all(r <= 0.07);
    fprintf('\nCHECK 3 (null calibration, %d seeds, 2,000 flips; %.0f s): share p < 0.05 STEP %.3f, SLOPE %.3f, LOWER %.3f, UPPER %.3f (each <= 0.07) -> %s\n', ...
        nS, toc(t), r, pf(ok));
end


function X = perm(Ls, seed)
% The five settings' columns permuted within each network, one randperm per network in order.
    rs = rng; rng(seed, 'twister'); X = Ls;
    for i = 1:size(Ls, 1), X(i, :) = Ls(i, randperm(5)); end
    rng(rs);
end


function s = flips(seed, G)
    rs = rng; rng(seed, 'twister'); s = 2 * (rand(G, 2000) < 0.5) - 1; rng(rs);
end


function p = flip_p(d, loc, s)
    Sr = accumarray(loc, d); T = abs(sum(Sr)); p = (1 + sum(abs(Sr' * s) >= T - 1e-12)) / (1 + size(s, 2));
end


% =========================================================================
% THE RECOUNT
% =========================================================================
function recount_lines(R, C)
% The four contrasts on the primary set with the mutants of the dialed connections left out, and
% the checks.
    fprintf('\nRECOUNT (primary set, toggles on the dialed connections excluded; Holm-adjusted p beside):\n');
    for j = 1:4
        fprintf('  %-5s %-14s %.4f [%.4f, %.4f]; sign-flip p %.4f (Holm %.4f)\n', R(j).name, R(j).def, R(j).e, R(j).ci, R(j).p, R(j).holm);
    end
    fprintf('  checks: 1 %s, 2 %s, 3 %s\n', pf(C.c1), pf(C.c2), pf(C.c3));
end


% =========================================================================
% ALSO REPORTED
% =========================================================================
function removed(K, D, E, P, R)
% Mutants left out, and the share of each contrast removed by leaving them out:
% 1 - recount / estimate with every mutant kept.
    p = D.prim; ex = E.ex(p, :); n = sum(ex, 2);
    fprintf('mutants excluded (primary set, %d networks): total %d of %d; per network median %g, range [%d, %d]; toggle off %d, toggle on %d; networks dropped (no kept mutant) %d\n', ...
        sum(p), sum(n), numel(ex), median(n), min(n), max(n), sum(E.off(p, :), 'all'), sum(E.on(p, :), 'all'), P.ndrop);
    fprintf('  (all networks in the set, %d: excluded %d)\n', sum(D.set), sum(E.ex(D.set, :), 'all'));
    fprintf('L_x by setting %s (primary, kept mutants): %s\n', mat2str(K.S), mat2str(mean(P.L, 1), 4));
    L0 = D.LR(p, :); o = [mean(L0(:, 3) - L0(:, 1)), mean(L0(:, 5) - L0(:, 2))];
    fprintf('share removed by the exclusion: STEP %.4f -> %.4f (%.3f removed); SLOPE %.4f -> %.4f (%.3f removed)\n', ...
        o(1), R(1).e, 1 - R(1).e / o(1), o(2), R(2).e, 1 - R(2).e / o(2));
end


function unsettled_too(K, D, E)
% Beside: also drop every mutant whose run from rest is not at rest by the cap at any setting.
    p = D.prim; keep = E.keep(p, :) & all(D.Z.rMOk(p, :, :), 3); nk = sum(keep, 2); u = nk > 0;
    L = Lx(D.lost(p, :, :), keep); rp = D.rep(p); R = contrasts4(K, L(u, :), rp(u));
    fprintf('beside (unsettled also dropped; %d kept mutants, %d networks): L_x %s\n', sum(keep(:)), sum(u), mat2str(mean(L(u, :), 1), 4));
    for j = 1:4
        fprintf('  %-5s %.4f [%.4f, %.4f]; sign-flip p %.4f\n', R(j).name, R(j).e, R(j).ci, R(j).p);
    end
end


function dev_rules(K, Hd, D, E)
% The two developmental rules at s = 0, 0.5, 1, with and without the left-out toggles. Set: kept
% at s = 1 and development test passed at the three settings; primary: fifths 2 to 5. MISMATCH is
% printed if the primary set is not 199 or L with every mutant kept differs from the fixed string
% (the value the stored response_dial_analysis.txt printed).
    Z = D.Z; I = D.I; GD = Hd.check_dev(K, Z, I.self);
    R = {'ENCOUNTER RULE', Z.eMut < 0.5, GD.keptE, '[0.02161 0.02513 0.07186]'; ...
         'CENSUS RULE', ~(Z.cMut > 0), GD.keptC, '[0.01173 0.01742 0.09447]'};
    for r = 1:2
        lost = R{r, 2}(:, :, K.dev); set = R{r, 3} & all(GD.pass(:, K.dev), 2); p = set & I.fifth >= 2; rp = D.rep(p);
        L0 = Lx(lost(p, :, :), true(sum(p), 30)); s0 = mat2str(mean(L0, 1), 4);
        mm = ~strcmp(s0, R{r, 4}) || sum(p) ~= 199;
        fprintf('%s (developmental, 0 / 0.5 / 1; primary %d): L no exclusion %s (response_dial_analysis.txt %s)%s\n', R{r, 1}, sum(p), s0, R{r, 4}, ...
            Hd.ternary(mm, ' MISMATCH', ''));
        keep = E.keep(p, :); u = sum(keep, 2) > 0; L1 = Lx(lost(p, :, :), keep);
        fprintf('  L excluded %s (%d networks)\n', mat2str(mean(L1(u, :), 1), 4), sum(u));
        dev_line('no exclusion', K, L0, rp); dev_line('excluded', K, L1(u, :), rp(u));
    end
end


function dev_line(tag, K, L, rp)
    d1 = L(:, 3) - L(:, 1); d2 = L(:, 3) - L(:, 2); [e1, c1] = cluster_t(d1, rp); [e2, c2] = cluster_t(d2, rp);
    fprintf('  %-12s STEP L(1) - L(0) %.4f [%.4f, %.4f] p %.4f; LOWER L(1) - L(0.5) %.4f [%.4f, %.4f] p %.4f\n', ...
        tag, e1, c1, signflip(K, d1, rp), e2, c2, signflip(K, d2, rp));
end


function alpha_table(K, D)
% alpha (alpha_and_beta, stored per setting) on the primary set, NaN omitted and counted: median,
% [25th, 75th] percentiles, mean, and the median per-network ratio alpha(s) / alpha(1); beside, the
% response slope and its median ratio slope(s) / slope(1), non-finite ratios dropped.
    p = D.prim; a = D.Z.alpha(p, :); sl = D.Z.slope(p, :);
    ra = a ./ a(:, K.i1); rs = sl ./ sl(:, K.i1); q = prctile(a, [25 75], 1);
    mr = @(x) arrayfun(@(k) median(x(isfinite(x(:, k)), k)), 1:5);
    fprintf('alpha = ||v|| (alpha_and_beta, stored per setting; primary %d; NaN omitted, count by setting %s) by setting %s: median %s; IQR %s; mean %s; median alpha(s)/alpha(1) %s (non-finite ratios dropped: %d)\n', ...
        sum(p), mat2str(sum(isnan(a), 1)), mat2str(K.S), mat2str(median(a, 1, 'omitnan'), 4), mat2str(q', 4), ...
        mat2str(mean(a, 1, 'omitnan'), 4), mat2str(mr(ra), 4), sum(~isfinite(ra(:))));
    fprintf('  beside: response slope median %s; median slope(s)/slope(1) %s (non-finite dropped: %d)\n', ...
        mat2str(median(sl, 1, 'omitnan'), 4), mat2str(mr(rs), 4), sum(~isfinite(rs(:))));
end


function decomp(Hd, D, E)
% The recounted STEP and SLOPE split by slot group and mutation type, each kept mutant weighted
% 1 / n_kept(n), so the groups of each split sum to the contrast.
    fprintf('\nWHERE THE RECOUNTED RISE COMES FROM (descriptive; kept mutants; 20,000-step classes)\n');
    Z = D.Z; p = D.prim; keep = E.keep(p, :); nk = sum(keep, 2); u = nk > 0; q = Z.q(p, :); tg = Z.tog(p, :);
    lo = D.lost(p, :, :); N = sum(u); pres = E.pres(p, :); w = keep(u, :) ./ nk(u);
    G = {'inputs to sensors (rows 3-4)', ismember(q, Hd.sensor_rows()); 'the four dialed connections', ismember(q, Hd.dialed_slots()); ...
        'inputs among C and D', ismember(q, Hd.cd_slots()); 'toggle off', tg & pres; 'toggle on', tg & ~pres; 'weight change', ~tg};
    C = {'STEP L(1) - L(0)', 3, 1; 'SLOPE L(2) - L(0.5)', 5, 2};
    for c = 1:2
        dl = (double(lo(u, :, C{c, 2})) - double(lo(u, :, C{c, 3}))) .* w; tot = sum(dl(:)) / N;
        fprintf('%s = %.4f (recount, %d networks)\n', C{c, 1}, tot, N);
        fprintf('  %-30s | part     | share  | mutants | + only higher | - only lower\n', 'group');
        for g = 1:size(G, 1)
            m = G{g, 2}(u, :) & keep(u, :); v = sum(dl(m)) / N;
            fprintf('  %-30s | %8.4f | %6.3f | %7d | %13d | %12d\n', G{g, 1}, v, v / tot, sum(m(:)), sum(dl(m) > 0), sum(dl(m) < 0));
            if g == 3 || g == 6
                fprintf('  %-30s | %8.4f |\n', '(sum of the three above)', sum(cellfun(@(x) sum(dl(x(u, :) & keep(u, :))), G(g - 2:g, 2))) / N);
            end
        end
    end
end
