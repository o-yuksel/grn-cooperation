function tipping_input_flag(stage, nPerArch)
% TIPPING_INPUT_FLAG  How many of the mutations that end a member's cooperation or retaliation does
% the tipping input catch? Finds each member's tipping inputs and flags a mutation when the input
% it adds passes the one of its gene and sign.
%
% Produces : main text, "Conditional cooperation stays within a mutation's reach of the
%            threshold"; Methods, "Mutation tests in the game"; Fig. 4A (the input that
%            fig4_tipping_input reads); SI S3.3, "The search" and "Treating a mutation as the
%            input it adds catches most losses" (printed by 'analyze')
% Inputs   : the members of sample_census_members (from publication/retaliation_fixed_opponent/
%              retaliation_fixed_opponent_census.mat, written by measure_retaliation_fixed_opponent)
%            'analyze': publication/tipping_input/tipping_input_<nPerArch>.mat (from 'run')
% Outputs  : 'run': publication/tipping_input/tipping_input_<nPerArch>.mat (S, the members; R, per
%              member the tipping inputs thrC and thrR and per draw the re-simulated outcomes, the
%              source gene's expression and the flags; timings; seed)
%            'analyze': prints the check of the re-simulated outcomes, the share of tipping inputs
%              not found, recall and precision with intervals, and per-member readings
% Run      : tipping_input_flag('run', 7000), then tipping_input_flag('analyze', 7000) (7000
%              exceeds every architecture's count, so every member with a worst-case mutant is
%              scored; 2 workers; long, not timed here: 3 members per architecture take about 30 s)
% Seed     : the simulations draw nothing; the sample's base seed is 20261016 (rng(20261016 + m,
%            'twister') for architecture m, in sample_census_members); the bootstraps use
%            rng(20261016 + 9, 'twister') pooled and 20261016 + 9 + m per architecture, around a
%            saved and restored stream
%
% Tipping input at gene i, of each sign: the smallest constant input that, added at every step
% from the membership test's fixed start (every gene at 0.005, 600 steps), ends cooperation (added
% in both networks of the clone pair) or retaliation (added to the member alone, against its
% worst-case mutant H); Inf when no input up to 10 does. Losses are the census's stored flags;
% 'analyze' first checks that re-simulating every mutant gives the same flags.

    O = grn_data_file('publication/tipping_input');
    if ~isfolder(O), mkdir(O); end
    f = fullfile(O, sprintf('tipping_input_%d.mat', nPerArch));
    switch stage
        case 'run', do_run(f, nPerArch);
        case 'analyze', do_analyze(f);
    end
end


function sp = sim_params(m)
% Architecture m = 1..4: graded separate, graded contrast, saturating separate, saturating contrast
    T = {'separate', Inf; 'contrast', Inf; 'separate', 1; 'contrast', 1};
    sp = struct('sensingMode', T{m,1}, 'hillCoefficient', T{m,2}, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), ...
        'externalSignal2', zeros(4, 1));
end


% =========================================================================
function do_run(f, n)
    S = sample_census_members(n, 20261016);
    start_pool();
    t0 = tic; R = cell(numel(S), 1); sec = nan(numel(S), 1);
    parfor k = 1:numel(S)
        tk = tic; R{k} = one_member(S(k), sim_params(S(k).arch)); sec(k) = toc(tk);
    end
    out = struct('S', S, 'R', {R}, 'sec', sec, 'wallSec', toc(t0), 'nPerArch', n, 'seed', 20261016);
    save(f, '-struct', 'out', '-v7.3');
    fprintf('tipping input: %d members, wall %.2f min, worker seconds %.0f\n', numel(S), out.wallSec / 60, sum(sec));
end


function o = one_member(s, sp)
% One member: its states with its clone (xc) and against its worst-case mutant H (xh), its
% tipping inputs at each gene and sign, and per recorded draw the re-simulated mutant, the source
% gene's expression x_j in the member and in the mutant, and the flag. The mutant is rebuilt as
% W with dW added at the drawn slot, as the census stored it.
    E = ones(4); W = s.W; K = numel(s.dW);
    H = W; H(s.hSlot) = H(s.hSlot) + s.hDelta;
    [~, fc] = simulate_pair(E, W, [], E, W, [], sp); xc = fc(1:4);
    [~, fh] = simulate_pair(E, W, [], E, H, [], sp); xh = fh(1:4);
    o.thrC = nan(4, 2); o.thrR = nan(4, 2);              % columns: sign +, sign -
    for i = 1:4
        for sgn = 1:2
            sv = 3 - 2 * sgn;                             % +1, -1
            o.thrC(i, sgn) = flip_threshold(@(d) flip_coop(W, sp, i, sv * d));
            o.thrR(i, sgn) = flip_threshold(@(d) flip_ret(W, H, sp, i, sv * d));
        end
    end
    o.lostRe = false(K, 1); o.retRe = false(K, 1); o.xjC = nan(K, 1); o.xjR = nan(K, 1);
    o.xjCm = nan(K, 1); o.xjRm = nan(K, 1); o.predC = false(K, 1); o.predR = false(K, 1);
    for t = 1:K
        q = s.slot(t); [i, j] = ind2sub([4 4], q); dl = s.dW(t);
        Wm = W; Wm(q) = W(q) + dl;
        [~, fm] = simulate_pair(E, Wm, [], E, Wm, [], sp); o.lostRe(t) = ~((fm(1) - fm(2)) > 0);
        [~, fmh] = simulate_pair(E, Wm, [], E, H, [], sp); o.retRe(t) = (fmh(1) - fmh(2)) > 0;
        o.xjC(t) = xc(j); o.xjR(t) = xh(j); o.xjCm(t) = fm(j); o.xjRm(t) = fmh(j);
        o.predC(t) = passes(dl * xc(j), o.thrC(i, :));
        o.predR(t) = passes(dl * xh(j), o.thrR(i, :));
    end
end


function p = passes(e, thr)
% A draw is flagged when the input it adds, e = x_j dW_ij, reaches the tipping input of its sign
% at gene i (thr(1) for e > 0, thr(2) for e < 0); e = 0 is never flagged
    if e == 0 || ~isfinite(e), p = false; return; end
    if e > 0, p = e >= thr(1); else, p = -e >= thr(2); end
end


function start_pool()
% A pool of 2 single-thread workers; on Windows at below-normal process priority (PowerShell),
% on other systems at the default priority
    cmd = 'powershell -NoProfile -Command "(Get-Process -Id %d).PriorityClass = ''BelowNormal''"';
    if ispc, system(sprintf(cmd, feature('getpid'))); end
    if isempty(gcp('nocreate')), parpool('Processes', 2); end
    spmd
        maxNumCompThreads(1);
        if ispc, system(sprintf(cmd, feature('getpid'))); end
    end
end


% =========================================================================
function do_analyze(f)
    L = load(f); S = L.S; R = L.R; n = numel(S); K = numel(S(1).dW);
    names = {'graded_separate','graded_contrast','saturating_separate','saturating_contrast'};
    ar = [S.arch]'; w = [S.count]';
    get = @(fn) cell2mat(cellfun(@(o) o.(fn)(:)', R, 'UniformOutput', false));
    lost = cell2mat(arrayfun(@(s) s.lost(:)', S, 'UniformOutput', false));
    retC = cell2mat(arrayfun(@(s) (s.rH(:) > 0)', S, 'UniformOutput', false));
    tog = cell2mat(arrayfun(@(s) s.isToggle(:)', S, 'UniformOutput', false));
    lostRe = get('lostRe'); retRe = get('retRe'); predC = get('predC'); predR = get('predR');
    fprintf('[SAMPLE] %d members (per architecture %s), seed %d, K = %d; wall %.2f min, worker seconds %.0f\n', ...
        n, mat2str(accumarray(ar, 1, [4 1])'), L.seed, K, L.wallSec / 60, sum(L.sec));
    % Check: the re-simulated mutants give the census's stored loss flags on every draw;
    % otherwise nothing further is read
    nC = sum(lostRe(:) ~= lost(:)); nR = sum(retRe(:) ~= retC(:));
    fprintf('[CONTROL] re-simulated mutants: cooperation lost (clone readout <= 0) vs res.lost %d disagreements; retaliation vs H (R(M''vs H) > 0) vs res.rH > 0 %d disagreements, over %d draws -> %s\n', ...
        nC, nR, numel(lost), pf(nC == 0 && nR == 0));
    if nC > 0 || nR > 0, fprintf('CONTROL FAILED: no tipping-input flag number is read.\n'); return; end
    thrC = cell2mat(cellfun(@(o) o.thrC(:)', R, 'UniformOutput', false));
    thrR = cell2mat(cellfun(@(o) o.thrR(:)', R, 'UniformOutput', false));
    fprintf('  thresholds Inf (no flip within |delta| <= 10): cooperation %.3f, retaliation %.3f of gene-sign cells\n', ...
        mean(isinf(thrC(:))), mean(isinf(thrR(:))));

    % Recall: the share of losses that are flagged; precision: the share of flagged draws that
    % are losses; both weighted by how many individuals carry the wiring (SI S3.3, "Treating a
    % mutation as the input it adds catches most losses").
    fprintf('\nRecall P(predicted | lost) and precision P(lost | predicted), clone-weighted; CI 2,000 member bootstraps within architecture (seed 20261016 + 9)\n');
    [rc, pc] = rp_boot(predC, lost, w, ar, 2000, 20261016 + 9);
    [rr, pr] = rp_boot(predR, retC, w, ar, 2000, 20261016 + 9);
    print_rp('cooperation', rc, pc, predC, lost, w, ar, names);
    print_rp('retaliation vs H', rr, pr, predR, retC, w, ar, names);

    fprintf('\n[also reported a] per member: predicted vs measured loss rate (fraction of its %d draws)\n', K);
    sides = {'cooperation', predC, lost; 'retaliation vs H', predR, retC};
    for s = 1:2
        pm = mean(sides{s, 2}, 2); mm = mean(sides{s, 3}, 2);
        rs = arrayfun(@(m) corr(pm(ar == m), mm(ar == m), 'Type', 'Spearman', 'Rows', 'complete'), 1:4);
        fprintf('  %-17s within-architecture Spearman %s; share of members with measured >= predicted %.3f (clone-weighted %.3f)\n', ...
            sides{s, 1}, mat2str(rs, 3), mean(mm >= pm), sum(w .* (mm >= pm)) / sum(w));
    end
    fprintf('\n[also reported b] misses (actual loss, not predicted): share on toggle draws and share where x_j changes by more than half in the mutant (same encounter), clone-weighted\n');
    W2 = repmat(w, 1, K);
    xc = get('xjC'); xcm = get('xjCm'); xr = get('xjR'); xrm = get('xjRm');
    bigC = abs(xcm - xc) > 0.5 * abs(xc); bigR = abs(xrm - xr) > 0.5 * abs(xr);
    mi = {lost & ~predC, bigC, 'cooperation'; retC & ~predR, bigR, 'retaliation vs H'};
    for s = 1:2
        M = mi{s, 1};
        fprintf('  %-17s misses %d draws: toggle %.3f (all draws %.3f); x_j moves > 1/2 %.3f (all draws %.3f)\n', mi{s, 3}, ...
            sum(M(:)), sum(W2(M) .* tog(M)) / sum(W2(M)), sum(W2(:) .* tog(:)) / sum(W2(:)), ...
            sum(W2(M) .* mi{s, 2}(M)) / sum(W2(M)), sum(W2(:) .* mi{s, 2}(:)) / sum(W2(:)));
    end
end


function [r, p] = rp_boot(P, Y, w, ar, B, seed)
% Weighted recall and precision of flags P against losses Y, with 95% percentile intervals over
% B bootstraps that resample members within each architecture
    W = repmat(w, 1, size(P, 2));
    a = sum(W .* (P & Y), 2); y = sum(W .* Y, 2); q = sum(W .* P, 2);
    r.p = sum(a) / sum(y); p.p = sum(a) / sum(q);
    rs = rng; rng(seed, 'twister'); bs = nan(B, 2);
    for b = 1:B
        c = zeros(numel(w), 1);
        for m = unique(ar)'
            ix = find(ar == m); c = c + accumarray(ix(randi(numel(ix), numel(ix), 1)), 1, [numel(w) 1]);
        end
        bs(b, :) = [c' * a / (c' * y), c' * a / (c' * q)];
    end
    rng(rs);
    r.ci = prctile(bs(:, 1), [2.5 97.5]); p.ci = prctile(bs(:, 2), [2.5 97.5]);
end


function print_rp(lab, r, p, P, Y, w, ar, names)
    fprintf('  %-17s pooled recall %.4f [%.4f, %.4f] | precision %.4f [%.4f, %.4f] | losses %d, predicted %d (draws)\n', ...
        lab, r.p, r.ci, p.p, p.ci, sum(Y(:)), sum(P(:)));
    for m = 1:4
        k = ar == m;
        [r1, p1] = rp_boot(P(k, :), Y(k, :), w(k), ar(k), 2000, 20261016 + 9 + m);
        fprintf('      %-20s recall %.4f [%.4f, %.4f] | precision %.4f [%.4f, %.4f] | losses %d\n', ...
            names{m}, r1.p, r1.ci, p1.p, p1.ci, sum(sum(Y(k, :))));
    end
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
