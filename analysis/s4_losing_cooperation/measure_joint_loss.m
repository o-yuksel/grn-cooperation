function out = measure_joint_loss(cfg)
% MEASURE_JOINT_LOSS  When one mutation of a member makes a mutant that defects on it, does the
% member still cooperate with that mutant more often than if the two were independent? Each
% member against its 30 drawn mutants under the runs' encounter rule, and beside them, for a
% sample, against its 160 single-weight grid mutants.
% Produces : no number the paper reports; a supporting analysis for SI S4.4, how members answer their own
%            mutants: printed block [1], J, A, B and X = J / (A B), and block [2], the escape rate
%            e over the drawn mutants and over the grid mutants and their ratio S
% Inputs   : cfg (optional struct; defaults in fill_defaults; 'reportOnly', true reprints the
%              saved report; 'pilot', true runs replicates 1 to 3 only and saves nothing)
%            publication/rho_by_architecture.mat   (written by per_mutation_census: per member
%              its count of individuals, clone readout and 30 draws of the mutation distribution)
%            publication/<architecture>_100reps.mat   (written by run_signaling_architectures;
%              the members' wirings are recovered from them)
% Outputs  : publication/later_defector_quality/joint_loss.mat and .txt (the printed report)
% Run      : measure_joint_loss()   (about 2 min)
% Seed     : all 'twister', offsets from 20260927 (cfg.rngSeed): + 1000 m + r, the start states
%            of architecture m, replicate r; + 50 + m, the sampled members; + 70000 + 1000 m + q,
%            their start states; + 90 + m, the membership_test check; + 7 and + 8, the bootstraps
%
% per_mutation_census stores each member's mutations but not its wiring; recover_parents finds
% it in the stored snapshot (rows left tied are dropped and counted). Main set: that census's
% 'ok' members (cooperation margin at least 1e-6, size of the settled state ||x*|| in [0.5, 14]);
% second set: every recovered member. Each member meets each mutant in K = 20 encounters, member
% first, the same start pairs for all its mutants. Weighted by individuals: A = P(mutant D),
% B = P(member C), J = P(mutant D and member C), X = J / (A B). The escape rate J / A over the 30
% drawn mutants ('e_distr') is set beside the same over the 160 grid mutants ('e_grid') for up
% to 150 members per architecture; S is their ratio. Intervals: 95% bootstrap over replicates
% within each architecture.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
    matFile = fullfile(cfg.outDir, 'joint_loss.mat');
    if cfg.reportOnly
        S = load(matFile, 'out'); out = S.out;
    else
        t0 = tic; out = run_all(cfg); out.minutes = toc(t0) / 60;
        if cfg.pilot
            evalc('print_report(out, cfg)'); fprintf('timing run: %.1f min, report code ran\n', out.minutes);
            return;
        end
        save(matFile, 'out', 'cfg', '-v7.3');
    end
    txt = fullfile(cfg.outDir, 'joint_loss.txt');
    if exist(txt, 'file'), delete(txt); end
    diary(txt); print_report(out, cfg); diary off;
end


% =========================================================================
function cfg = fill_defaults(cfg)
    pub = fullfile(grn_data_file('publication'));
    def = struct('dataDir', pub, 'outDir', '', 'rngSeed', 20260927, 'K', 20, 'nBoot', 2000, ...
        'nGridPerArch', 150, 'nCssCheck', 50, 'pilot', false, 'reportOnly', false);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'later_defector_quality'); end
end


% =========================================================================
function out = run_all(cfg)
    R = load(fullfile(cfg.dataDir, 'rho_by_architecture.mat'), 'W_all', 'modeNames', 'metadata');
    rc = R.metadata.cfg;
    modes = arch_table();
    out.modes = R.modeNames; out.arch = cell(1, 4);
    for m = 1:4
        assert(strcmp(R.modeNames{m}, modes{m, 1}), 'mode order');
        rows = R.W_all{m}; R.W_all{m} = [];
        if cfg.pilot, rows = rows([rows.rep] <= 3); end
        [sp, cp] = params(modes(m, :));
        t0 = tic;
        [Wp, rec] = recover_parents(rows, modes{m, 1}, rc, cfg);
        A = encounter_all(rows, Wp, sp, cfg, m);
        A.recovery = rec;
        A.cssCheck = css_check(rows, Wp, cp, sp, cfg, m);
        A.grid = grid_arm(rows, Wp, sp, cfg, m);
        A.name = modes{m, 1}; A.minutes = toc(t0) / 60;
        out.arch{m} = A;
        fprintf('%s: %d wirings, %.1f min\n', modes{m, 1}, numel(rows), A.minutes);
    end
end


function T = arch_table()
    T = {'graded_separate',     'separate',     Inf
         'graded_contrast',     'contrast', Inf
         'saturating_separate', 'separate',     1
         'saturating_contrast', 'contrast', 1};
end


function [sp, cp] = params(t)
% sp: the development parameters of the runs and of simulate_pair (as per_mutation_census's
% sim_params); cp: membership_test's (as its params).
    sp = struct('sensingMode', t{2}, 'hillCoefficient', t{3}, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), ...
        'externalSignal2', zeros(4, 1));
    cp = struct('sensingMode', t{2}, 'hillCoefficient', t{3}, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, ...
        'developmentSteps', 600, 'backend', 'simulate_pair');
end


% =========================================================================
function [Wp, rec] = recover_parents(rows, mode, rc, cfg)
% Parent wiring of every census row: per_mutation_census's snapshot selection, then a
% match on clone count and clone readout in ascending wiring order. Check (a).
    d = load(fullfile(cfg.dataDir, sprintf('%s_100reps.mat', mode)), 'results');
    sp = params_by_name(mode);
    nR = numel(rows); Wp = nan(nR, 16);
    rep = [rows.rep]; gen = [rows.gen];
    nAmb = 0; nMiss = 0; nTie = 0; maxDr = 0;
    for r = unique(rep)
        R = d.results{r};
        sg = double(R.historyGenerations(:))';
        ks = find(sg >= rc.steadyStateFrom);
        ks = ks(round(linspace(1, numel(ks), rc.nSnaps)));
        for k = ks
            idx = find(rep == r & gen == sg(k));
            if isempty(idx), continue; end
            M = double(R.historyM(:,:,:,k)); B = double(R.historyB(:,:,:,k)); N = size(M, 3);
            [Wu, ~, ic] = unique(reshape(M .* B, 16, N)', 'rows');
            counts = accumarray(ic, 1);
            last = 0;
            for i = idx
                cand = find(counts == rows(i).count & (1:numel(counts))' > last);
                dr = inf(numel(cand), 1);
                for c = 1:numel(cand)
                    dr(c) = abs(clone_r(Wu(cand(c), :), sp) - rows(i).r);
                end
                hit = cand(dr < 1e-9);
                if isempty(hit), nMiss = nMiss + 1; continue; end
                maxDr = max(maxDr, min(dr(dr < 1e-9)));
                if numel(hit) > 1                       % a tie: the draws, then membership
                    nTie = nTie + 1;
                    hit = hit(arrayfun(@(h) draws_consistent(Wu(h, :), rows(i)), hit));
                    if numel(hit) > 1
                        hit = hit(arrayfun(@(h) is_member_emu(Wu(h, :), sp), hit));
                    end
                    if numel(hit) ~= 1
                        nAmb = nAmb + 1; last = min([cand(dr < 1e-9); Inf]) ; continue;
                    end
                end
                Wp(i, :) = Wu(hit, :); last = hit;
            end
        end
    end
    rec = struct('n', nR, 'missing', nMiss, 'ambiguous', nAmb, 'ties', nTie, 'maxDr', maxDr, ...
        'toggleBad', toggle_check(rows, Wp));
end


function ok = draws_consistent(w, row)
% A candidate parent must explain every recorded draw: a toggle on a nonzero
% weight sets it to zero (dW = -W); a weight step on a zero weight (M = 0) changes nothing
% (dW = 0), and on a nonzero weight changes it (dW ~= 0).
    s = row.slot(:); t = row.isToggle(:); dW = row.dW(:); ws = w(s)'; ws = ws(:);
    okT = ~t | ws == 0 | abs(dW + ws) < 1e-12;
    okW = t | ((ws == 0) == (dW == 0));
    ok = all(okT & okW);
end


function m = is_member_emu(w, sp)
% The last tie-breaker: membership by the batch emulation of the membership test (fixed start
% 0.005, 600 steps; C1, and no grid mutant defects on the wiring while it cooperates).
    W = reshape(w, 4, 4); D = grid_deltas(); nN = size(D, 1); X0 = 0.005 * ones(4, 1);
    Wn = repmat(W, 1, 1, nN);
    for q = 1:nN, Wn(D(q,1), D(q,2), q) = W(D(q,1), D(q,2)) + D(q,3); end
    [a1, a2] = encounter_batch(repmat(W, 1, 1, nN + 1), cat(3, W, Wn), repmat(X0, 1, nN + 1), ...
        repmat(X0, 1, nN + 1), sp, 600, 'simpair');
    m = a1(1) && a2(1) && ~any(a1(2:end) & ~a2(2:end));
end


function sp = params_by_name(mode)
    T = arch_table(); sp = params(T(strcmp(T(:, 1), mode), :));
end


function r = clone_r(w, sp)
    W = reshape(w, 4, 4); E = ones(4);
    [~, fs, ~] = simulate_pair(E, W, [], E, W, [], sp);
    r = fs(1) - fs(2);
end


function nBad = toggle_check(rows, Wp)
% Control (b): a toggle on a nonzero weight must set it to zero.
    nBad = 0;
    for i = 1:numel(rows)
        if any(isnan(Wp(i, :))), continue; end
        t = rows(i).isToggle; s = rows(i).slot; dW = rows(i).dW;
        live = t & Wp(i, s(:))' ~= 0;
        nBad = nBad + sum(abs(dW(live) + Wp(i, s(live))') > 1e-12);
    end
end


% =========================================================================
function A = encounter_all(rows, Wp, sp, cfg, m)
% Parent against each of its 30 mutants, K draws under the runs' encounter rule, the same start
% pairs for all of a parent's mutants.
    nR = numel(rows); nK = numel(rows(1).dW); K = cfg.K;
    J = nan(nR, nK); Am = nan(nR, nK); Bp = nan(nR, nK);
    rep = [rows.rep];
    for r = unique(rep)
        rng(cfg.rngSeed + 1000 * m + r, 'twister');
        idx = find(rep == r);
        idx = idx(all(isfinite(Wp(idx, :)), 2));
        if isempty(idx), continue; end
        nP = numel(idx) * nK;
        W1 = zeros(4, 4, nP * K); W2 = W1; X1 = zeros(4, nP * K); X2 = X1;
        p = 0;
        for i = idx
            W = reshape(Wp(i, :), 4, 4);
            x1 = rand(4, K) * 0.01; x2 = rand(4, K) * 0.01;
            for t = 1:nK
                Wm = W; Wm(rows(i).slot(t)) = W(rows(i).slot(t)) + rows(i).dW(t);
                cols = p * K + (1:K); p = p + 1;
                W1(:, :, cols) = repmat(W, 1, 1, K); W2(:, :, cols) = repmat(Wm, 1, 1, K);
                X1(:, cols) = x1; X2(:, cols) = x2;
            end
        end
        [a1, a2] = encounter_batch(W1, W2, X1, X2, sp, 300, 'encounter_rule');
        a1 = reshape(a1, K, nK, numel(idx)); a2 = reshape(a2, K, nK, numel(idx));
        J(idx, :) = squeeze(mean(a1 & ~a2, 1))';
        Am(idx, :) = squeeze(mean(~a2, 1))';
        Bp(idx, :) = squeeze(mean(a1, 1))';
    end
    A = struct('J', J, 'A', Am, 'B', Bp, 'rep', rep(:), 'w', [rows.count]', 'ok', [rows.ok]', ...
        'isToggle', [rows.isToggle]', 'absdW', abs([rows.dW]'), 'recovered', all(isfinite(Wp), 2));
end


% =========================================================================
function c = css_check(rows, Wp, cp, sp, cfg, m)
% Check (c): on a sample of recovered mutants, membership_test must reproduce the stored 'ejected' flag
% (the mutation ends membership) and simulate_pair's clone readout the 'lostReadout' flag (it ends
% cooperation).
    rng(cfg.rngSeed + 90 + m, 'twister');
    ok = find(all(isfinite(Wp), 2)); nK = numel(rows(1).dW);
    pick = ok(randi(numel(ok), cfg.nCssCheck, 1)); tt = randi(nK, cfg.nCssCheck, 1);
    mE = 0; mR = 0; E = ones(4);
    for q = 1:cfg.nCssCheck
        i = pick(q); t = tt(q);
        W = reshape(Wp(i, :), 4, 4); Wm = W; Wm(rows(i).slot(t)) = W(rows(i).slot(t)) + rows(i).dW(t);
        mE = mE + ((~membership_test(Wm, cp).is_css) ~= rows(i).ejected(t));
        [~, fm, ~] = simulate_pair(E, Wm, [], E, Wm, [], sp);
        mR = mR + ((~((fm(1) - fm(2)) > 0)) ~= rows(i).lostReadout(t));
    end
    c = struct('n', cfg.nCssCheck, 'endsMembershipMismatch', mE, 'readoutMismatch', mR);
end


% =========================================================================
function G = grid_arm(rows, Wp, sp, cfg, m)
% The grid arm: the escape rate over the 30 drawn mutants ('e_distr') against the same over the
% 160 grid mutants ('e_grid'), on a seeded subsample of 'ok' parents, the same start pairs for both.
    rng(cfg.rngSeed + 50 + m, 'twister');
    ok = find([rows.ok]' & all(isfinite(Wp), 2));
    sub = sort(ok(randperm(numel(ok), min(cfg.nGridPerArch, numel(ok)))));
    D = grid_deltas(); nN = size(D, 1); K = cfg.K; nK = numel(rows(1).dW);
    G = struct('rows', sub, 'rep', [rows(sub).rep]', 'w', [rows(sub).count]', ...
        'kJ', nan(numel(sub), 1), 'kA', nan(numel(sub), 1), 'gJ', nan(numel(sub), 1), ...
        'gA', nan(numel(sub), 1), 'gJcls', nan(numel(sub), 1), 'gAcls', nan(numel(sub), 1));
    for q = 1:numel(sub)
        i = sub(q);
        rng(cfg.rngSeed + 70000 + 1000 * m + q, 'twister');
        W = reshape(Wp(i, :), 4, 4);
        x1 = rand(4, K) * 0.01; x2 = rand(4, K) * 0.01;
        Wn = zeros(4, 4, nN + nK);
        for n = 1:nN, Wn(:, :, n) = W; Wn(D(n,1), D(n,2), n) = W(D(n,1), D(n,2)) + D(n,3); end
        for t = 1:nK
            Wm = W; Wm(rows(i).slot(t)) = W(rows(i).slot(t)) + rows(i).dW(t); Wn(:, :, nN + t) = Wm;
        end
        nO = nN + nK;
        [a1, a2] = encounter_batch(repmat(W, 1, 1, nO * K), Wn(:, :, repelem(1:nO, K)), ...
            repmat(x1, 1, nO), repmat(x2, 1, nO), sp, 300, 'encounter_rule');
        a1 = reshape(a1, K, nO); a2 = reshape(a2, K, nO);
        jj = mean(a1 & ~a2, 1); aa = mean(~a2, 1);
        G.gJ(q) = mean(jj(1:nN)); G.gA(q) = mean(aa(1:nN));
        G.kJ(q) = mean(jj(nN+1:end)); G.kA(q) = mean(aa(nN+1:end));
        X0 = 0.005 * ones(4, 1);                 % grid mutants that defect from the membership test's start
        [~, c2] = encounter_batch(repmat(W, 1, 1, nN), Wn(:, :, 1:nN), repmat(X0, 1, nN), ...
            repmat(X0, 1, nN), sp, 600, 'simpair');
        dc = ~c2;
        G.gJcls(q) = sum(jj(dc)) / nN; G.gAcls(q) = sum(aa(dc)) / nN;
    end
end


function D = grid_deltas()
% membership_test's grid: linspace(-3*0.2, 3*0.2, 11) with the exact zero skipped.
    d = linspace(-0.6, 0.6, 11); d = d(d ~= 0);
    [dd, jj, ii] = ndgrid(d, 1:4, 1:4);
    D = [ii(:), jj(:), dd(:)];
end


% =========================================================================
% REPORT
% =========================================================================
function print_report(out, cfg)
    fprintf('JOINT LOSS: does one mutation from the mutation distribution make a defector its member parent does not punish?\n');
    fprintf('K = %d draws per pairing under the runs'' encounter rule; nBoot %d; %.1f min.\n', ...
        cfg.K, cfg.nBoot, out.minutes);
    print_controls(out);
    print_main_measure(out, cfg);
    print_other_measures(out, cfg);
    print_split(out);
end


% Checks per architecture: (a) no row unmatched and at most 1% of rows left tied; (b) no toggle
% draw on a nonzero parent weight with dW ~= -W; (c) no mismatch in the membership_test sample.
function print_controls(out)
    fprintf('\n[0] CONTROLS\n');
    for m = 1:4
        A = out.arch{m}; r = A.recovery; c = A.cssCheck;
        o1 = r.missing == 0 && r.ambiguous <= 0.01 * r.n;
        o2 = r.toggleBad == 0; o3 = c.endsMembershipMismatch == 0 && c.readoutMismatch == 0;
        fprintf('  %-20s (a) recovered %d/%d rows (%d count+readout ties resolved by the draws, %d unresolved and excluded), max|dr| %.1e -> %s | (b) toggle mismatches %d -> %s | (c) membership_test on %d mutants: mismatches in ending membership %d, readout %d -> %s\n', ...
            A.name, r.n - r.missing - r.ambiguous, r.n, r.ties - r.ambiguous, r.ambiguous, r.maxDr, pf(o1), r.toggleBad, pf(o2), ...
            c.n, c.endsMembershipMismatch, c.readoutMismatch, pf(o3));
    end
end


function print_main_measure(out, cfg)
    fprintf('\n[1] J = P(mutant D & parent C), A = P(mutant D), B = P(parent C); X = J/(A*B)\n');
    fprintf('    clone-weighted over draws; 95%% cluster bootstrap over replicates within architecture\n');
    fprintf('  %-20s %-5s %7s %8s | %8s %8s %8s | %24s | %8s | %8s\n', 'architecture', 'set', 'wirings', 'draws', ...
        'J', 'A', 'B', 'X = J/(A*B) [95% CI]', 'ceil 1/B', 'J/A');
    for set = 1:2
        for m = [1:4, 0]
            [s, nW, nD] = main_stats(out, m, set, cfg);
            fprintf('  %-20s %-5s %7d %8d | %8.5f %8.5f %8.5f | %.3f [%.3f, %.3f]     | %8.3f | %8.4f\n', ...
                mname(out, m), setname(set), nW, nD, s.J, s.A, s.B, s.X, s.Xci, 1 / s.B, s.J / s.A);
        end
    end
end


function [s, nW, nD] = main_stats(out, m, set, cfg)
    if m == 0, ms = 1:4; else, ms = m; end
    parts = cell(numel(ms), 1);
    for q = 1:numel(ms)
        A = out.arch{ms(q)};
        keep = A.recovered & (set == 2 | A.ok);
        parts{q} = struct('J', A.J(keep, :), 'A', A.A(keep, :), 'B', A.B(keep, :), 'w', A.w(keep), 'rep', A.rep(keep));
    end
    f = @(P) ratio_of(P);
    s = f(parts); nW = sum(cellfun(@(p) numel(p.w), parts)); nD = sum(cellfun(@(p) numel(p.J), parts));
    rng(cfg.rngSeed + 7, 'twister');
    Xb = nan(cfg.nBoot, 1);
    for b = 1:cfg.nBoot, Xb(b) = getfield(f(resample(parts)), 'X'); end %#ok<GFLD>
    s.Xci = quantile(Xb, [0.025 0.975]);
end


function s = ratio_of(parts)
    sw = 0; sJ = 0; sA = 0; sB = 0;
    for q = 1:numel(parts)
        p = parts{q}; w = p.w;
        sw = sw + sum(w) * size(p.J, 2);
        sJ = sJ + sum(w .* sum(p.J, 2)); sA = sA + sum(w .* sum(p.A, 2)); sB = sB + sum(w .* sum(p.B, 2));
    end
    s.J = sJ / sw; s.A = sA / sw; s.B = sB / sw; s.X = s.J / (s.A * s.B);
end


function parts = resample(parts)
% Cluster bootstrap: replicates resampled with replacement within each architecture.
    for q = 1:numel(parts)
        p = parts{q}; u = unique(p.rep); pick = u(randi(numel(u), numel(u), 1));
        idx = cell2mat(arrayfun(@(r) find(p.rep == r), pick, 'uni', 0));
        parts{q} = struct('J', p.J(idx, :), 'A', p.A(idx, :), 'B', p.B(idx, :), 'w', p.w(idx), 'rep', p.rep(idx));
    end
end


function print_other_measures(out, cfg)
    fprintf('\n[2] also reported: escape e = P(parent C | opponent D), mutants from the mutation distribution (distr) vs the parent''s 160 grid mutants\n');
    fprintf('    same %d-parent-per-architecture subsample, same start pairs; S = e_distr / e_grid, cluster bootstrap\n', cfg.nGridPerArch);
    fprintf('  %-20s %5s | %9s %9s | %9s %9s | %24s | %9s\n', 'architecture', 'n', 'e_distr', 'A_distr', 'e_grid', ...
        'A_grid', 'S [95% CI]', 'e_gridCls');
    for m = [1:4, 0]
        s = other_stats(out, m, cfg);
        fprintf('  %-20s %5d | %9.4f %9.4f | %9.4f %9.4f | %.3f [%.3f, %.3f]     | %9.4f\n', mname(out, m), s.n, ...
            s.ek, s.Ak, s.eg, s.Ag, s.S, s.Sci, s.egc);
    end
    fprintf('  (e_gridCls: escape against only the grid mutants that defect under the membership test)\n');
end


function s = other_stats(out, m, cfg)
    if m == 0, ms = 1:4; else, ms = m; end
    G = cellfun(@(A) A.grid, out.arch(ms), 'uni', 0);
    f = @(G) sec_of(G);
    s = f(G); s.n = sum(cellfun(@(g) numel(g.w), G));
    rng(cfg.rngSeed + 8, 'twister');
    Sb = nan(cfg.nBoot, 1);
    for b = 1:cfg.nBoot
        Gb = G;
        for q = 1:numel(Gb)
            g = Gb{q}; u = unique(g.rep); pick = u(randi(numel(u), numel(u), 1));
            idx = cell2mat(arrayfun(@(r) find(g.rep == r), pick, 'uni', 0));
            Gb{q} = struct('w', g.w(idx), 'kJ', g.kJ(idx), 'kA', g.kA(idx), 'gJ', g.gJ(idx), 'gA', g.gA(idx), ...
                'gJcls', g.gJcls(idx), 'gAcls', g.gAcls(idx), 'rep', g.rep(idx));
        end
        Sb(b) = getfield(f(Gb), 'S'); %#ok<GFLD>
    end
    s.Sci = quantile(Sb, [0.025 0.975]);
end


function s = sec_of(G)
    w = []; kJ = []; kA = []; gJ = []; gA = []; cJ = []; cA = [];
    for q = 1:numel(G)
        g = G{q}; w = [w; g.w]; kJ = [kJ; g.kJ]; kA = [kA; g.kA]; gJ = [gJ; g.gJ]; gA = [gA; g.gA]; %#ok<AGROW>
        cJ = [cJ; g.gJcls]; cA = [cA; g.gAcls]; %#ok<AGROW>
    end
    s.ek = (w' * kJ) / (w' * kA); s.eg = (w' * gJ) / (w' * gA); s.S = s.ek / s.eg;
    s.Ak = (w' * kA) / sum(w); s.Ag = (w' * gA) / sum(w); s.egc = (w' * cJ) / max(w' * cA, eps);
end


function print_split(out)
    fprintf('\n[3] J PER MUTATION BY MUTATION TYPE (primary set, pooled, clone-weighted)\n');
    fprintf('  %-26s %9s %9s %9s %9s\n', 'type', 'share', 'J', 'A', 'J/A');
    labs = {'toggle', 'weight |dW| <= 0.6', 'weight |dW| > 0.6'};
    for c = 1:3
        [sw, sJ, sA, sT] = deal(0);
        for m = 1:4
            A = out.arch{m}; keep = A.recovered & A.ok;
            T = A.isToggle(keep, :); dW = A.absdW(keep, :); w = repmat(A.w(keep), 1, size(T, 2));
            switch c
                case 1, sel = T;
                case 2, sel = ~T & dW <= 0.6;
                case 3, sel = ~T & dW > 0.6;
            end
            Jm = A.J(keep, :); Am = A.A(keep, :);
            sw = sw + sum(w(sel)); sJ = sJ + sum(w(sel) .* Jm(sel)); sA = sA + sum(w(sel) .* Am(sel));
            sT = sT + sum(w(:));
        end
        fprintf('  %-26s %9.4f %9.5f %9.5f %9.4f\n', labs{c}, sw / sT, sJ / sw, sA / sw, sJ / sA);
    end
end


function s = mname(out, m)
    if m == 0, s = 'POOLED'; else, s = out.arch{m}.name; end
end


function s = setname(set)
    if set == 1, s = 'ok'; else, s = 'all'; end
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
