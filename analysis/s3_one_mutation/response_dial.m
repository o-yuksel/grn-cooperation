function varargout = response_dial(stage, nWorkers)
% RESPONSE_DIAL  The response dial of SI S3.5: in each of the 500 networks of
% response_dial_networks, the four connections from the sensor genes to C and D are scaled by
% s = 0, 0.5, 1, 1.5 or 2, while a constant input on C and D keeps the resting state with the
% clone, x*, unchanged. For each network and its 30 mutants at each setting it stores the outcomes
% every later reading uses.
%
% Produces : as input only, for main text, "Conditional cooperation stays within a mutation's
%            reach of the threshold"; Methods, "Mutation tests in the game"; Fig. 4; SI S3.4;
%            SI S3.5; Table S8 (response_dial_inputs.mat and chunks/, read through 'handles' by
%            dial_loss_shares, fig4b_retaliation_from_rest and the response_dial_* files)
% Inputs   : stage: 'checks', 'timing_run' (also accepted as 'pilot'), 'full' or 'handles'; nWorkers (default 3), the pool size
%            publication/conditionality_score/natural_partners/pool.mat and chunks/
%              (networks, mutations, reactions; response_dial_networks)
%            publication/<architecture>_100reps.mat   (saturating contrast's run config and the
%              networks' M and B; the other three are loaded unused; run_signaling_architectures)
% Outputs  : under publication/conditionality_score/response_dial/: 'checks' response_dial_inputs.mat (inputs
%            I, membership-test class cls), checks_a2.mat (pass flags), log_checks.txt; 'timing_run'
%            timing_run_a1/, timing_run_manifest.mat, timing_run_losses.txt, log_timing_run.txt; 'full'
%            chunks/chunk_<k>.mat (25 networks each), log_full.txt
% Run      : response_dial('checks')   builds the inputs, runs the membership test on each network
%              and checks 1 to 4 and R0; then ('timing_run'), a trial on networks 1 to 50 that projects
%              the time for 500; then ('full'), all 500, refused unless the checks passed (each a
%              few minutes; the stored 'full' took 2.8 min on 4 workers). ('handles') returns {K, helpers},
%              the helpers the response_dial_* files share.
% Seed     : K.seed = 20261028, Mersenne twister, the caller's stream restored after each use:
%            mutations re-drawn from each network's stored seed; signs of the small nudge to x*
%            (the return test, see check_rest) K.seed + 13000 + f for network f; check 3 K.seed +
%            12004; check R0 K.seed + 13001; encounters under rng(20261022) inside
%            mutate_with_refusal and encounter_rule_raw

    if nargin < 2, nWorkers = 3; end
    [K, H, NP] = conditional_cooperation_measures('np_handles');
    K = consts(K);
    switch stage
        case 'checks',  logged(K, 'checks', @() do_checks(K, H, NP, nWorkers));
        case {'timing_run', 'pilot'},   logged(K, 'timing_run', @() do_timing_run(K, H, nWorkers));
        case 'full',    logged(K, 'full', @() do_full(K, H, nWorkers));
        case 'handles', varargout = {K, helpers()};
        otherwise, error('unknown stage %s', stage);
    end
end


function K = consts(K)
% Settings S (i1 = the index of s = 1; dev = s = 0, 0.5, 1, the developmental rules' settings);
% 30 mutants per network; chunks of 25; the trial run's 50 networks; the encounters' seed; the 1%
% tolerance of the return and development tests and their 0.01 own-copy cooperation tolerance.
% rest: a run is at rest when no gene changes by tol in a step (cap
% steps at most); the return test's nudge pert x max(x*) per gene; the slope's input e, its stop
% slopeTol and cap slopeCap. K.version is stored in every chunk and compared on load.
    K.rd = fullfile(K.out, 'response_dial'); if ~isfolder(K.rd), mkdir(K.rd); end
    K.S = [0 0.5 1 1.5 2]; K.i1 = 3; K.dev = 1:3; K.nMut = 30; K.chunkN = 25; K.nPilot = 50;
    K.encSeed = 20261022; K.tol = 0.01; K.engTol = 0.01;
    K.rest = struct('tol', 1e-10, 'cap', 20000, 'pert', 1e-4, 'pertSeed', K.seed + 13000, ...
        'e', 1e-6, 'slopeTol', 1e-14, 'slopeCap', 1e6);
    K.version = 'RD_VERSION 2 (amendment 1, 2026-10-28)';
end


function H = helpers()
% Handles to the local functions the analysis files call.
    H = struct('check_dev', @check_dev, 'check_rest', @check_rest, 'load_set', @load_set, ...
        'ternary', @ternary, 'sensor_rows', @sensor_rows, 'dialed_slots', @dialed_slots, 'cd_slots', @cd_slots, ...
        'census_params', @census_params, 'rest_pair', @rest_pair, 'dial', @dial, 'apply_mut', @apply_mut, ...
        'readout', @readout, 'open_workers', @open_workers, 'rest_data', @rest_data, 'exclusion', @exclusion);
end


function logged(K, tag, fn)
% Runs fn with the console copied to log_<tag>.txt.
    diary(fullfile(K.rd, ['log_' tag '.txt'])); c = onCleanup(@() diary('off'));
    fprintf('\n===== RESPONSE DIAL %s %s UTC =====\n', tag, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('version key (stored in every chunk): %s\n', K.version);
    fn();
end


function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end


function [cp, sp] = census_params()
% The membership test's settings (cp, for membership_test) and simulate_pair's (sp): contrast sensing,
% Hill coefficient 1, 600 development steps, simulate_pair's default stopping rule; the constant
% inputs externalSignal1/2 are set per call (zero for the default defector).
    cp = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'backend', 'simulate_pair', 'sigma_mut', 0.2);
    sp = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), 'externalSignal2', zeros(4, 1));
end


function c = encounter_rule_config(K, H)
% The saturating-contrast run config, checked for what mutate_with_refusal('coop') assumes.
    cfg = H.load_configs(K); c = cfg{4};
    assert(strcmp(stored_value('sensingMode', c.sensingMode), 'contrast') && c.hillCoefficient == 1, 'config 4 is not saturating contrast');
    assert(c.numOpponents == 1 && strcmp(c.payoffMode, 'focal') && ~c.requireStability, 'coop assumptions');
    assert(~isfield(c, 'geneInputBias'), 'stored config already carries geneInputBias');
end


function close_workers()
    p = gcp('nocreate'); if ~isempty(p), delete(p); end
end


% =========================================================================
% INPUTS and membership classes
% =========================================================================
function I = build_inputs(K, NP)
% The 500 networks with their stored mutants (pool.mat), their own-copy cooperation and reaction R
% to the 80 partners (natural-partner chunks), the parents' M and B, the re-drawn mutations, the
% reaction fifths and the return test's nudge signs.
    S = load(fullfile(K.np, 'pool.mat'), 'Pt', 'F'); F = S.F; n = numel(F.ind);
    Z = NP.load_chunks(K, S, (1:n)', fullfile(K.np, 'chunks')); A = NP.measures(K, Z, F);
    [Mf, Bf] = focal_parents(K, F);
    D = redraw(F);
    I = struct('W', F.W, 'Wmut', F.Wmut, 'slotS', F.slot, 'dWS', F.dW, 'togS', F.isTog, 'mutSeed', F.mutSeed, ...
        'rep', F.rep, 'gen', F.gen, 'M', Mf, 'B', Bf, 'q', D.q, 'tog', D.tog, 'd', D.d, ...
        'self', Z.self, 'selfM', Z.selfM, 'R', A.R, 'Rdef', A.Rdef);
    I.fifth = fifths(A.R, A.Rdef);
    I.pert = perturbations(K, n);
end


function P = perturbations(K, n)
% The return test's nudge signs (+1 or -1 per gene), one vector per network f under rng(K.seed + 13000 +
% f); the caller's stream restored.
    P = zeros(n, 4); s0 = rng;
    for f = 1:n, rng(K.rest.pertSeed + f, 'twister'); P(f, :) = 2 * (rand(1, 4) < 0.5) - 1; end
    rng(s0);
end


function [Mf, Bf] = focal_parents(K, F)
% Each network's M and B: those of the first individual of its replicate and generation whose
% M .* B equals the stored W (asserted).
    d = load(fullfile(K.pub, 'saturating_contrast_100reps.mat'), 'results'); R = d.results; clear d;
    n = numel(F.ind); Mf = zeros(4, 4, n); Bf = zeros(4, 4, n);
    for f = 1:n
        rep = R{F.rep(f)}; k = find(double(rep.historyGenerations) == F.gen(f), 1);
        M = double(rep.historyM(:, :, :, k)); B = double(rep.historyB(:, :, :, k));
        Wall = reshape(M .* B, 16, size(M, 3))'; i = find(ismember(Wall, F.W(f, :), 'rows'), 1);
        assert(~isempty(i), 'network %d: W not in its snapshot', f);
        Mf(:, :, f) = M(:, :, i); Bf(:, :, f) = B(:, :, i);
        assert(isequal(reshape(Mf(:, :, f) .* Bf(:, :, f), 1, 16), F.W(f, :)), 'network %d: M.*B differs from W', f);
    end
    fprintf('network M and B (from the first carrier): %d networks, M.*B equals the stored W for all (asserted)\n', n);
end


function D = redraw(F)
% The 30 mutations of each network re-drawn from its stored seed F.mutSeed with the draw lines of
% the natural-partner mutations: slot (row, column), a toggle with probability 0.2, else a weight
% change 0.2 N(0, 1). The draws do not depend on M or B. The drawn change to B is kept, since the
% stored dW is the change in W (0 for a weight change of an absent edge). Check 1 compares them.
    n = numel(F.ind); D = struct('q', zeros(n, 30), 'tog', false(n, 30), 'd', zeros(n, 30));
    s0 = rng;
    for f = 1:n
        rng(F.mutSeed(f), 'twister');
        for t = 1:30
            g1 = randi(4); g2 = randi(4); tog = rand() < 0.2; dd = 0;
            if ~tog, dd = randn() * 0.2; end
            D.q(f, t) = sub2ind([4 4], g1, g2); D.tog(f, t) = tog; D.d(f, t) = dd;
        end
    end
    rng(s0);
end


function fifth = fifths(R, def)
% Fifths of the reaction R: edges at the 0th, 20th, ..., 100th percentiles of the defined R (ties
% joined); fifth 1 reacts least; a network with R undefined gets 0.
    ed = unique(prctile(R(def), 0:20:100)); ed(end) = Inf; fifth = zeros(size(R));
    for b = 1:numel(ed) - 1, fifth(def & R >= ed(b) & R < ed(b + 1)) = b; end
end


function D = rest_data(K, I, dirC)
% The from-rest set of the analyses (kept at s = 1 and the return test, check_rest, passed at every
% setting; asserted 338) and its subset in reaction fifths 2 to 5 (field prim; asserted 267); lost = the
% readout from rest is not positive.
    Z = load_set(K, I, (1:500)', dirC); G = check_rest(K, Z);
    kp = G.kept; nK = sum(kp); use = 1:5;
    if sum(G.pass(kp, 5)) < nK / 2 && sum(G.pass(kp, 4)) >= nK / 2, use = 1:4; end
    set = kp & all(G.pass(:, use), 2); prim = set & I.fifth >= 2;
    assert(sum(set) == 338 && sum(prim) == 267 && isequal(use, 1:5), 'set differs from the analysis');
    lost = ~(Z.rM > 0);
    D = struct('Z', Z, 'I', I, 'set', set, 'prim', prim, 'lost', lost, 'LR', squeeze(mean(lost, 2)), ...
        'rep', I.rep, 'dir', dirC);
end


function E = exclusion(D)
% Left out at every setting: toggles of the four dialed connections, whose size scales with s.
% pres: the mutated edge is present in M (toggle off removes it; toggle on adds an absent edge).
    Z = D.Z; n = numel(Z.r1); pres = false(size(Z.q));
    for i = 1:n, Mi = Z.M(:, :, i); pres(i, :) = Mi(Z.q(i, :)) ~= 0; end
    ex = ismember(Z.q, dialed_slots()) & Z.tog;
    E = struct('ex', ex, 'keep', ~ex, 'off', ex & pres, 'on', ex & ~pres, 'pres', pres);
end


function [Mm, Bm] = apply_mut(M, Bs, q, tog, d)
% A mutation of (M, B_s): a toggle flips M at the slot (a restored edge returns with its B_s
% value); a weight change adds the drawn change to B_s at the slot.
    Mm = M; Bm = Bs;
    if tog, Mm(q) = ~Mm(q); else, Bm(q) = Bm(q) + d; end
end


function cls = classes(I)
% The membership test (membership_test) on each network: 1 member (is_css); 2 cooperates with its clone
% but does not retaliate (is_weak, not is_css; printed as 'cooperate_only'); 3 other.
    cp = census_params(); n = size(I.W, 1); isC = false(n, 1); isW = false(n, 1); W = I.W;
    parfor i = 1:n
        r = membership_test(reshape(W(i, :), 4, 4), cp); isC(i) = r.is_css; isW(i) = r.is_weak;
    end
    cls = 3 * ones(n, 1); cls(isW) = 2; cls(isC) = 1;
end


% =========================================================================
% DYNAMICS: the resting state x*, runs to rest, the response slope
% =========================================================================
function [x, it, ok] = refine(W, bias, x, sp)
% x*: the clone pair's Euler step (the partner equal to x itself) from x until no gene changes by
% 1e-10 in a step, at most 1e6 steps; ok is false if it does not converge.
    km = sp.michaelisConstant; mu = sp.degradationRate; dt = sp.timeStep; ok = false;
    for it = 1:1e6
        u = max(W * x + partner_signal(x, sp) + bias, 0);
        xn = max(x + dt * (u ./ (km + u) - mu * x), 0);
        ch = max(abs(xn - x)); x = xn;
        if ch < 1e-10, ok = true; return; end
    end
end


function s = partner_signal(x, sp)
% simulate_pair's sensing (Hill, contrast) with the partner equal to x.
    n = sp.hillCoefficient; cn = x(1)^n; dn = x(2)^n; den = cn + dn + 1e-12;
    g = sp.sensorIntegrationRate * (cn / den - dn / den); s = [0; 0; g; -g];
end


function [sl, ok, it] = slope_fd(W, b, xstar, sp, R)
% Response slope: the partner signal held at partner_signal(x*); inputs [0; 0; +e; -e] and
% [0; 0; -e; +e] on the sensors, the two runs stepped together from x* until neither changes by
% slopeTol in a step (one stop for both, so a network blind to its sensors gives exactly 0), at most
% slopeCap steps; slope = (R(+e) - R(-e)) / (2 e), NaN if not converged.
    km = sp.michaelisConstant; mu = sp.degradationRate; dt = sp.timeStep; e = R.e;
    s0 = partner_signal(xstar, sp) + b; X = [xstar xstar]; E = [0 0; 0 0; e -e; -e e]; ok = false;
    for it = 1:R.slopeCap
        U = max(W * X + s0 + E, 0); Xn = max(X + dt * (U ./ (km + U) - mu * X), 0);
        ch = max(abs(Xn(:) - X(:))); X = Xn;
        if ch < R.slopeTol, ok = true; break; end
    end
    sl = ((X(1, 1) - X(2, 1)) - (X(1, 2) - X(2, 2))) / (2 * e);
    if ~ok, sl = NaN; end
end


function Wd = defector()
% The default defector: its only connection is the defection gene's activation of itself.
    Wd = zeros(4); Wd(2, 2) = 1;
end


% =========================================================================
% ONE NETWORK
% =========================================================================
function R = dial_focal(M, B, q, tog, d, pert, S, cfgE, sp, Rc)
% Every quantity for one network and its 30 mutants at each setting in S: the census snapshot x0
% of the clone pair and its refinement x*; alpha (alpha_and_beta at x*, readout gradient
% [1; -1; 0; 0]) and the response slope; the return test's run from the nudged x*; the parent from x*
% (rr0) and against the default defector (rDr; the defector at its own rest with its clone, from
% 0.005); each mutant from x* with b_s on both copies, run to rest (rM: lost if not positive).
% At s <= 1 also the developmental quantities and each mutant's census readout (cMut) and own-copy
% cooperation over 100 encounters (eMut).
    t0 = tic; mutate_with_refusal('reset'); W = M .* B; nS = numel(S); nM = numel(q); z4 = zeros(4, 1);
    f0 = readout(W, W, z4, z4, sp); x0 = f0(1:4); [xstar, it, ok] = refine(W, z4, x0, sp);
    R = struct('x0', x0, 'r1', x0(1) - x0(2), 'xstar', xstar, 'refIt', it, 'refOk', ok, ...
        'refMove', max(abs(xstar - x0)) / max(xstar), 'xs', nan(4, nS), 'rs', nan(1, nS), 'xr', nan(4, nS), ...
        'xrOk', false(1, nS), 'rD', nan(1, nS), 'alpha', nan(1, nS), 'bias', nan(2, nS), ...
        'eSelf', nan(1, nS), 'eMut', nan(nM, nS), 'cMut', nan(nM, nS), ...
        'xg', nan(4, nS), 'gOk', false(1, nS), 'rr0', nan(1, nS), 'rr0Ok', false(1, nS), ...
        'rM', nan(nM, nS), 'rMOk', false(nM, nS), 'rDr', nan(1, nS), 'rDrOk', false(1, nS), ...
        'slope', nan(1, nS), 'slopeOk', false(1, nS), 'slopeLin', nan(1, nS), 'sec', NaN);
    Wd = defector(); [xd, ~] = rest_pair(Wd, z4, 0.005 * ones(4, 1), Wd, z4, 0.005 * ones(4, 1), sp, Rc.tol, Rc.cap);
    xp = max(xstar + Rc.pert * max(xstar) * pert(:), 0);
    for k = 1:nS
        [Ws, bias, Bs] = dial(M, B, W, xstar, S(k)); R.bias(:, k) = bias(1:2);
        th = alpha_and_beta(Ws, xstar, partner_signal(xstar, sp) + bias, sp.degradationRate, sp.michaelisConstant, [1; -1; 0; 0]);
        R.alpha(k) = th.alpha; R.slopeLin(k) = th.v' * [0; 0; 1; -1];
        [R.slope(k), R.slopeOk(k)] = slope_fd(Ws, bias, xstar, sp, Rc);
        [R.xg(:, k), ~, ~, R.gOk(k)] = rest_pair(Ws, bias, xp, Ws, bias, xp, sp, Rc.tol, Rc.cap);
        [x1, ~, ~, R.rr0Ok(k)] = rest_pair(Ws, bias, xstar, Ws, bias, xstar, sp, Rc.tol, Rc.cap); R.rr0(k) = x1(1) - x1(2);
        [x1, ~, ~, R.rDrOk(k)] = rest_pair(Ws, bias, xstar, Wd, z4, xd, sp, Rc.tol, Rc.cap); R.rDr(k) = x1(1) - x1(2);
        dev = S(k) <= 1;
        if dev, R = dev_setting(R, k, Ws, bias, sp, cfgE, Wd); end
        cfg = cfgE; cfg.geneInputBias = bias;
        for t = 1:nM
            [Mm, Bm] = apply_mut(M, Bs, q(t), tog(t), d(t)); Wm = Mm .* Bm;
            [x1, ~, ~, R.rMOk(t, k)] = rest_pair(Wm, bias, xstar, Wm, bias, xstar, sp, Rc.tol, Rc.cap); R.rM(t, k) = x1(1) - x1(2);
            if dev
                fm = readout(Wm, Wm, bias, bias, sp); R.cMut(t, k) = fm(1) - fm(2);
                R.eMut(t, k) = mutate_with_refusal('coop', Wm, Wm, cfg);
            end
        end
    end
    mutate_with_refusal('reset'); R.sec = toc(t0);
end


function R = dev_setting(R, k, Ws, bias, sp, cfgE, Wd)
% Developmental quantities at one setting: the census snapshot of the dialed clone pair (xs, rs),
% its refinement (xr), its readout against the default defector (rD), and own-copy cooperation over
% 100 encounters under the runs' encounter rule with the constant input (eSelf).
    f = readout(Ws, Ws, bias, bias, sp); R.xs(:, k) = f(1:4); R.rs(k) = f(1) - f(2);
    [R.xr(:, k), ~, R.xrOk(k)] = refine(Ws, bias, f(1:4), sp);
    fd = readout(Ws, Wd, bias, zeros(4, 1), sp); R.rD(k) = fd(1) - fd(2);
    cfg = cfgE; cfg.geneInputBias = bias; R.eSelf(k) = mutate_with_refusal('coop', Ws, Ws, cfg);
end


% =========================================================================
% CHECKS 1-4 and R0
% =========================================================================
function do_checks(K, H, NP, nW)
% Builds and saves the inputs and membership classes, runs checks 1 to 4 and R0, and saves their
% pass flags; 'full' and the analysis read them.
    I = build_inputs(K, NP); cfgE = encounter_rule_config(K, H); [~, sp] = census_params();
    save(fullfile(K.rd, 'response_dial_inputs.mat'), 'I', '-v7.3');
    ok1 = check1(I);
    nW = open_workers(nW); c = onCleanup(@close_workers);
    cls = classes(I); save(fullfile(K.rd, 'response_dial_inputs.mat'), 'cls', '-append');
    fprintf('membership classes (membership_test): member %d, cooperate_only %d, other %d\n', ...
        sum(cls == 1), sum(cls == 2), sum(cls == 3));
    t = tic; ok2 = check2(K, I, cls, cfgE, sp); fprintf('  check 2 time %.1f min on %d workers\n', toc(t) / 60, nW);
    ok3 = check3(K, I, cfgE);
    t = tic; ok4 = check4(K, I, sp); fprintf('  check 4 time %.1f min\n', toc(t) / 60);
    okR = checkR0(K, I, sp);
    C = struct('check1', ok1, 'check2', ok2, 'check3', ok3, 'check4', ok4, 'checkR0', okR, 'version', K.version, ...
        'utc', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    save(fullfile(K.rd, 'checks_a2.mat'), '-struct', 'C');
    fprintf('\nCHECKS 1-4, R0: %s %s %s %s %s\n', pf(ok1), pf(ok2), pf(ok3), pf(ok4), pf(okR));
end


function ok = check1(I)
% Check 1: the re-drawn mutations applied to (M, B) give the stored mutant W exactly, and the
% re-drawn slot, toggle and change in W equal the stored ones, for every mutant.
    n = size(I.W, 1); nOk = 0; okD = 0;
    for f = 1:n
        M = I.M(:, :, f); B = I.B(:, :, f); W = M .* B;
        for t = 1:30
            [Mm, Bm] = apply_mut(M, B, I.q(f, t), I.tog(f, t), I.d(f, t)); Wm = Mm .* Bm;
            nOk = nOk + isequal(reshape(Wm, 1, 16), reshape(I.Wmut(f, t, :), 1, 16));
            q = I.q(f, t);
            okD = okD + (q == I.slotS(f, t) && I.tog(f, t) == I.togS(f, t) && Wm(q) - W(q) == I.dWS(f, t));
        end
    end
    ok = nOk == 30 * n && okD == 30 * n;
    fprintf('\nCHECK 1 (mutation replay): replayed W equals the stored mutant W exactly in %d of %d; re-drawn slot, toggle and dW equal the stored ones in %d of %d -> %s\n', ...
        nOk, 30 * n, okD, 30 * n, pf(ok));
end


function ok = check2(K, I, cls, cfgE, sp)
% Check 2, the dial code at s = 1: own-copy cooperation of each network and mutant equals the
% stored values exactly; census readouts equal the undialed simulate_pair path; and the numbers and
% mean census loss shares of members and of networks that cooperate but do not retaliate (census
% readout > 0) equal the fixed values in the code.
    n = size(I.W, 1); eS = nan(n, 1); eM = nan(n, 30); cM = nan(n, 30); r1 = nan(n, 1); pl = false(n, 1);
    M = I.M; B = I.B; q = I.q; tg = I.tog; d = I.d; Wm = I.Wmut; Wf = I.W; P = I.pert; Rc = K.rest;
    parfor i = 1:n
        R = dial_focal(M(:, :, i), B(:, :, i), q(i, :), tg(i, :), d(i, :), P(i, :), 1, cfgE, sp, Rc);
        eS(i) = R.eSelf; eM(i, :) = R.eMut'; cM(i, :) = R.cMut'; r1(i) = R.r1;
        pl(i) = plain_equal(Wf(i, :), squeeze(Wm(i, :, :)), R, sp);
    end
    okE = isequal(eS, I.self) && isequal(eM, I.selfM);
    fprintf('\nCHECK 2 (control at s = 1 through the dial code)\n');
    fprintf('  own-copy cooperation under the runs'' encounter rule equal to the stored values: networks %d of %d, mutants %d of %d; max |diff| %.3g\n', ...
        sum(eS == I.self), n, sum(eM(:) == I.selfM(:)), numel(eM), max(abs([eS - I.self; eM(:) - I.selfM(:)])));
    fprintf('  census readouts equal to the plain undialed path (zero externalSignal) for the network and all 30 mutants: %d of %d networks\n', sum(pl), n);
    room = r1 > 0; Lc = mean(~(cM > 0), 2); im = room & cls == 1; in = room & cls == 2;
    lm = mean(Lc(im)); ln = mean(Lc(in));
    okC = sum(im) == 210 && sum(in) == 214 && round(lm, 4) == 0.1151 && round(ln, 4) == 0.0871;
    fprintf('  census rule, networks with room > 0: member %d (mean loss share %.4f, %.6f), cooperate_only %d (%.4f, %.6f); expected: member 210 0.1151, cooperate_only 214 0.0871\n', ...
        sum(im), lm, lm, sum(in), ln, ln);
    ok = okE && all(pl) && okC;
    fprintf('CHECK 2: %s (encounter rule exact %s; census plain-path exact %s; membership shares %s)\n', pf(ok), pf(okE), pf(all(pl)), pf(okC));
end


function ok = plain_equal(Wf, Wm, R, sp)
% True if the undialed simulate_pair readouts of the network and its 30 mutants equal dial_focal's.
    E = ones(4); W = reshape(Wf, 4, 4);
    [~, fs] = simulate_pair(E, W, [], E, W, [], sp); ok = (fs(1) - fs(2)) == R.r1;
    for t = 1:30
        Wt = reshape(Wm(t, :), 4, 4); [~, fm] = simulate_pair(E, Wt, [], E, Wt, [], sp);
        ok = ok && (fm(1) - fm(2)) == R.cMut(t);
    end
end


function ok = check3(K, I, cfgE)
% Check 3: on 200 networks, the encounter output with geneInputBias absent equals that with
% zeros(4, 1), with the own copy and against the default defector; and a bias [0.5; 0; 0; 0]
% changes it for at least one network (so the check can fail).
    rng(K.seed + 12004, 'twister'); pick = randperm(size(I.W, 1), 200); Wd = defector();
    same = false(200, 2); moved = false(200, 1); cz = cfgE; cz.geneInputBias = zeros(4, 1);
    cb = cfgE; cb.geneInputBias = [0.5; 0; 0; 0];
    for j = 1:200
        W = reshape(I.W(pick(j), :), 4, 4);
        same(j, 1) = isequaln(encounter_rule_raw(W, W, cfgE), encounter_rule_raw(W, W, cz));
        same(j, 2) = isequaln(encounter_rule_raw(W, Wd, cfgE), encounter_rule_raw(W, Wd, cz));
        moved(j) = ~isequaln(encounter_rule_raw(W, W, cfgE), encounter_rule_raw(W, W, cb));
    end
    ok = all(same(:)) && any(moved);
    fprintf('\nCHECK 3 (geneInputBias of the runs inert when unused; 200 networks, seed %d)\n', K.seed + 12004);
    fprintf('  absent vs zeros identical (fitness, rate, per-individual cooperation): own copy %d of 200, against the built-in defector %d of 200\n', ...
        sum(same(:, 1)), sum(same(:, 2)));
    fprintf('  can it fail: bias [0.5; 0; 0; 0] changes the own-copy output in %d of 200\n', sum(moved));
    fprintf('CHECK 3: %s\n', pf(ok));
end


function R = encounter_rule_raw(W, P, config)
% The population mutate_with_refusal('coop') builds (100 copies of W with each other, or 100 of W
% against 100 of P), evaluated by evolve_population under rng(20261022): its full output, not
% cached; the caller's stream restored.
    NENC = 100; s = rng; rng(20261022, 'twister');
    if isequal(W, P)
        n = NENC; Wp = repmat(W, 1, 1, n); opp = [n/2 + (1:n/2), 1:n/2]';
    else
        n = 2 * NENC; Wp = cat(3, repmat(W, 1, 1, NENC), repmat(P, 1, 1, NENC)); opp = [NENC + (1:NENC), 1:NENC]';
    end
    cfg = config; cfg.precomputedOpponents = opp; cfg.populationSize = n;
    Kk = ones(config.numGenes, n, 'single') * config.michaelisConstant;
    R = evolve_population('__test_evaluateFitness__', true(config.numGenes, config.numGenes, n), single(Wp), Kk, cfg);
    rng(s);
end


function ok = check4(K, I, sp)
% Check 4, at s = 0 (the sensors cut off from C and D): with no stopping rule (all 600 steps) the
% readout against the default defector equals that with the clone (G = 0) and mutations of the
% sensors' inputs (rows 3 and 4) leave the census readout unchanged, to 1e-9; from rest, each such
% mutant's readout equals the parent run from x* for the mutant's own step count, to 1e-12.
% Printed beside, not part of the check: the shift against the parent run to its own rest.
    n = size(I.W, 1); g = nan(n, 1); dm = nan(n, 1); de = nan(n, 1); dr = nan(n, 1); nS = zeros(n, 1);
    nNo = zeros(n, 1); kp = false(n, 1);
    M = I.M; B = I.B; q = I.q; tg = I.tog; d = I.d; P = I.pert; sp0 = sp; sp0.stabilityThreshold = 0;
    Wd = defector(); z4 = zeros(4, 1); Rc = K.rest; tol = K.tol;
    parfor i = 1:n
        Mi = M(:, :, i); Bi = B(:, :, i); W = Mi .* Bi;
        f0 = readout(W, W, z4, z4, sp); [xstar, ~, okR] = refine(W, z4, f0(1:4), sp);
        kp(i) = kept_rest(Mi, Bi, W, xstar, okR, P(i, :), sp, Rc, tol);
        [Ws, bias, Bs] = dial(Mi, Bi, W, xstar, 0);
        f = readout(Ws, Ws, bias, bias, sp0); fd = readout(Ws, Wd, bias, z4, sp0);
        r0 = f(1) - f(2); g(i) = r0 - (fd(1) - fd(2));
        xp = rest_pair(Ws, bias, xstar, Ws, bias, xstar, sp, Rc.tol, Rc.cap); rp0 = xp(1) - xp(2);
        dmax = 0; demax = 0; drmax = 0; ns = 0; nno = 0;
        for t = find(ismember(q(i, :), sensor_rows()))
            [Mm, Bm] = apply_mut(Mi, Bs, q(i, t), tg(i, t), d(i, t)); Wm = Mm .* Bm;
            fm = readout(Wm, Wm, bias, bias, sp0); dmax = max(dmax, abs((fm(1) - fm(2)) - r0)); ns = ns + 1;
            [xm, ~, itm, okm] = rest_pair(Wm, bias, xstar, Wm, bias, xstar, sp, Rc.tol, Rc.cap);
            xq = rest_pair(Ws, bias, xstar, Ws, bias, xstar, sp, -1, itm);
            demax = max(demax, abs((xm(1) - xm(2)) - (xq(1) - xq(2))));
            drmax = max(drmax, abs((xm(1) - xm(2)) - rp0)); nno = nno + ~okm;
        end
        dm(i) = dmax; de(i) = demax; dr(i) = drmax; nS(i) = ns; nNo(i) = nno;
    end
    ok = all(abs(g) <= 1e-9) && all(dm <= 1e-9) && all(de <= 1e-12);
    fprintf('\nCHECK 4 (cut is exact at s = 0)\n');
    fprintf('  fixed-start part (stabilityThreshold 0, all 600 steps): |G(n, 0)| max %.3g over %d networks, > 1e-9 in %d\n', max(abs(g)), n, sum(abs(g) > 1e-9));
    fprintf('    sensor-input mutants (rows 3-4): %d; max |census readout shift| %.3g; > 1e-9 in %d networks\n', sum(nS), max(dm), sum(dm > 1e-9));
    fprintf('  from-rest part (parent run for the mutant''s step count): max |shift| %.3g; > 1e-12 in %d networks; mutant runs not at rest by the cap %d\n', ...
        max(de), sum(de > 1e-12), sum(nNo));
    fprintf('  beside, not gating (parent run to its own rest): max |shift| %.3g; > 1e-9 in %d networks (kept %d, not kept %d; kept = passes the return test at s = 1, %d of %d)\n', ...
        max(dr), sum(dr > 1e-9), sum(dr > 1e-9 & kp), sum(dr > 1e-9 & ~kp), sum(kp), n);
    fprintf('CHECK 4: %s\n', pf(ok));
end


function k = kept_rest(M, B, W, xstar, okR, pv, sp, Rc, tol)
% Kept at s = 1: x* converged, r* > 0, and the run from the perturbed x* rests within 1% of x* in
% state and readout.
    rs = xstar(1) - xstar(2); k = false;
    if ~okR || ~(rs > 0), return; end
    [Ws, b] = dial(M, B, W, xstar, 1); xp = max(xstar + Rc.pert * max(xstar) * pv(:), 0);
    x = rest_pair(Ws, b, xp, Ws, b, xp, sp, Rc.tol, Rc.cap);
    k = max(abs(x - xstar)) / max(xstar) <= tol && abs((x(1) - x(2)) - rs) <= tol * abs(rs);
end


function ok = checkR0(K, I, sp)
% Check R0: rest_pair with no stopping for 600 steps equals simulate_pair bit for bit (stability
% threshold 0, the same start and constant inputs), on 200 (network, setting) draws, with the own
% copy and against the default defector (from 0.005).
    rng(K.seed + 13001, 'twister'); pick = randi(size(I.W, 1), 200, 1); ks = randi(numel(K.S), 200, 1);
    sp0 = sp; sp0.stabilityThreshold = 0; z4 = zeros(4, 1); Wd = defector(); same = false(200, 2);
    for j = 1:200
        i = pick(j); M = I.M(:, :, i); B = I.B(:, :, i); W = M .* B;
        f0 = readout(W, W, z4, z4, sp); xstar = refine(W, z4, f0(1:4), sp);
        [Ws, b] = dial(M, B, W, xstar, K.S(ks(j))); xd = 0.005 * ones(4, 1);
        s1 = sp0; s1.initialState1 = xstar; s1.initialState2 = xstar; f = readout(Ws, Ws, b, b, s1);
        [a1, a2] = rest_pair(Ws, b, xstar, Ws, b, xstar, sp, -1, 600); same(j, 1) = isequal(f, [a1; a2]);
        s2 = sp0; s2.initialState1 = xstar; s2.initialState2 = xd; f = readout(Ws, Wd, b, z4, s2);
        [a1, a2] = rest_pair(Ws, b, xstar, Wd, z4, xd, sp, -1, 600); same(j, 2) = isequal(f, [a1; a2]);
    end
    ok = all(same(:));
    fprintf('\nCHECK R0 (instrument; seed %d): rest_pair equals simulate_pair bit for bit after 600 steps: own copy %d of 200, against the built-in defector %d of 200 -> %s\n', ...
        K.seed + 13001, sum(same(:, 1)), sum(same(:, 2)), pf(ok));
end


function qs = sensor_rows()
% Linear slots of W(3:4, :): inputs to the sensor genes.
    [r, c] = ndgrid(3:4, 1:4); qs = sub2ind([4 4], r(:), c(:))';
end


function qs = dialed_slots()
% Linear slots of W(1:2, 3:4): the four dialed connections (sensors to C and D).
    [r, c] = ndgrid(1:2, 3:4); qs = sub2ind([4 4], r(:), c(:))';
end


function qs = cd_slots()
% Linear slots of W(1:2, 1:2): the inputs among C and D.
    [r, c] = ndgrid(1:2, 1:2); qs = sub2ind([4 4], r(:), c(:))';
end


% =========================================================================
% SCORING in chunks of 25 networks
% =========================================================================
function [Z, wall, ws, nW] = score(K, H, fset, dirC, nW)
% Scores the networks fset chunk by chunk into dirC; a chunk already there is loaded after its
% stored inputs are checked against the current ones (a stale chunk stops the run).
    L = load(fullfile(K.rd, 'response_dial_inputs.mat'), 'I'); I = L.I; cfgE = encounter_rule_config(K, H); [~, sp] = census_params();
    if ~isfolder(dirC), mkdir(dirC); end
    nW = open_workers(nW); c = onCleanup(@close_workers);
    nC = ceil(numel(fset) / K.chunkN); Z = cell(nC, 1); t0 = tic;
    fprintf('settings %s; encounter seed %d; %d networks in %d chunks\n', mat2str(K.S), K.encSeed, numel(fset), nC);
    for k = 1:nC
        id = fset((k - 1) * K.chunkN + 1:min(k * K.chunkN, numel(fset)));
        f = fullfile(dirC, sprintf('chunk_%03d.mat', k)); want = chunk_inputs(K, I, id, sp);
        if isfile(f)
            Zk = load(f); assert_chunk(Zk, want, f); fprintf('  chunk %d/%d loaded (resumed)\n', k, nC);
        else
            Zk = run_chunk(want, cfgE, sp); save(f, '-struct', 'Zk');
        end
        Z{k} = Zk; fprintf('  chunk %d/%d (%d networks) done, elapsed %.1f min\n', k, nC, numel(id), toc(t0) / 60);
    end
    wall = toc(t0); Z = cat_chunks(Z); ws = sum(Z.sec);
end


function in = chunk_inputs(K, I, id, sp)
% Everything a chunk depends on, stored in it: networks, mutations, perturbations, settings, the
% encounters' seed, simulate_pair's settings, the rest constants and K.version.
    in = struct('focal', id(:), 'W', I.W(id, :), 'M', I.M(:, :, id), 'B', I.B(:, :, id), 'q', I.q(id, :), ...
        'tog', I.tog(id, :), 'd', I.d(id, :), 'pert', I.pert(id, :), 'mutSeed', I.mutSeed(id), 'S', K.S, ...
        'encSeed', K.encSeed, 'sp', sp, 'rest', K.rest, 'version', K.version);
end


function assert_chunk(Z, want, f)
    % a stored chunk holds the older spelling of sp.sensingMode (stored_value)
    if isfield(Z, 'sp') && isfield(Z.sp, 'sensingMode'), Z.sp.sensingMode = stored_value('sensingMode', Z.sp.sensingMode); end
    for fn = fieldnames(want)'
        assert(isfield(Z, fn{1}) && isequal(Z.(fn{1}), want.(fn{1})), 'chunk %s: %s does not match (stale)', f, fn{1});
    end
end


function Z = run_chunk(in, cfgE, sp)
    n = numel(in.focal); C = cell(n, 1); M = in.M; B = in.B; q = in.q; tg = in.tog; d = in.d; P = in.pert;
    S = in.S; Rc = in.rest;
    parfor i = 1:n
        C{i} = dial_focal(M(:, :, i), B(:, :, i), q(i, :), tg(i, :), d(i, :), P(i, :), S, cfgE, sp, Rc);
    end
    Z = in; R = [C{:}]; nS = numel(S);
    for fn = {'x0', 'xstar'}, Z.(fn{1}) = [R.(fn{1})]'; end
    for fn = {'r1', 'refIt', 'refOk', 'refMove', 'sec'}, Z.(fn{1}) = [R.(fn{1})]'; end
    for fn = {'xs', 'xr', 'xg', 'bias'}, Z.(fn{1}) = permute(cat(3, R.(fn{1})), [3 1 2]); end     % n x 4 (2) x nS
    for fn = {'rs', 'xrOk', 'rD', 'alpha', 'eSelf', 'gOk', 'rr0', 'rr0Ok', 'rDr', 'rDrOk', 'slope', 'slopeOk', 'slopeLin'}
        Z.(fn{1}) = reshape([R.(fn{1})], nS, n)';
    end
    for fn = {'eMut', 'cMut', 'rM', 'rMOk'}, Z.(fn{1}) = permute(cat(3, R.(fn{1})), [3 1 2]); end   % n x 30 x nS
end


function Z = cat_chunks(C)
    Z = C{1}; keep = {'S', 'encSeed', 'sp', 'rest', 'version'};
    for k = 2:numel(C)
        for f = fieldnames(Z)'
            if ismember(f{1}, keep), continue; end
            if ismember(f{1}, {'M', 'B'}), Z.(f{1}) = cat(3, Z.(f{1}), C{k}.(f{1}));
            else, Z.(f{1}) = cat(1, Z.(f{1}), C{k}.(f{1})); end
        end
    end
end


function Z = load_set(K, I, fset, dirC)
    [~, sp] = census_params(); nC = ceil(numel(fset) / K.chunkN); C = cell(nC, 1);
    for k = 1:nC
        id = fset((k - 1) * K.chunkN + 1:min(k * K.chunkN, numel(fset)));
        f = fullfile(dirC, sprintf('chunk_%03d.mat', k)); assert(isfile(f), 'missing %s', f);
        C{k} = load(f); assert_chunk(C{k}, chunk_inputs(K, I, id, sp), f);
    end
    Z = cat_chunks(C);
end


% =========================================================================
% THE DEVELOPMENT TEST AND THE RETURN TEST
% =========================================================================
function G = check_dev(K, Z, self1)
% The development test at s = 0, 0.5, 1 (does the dialed clone pair still develop to x*?). Census
% part: the dialed clone pair's state after development from the fixed start lies within 1% of x*
% (relative to max x*) and its readout within 1% of the undialed census readout r1. Encounter-rule
% part: own-copy cooperation under the runs' encounter rule within 0.01 of its value at s = 1.
% Kept at s = 1: own-copy cooperation >= 0.99 (keptE), r1 > 0 (keptC). refMatch (computed, not
% printed here): the refined dialed state meets the census part's two 1% tolerances against x* and r*.
    n = numel(Z.r1); nS = numel(K.dev); xs = Z.xstar; mx = max(xs, [], 2); rstar = xs(:, 1) - xs(:, 2);
    G.state = false(n, nS); G.read = false(n, nS); G.enc_rule = false(n, nS); G.refMatch = false(n, nS);
    for k = K.dev
        G.state(:, k) = max(abs(Z.xs(:, :, k) - xs), [], 2) ./ mx <= K.tol;
        G.read(:, k) = abs(Z.rs(:, k) - Z.r1) <= K.tol * abs(Z.r1);
        G.enc_rule(:, k) = abs(Z.eSelf(:, k) - Z.eSelf(:, K.i1)) <= K.engTol;
        xr = Z.xr(:, :, k); rr = xr(:, 1) - xr(:, 2);
        G.refMatch(:, k) = max(abs(xr - xs), [], 2) ./ mx <= K.tol & abs(rr - rstar) <= K.tol * abs(rstar);
    end
    G.census = G.state & G.read; G.pass = G.census & G.enc_rule;
    G.keptE = Z.eSelf(:, K.i1) >= 0.99; G.keptC = Z.r1 > 0;
    G.okSelf = isequal(Z.eSelf(:, K.i1), self1);
end


function G = check_rest(K, Z)
% The return test at each setting (does the clone pair return to x* after a small nudge?): the run
% from the nudged x* rests within 1% of x* in state (relative to max x*) and in readout (relative to
% r* = x*_C - x*_D). Kept: x* converged, r* > 0 and the test passed at s = 1.
    n = numel(Z.r1); nS = numel(K.S); xs = Z.xstar; mx = max(xs, [], 2); rstar = xs(:, 1) - xs(:, 2);
    G.state = false(n, nS); G.read = false(n, nS);
    for k = 1:nS
        xg = Z.xg(:, :, k); rg = xg(:, 1) - xg(:, 2);
        G.state(:, k) = max(abs(xg - xs), [], 2) ./ mx <= K.tol;
        G.read(:, k) = abs(rg - rstar) <= K.tol * abs(rstar);
    end
    G.pass = G.state & G.read; G.rstar = rstar;
    G.kept = Z.refOk & rstar > 0 & G.pass(:, K.i1);
end


% =========================================================================
% 'timing_run' AND 'full'
% =========================================================================
function do_timing_run(K, H, nW)
% The trial run: networks 1 to 50 into timing_run_a1/, with the return test, rest counts, slope table and
% development test printed; the projected time for 500 (the larger of the wall-time and worker-time
% scalings) is printed.
    L0 = load(fullfile(K.rd, 'response_dial_inputs.mat'), 'I'); I = L0.I; fset = (1:K.nPilot)';
    [Z, wall, ws, nW] = score(K, H, fset, fullfile(K.rd, 'timing_run_a1'), nW);
    timing_run_report(K, Z, check_rest(K, Z), check_dev(K, Z, I.self(fset)));
    hW = wall / K.nPilot * 500 / 3600; hS = ws / nW / K.nPilot * 500 / 3600; h = max(hW, hS);
    fprintf('\nTIMING: %d networks, wall %.1f min on %d workers, worker-seconds %.0f (%.1f s per network)\n', ...
        K.nPilot, wall / 60, nW, ws, ws / K.nPilot);
    fprintf('  projection for 500: %.2f h (wall-scaled), %.2f h (worker-seconds / %d); read %.2f h; on 4 workers %.2f h, on 2 %.2f h\n', ...
        hW, hS, nW, h, h * nW / 4, h * nW / 2);
    T = struct('wall', wall, 'ws', ws, 'nW', nW, 'h', h, 'version', K.version);
    save(fullfile(K.rd, 'timing_run_manifest.mat'), 'T');
    write_timing_run_losses(K, Z);
end


function timing_run_report(K, Z, G, D)
% s = 2 is dropped (top 1.5) if fewer than half the kept networks pass the return test there but at least
% half pass at 1.5.
    n = numel(Z.r1); S = K.S; kp = G.kept;
    fprintf('\nx* REFINEMENT: converged %d of %d; relative move median %.3g, max %.3g\n', sum(Z.refOk), n, median(Z.refMove), max(Z.refMove));
    fprintf('\nRETURN TEST. Kept at s = 1 (refined, r* > 0, passes at s = 1): %d of %d (refined %d, r* > 0 %d)\n', ...
        sum(kp), n, sum(Z.refOk), sum(G.rstar > 0));
    fprintf('  %-5s %-22s %-22s %-14s\n', 's', 'pass all 50 (st/rd)', 'pass kept (st/rd)', 'test run at rest');
    for k = 1:numel(S)
        fprintf('  %-5.1f %3d (%3d/%3d)          %3d (%3d/%3d)          %3d of %d\n', S(k), sum(G.pass(:, k)), sum(G.state(:, k)), ...
            sum(G.read(:, k)), sum(G.pass(kp, k)), sum(G.state(kp, k)), sum(G.read(kp, k)), sum(Z.gOk(:, k)), n);
    end
    p2 = sum(G.pass(kp, 5)); p15 = sum(G.pass(kp, 4)); nK = sum(kp);
    fprintf('  pass at all five: all %d, kept %d; fallback (kept: < half at 2 and >= half at 1.5): %d / %d at 2, %d / %d at 1.5 -> %s\n', ...
        sum(all(G.pass, 2)), sum(all(G.pass(kp, :), 2)), p2, nK, p15, nK, ternary(p2 < nK / 2 && p15 >= nK / 2, 'APPLIES (top 1.5)', 'does not apply'));
    if p2 < nK / 2 && p15 >= nK / 2, use = 1:4; else, use = 1:5; end
    fprintf('  set (kept, pass at every setting used %s): %d\n', mat2str(S(use)), sum(kp & all(G.pass(:, use), 2)));
    fprintf('\nRUN TO REST, not at rest by the cap (20,000 steps), per setting %s:\n', mat2str(S));
    fprintf('  parent from x*: %s; parent gate run: %s; parent vs defector: %s; mutants (of %d): %s\n', ...
        mat2str(sum(~Z.rr0Ok, 1)), mat2str(sum(~Z.gOk, 1)), mat2str(sum(~Z.rDrOk, 1)), n * 30, mat2str(squeeze(sum(sum(~Z.rMOk, 1), 2))'));
    fprintf('  mutants not at rest, kept networks only (of %d): %s\n', nK * 30, mat2str(squeeze(sum(sum(~Z.rMOk(kp, :, :), 1), 2))'));
    resp_table(K, Z, G, kp & all(G.pass(:, use), 2));
    dev_table(K, D);
end


function resp_table(K, Z, G, set)
    sr = Z.slope ./ Z.slope(:, K.i1); Gr = Z.rr0 - Z.rDr; gr = Gr ./ Gr(:, K.i1);
    fprintf('\nRESPONSE SLOPE (partner held at x*, e = 1e-6) and G FROM REST, medians (kept %d, set %d)\n', sum(G.kept), sum(set));
    fprintf('  slope runs not converged per setting: %s; max |slope(s = 0)| %.3g\n', mat2str(sum(~Z.slopeOk, 1)), max(abs(Z.slope(:, 1))));
    d = abs(Z.slope - Z.slopeLin) ./ abs(Z.slopeLin); d = d(isfinite(d));
    fprintf('  finite-difference vs linearized slope (alpha_and_beta v''[0;0;1;-1]): median rel diff %.3g, p90 %.3g\n', median(d), prctile(d, 90));
    fprintf('  %-5s %-22s %-22s %-14s %-22s %-12s\n', 's', 'slope/slope(1), kept', 'slope/slope(1), set', 'slope, kept', 'G_rest/G_rest(1), kept', 'G_rest, kept');
    for k = 1:numel(K.S)
        a = sr(G.kept, k); b = sr(set, k); g = gr(G.kept, k);
        fprintf('  %-5.1f %8.4f (%3d)         %8.4f (%3d)         %8.4f       %8.4f (%3d)          %8.4f\n', K.S(k), md(a), sum(isfinite(a)), ...
            md(b), sum(isfinite(b)), md(Z.slope(G.kept, k)), md(g), sum(isfinite(g)), md(Gr(G.kept, k)));
    end
end


function m = md(x)
    x = x(isfinite(x)); if isempty(x), m = NaN; else, m = median(x); end
end


function dev_table(K, D)
    fprintf('\nDEVELOPMENT TEST at 0 / 0.5 / 1: kept by the encounter rule %d, kept by the census %d\n', sum(D.keptE), sum(D.keptC));
    fprintf('  %-5s %-28s %-28s %-28s\n', 's', 'census part all/E/C', 'encounter part all/E/C', 'joint all/E/C');
    for k = K.dev
        fprintf('  %-5.1f %3d / %3d / %3d              %3d / %3d / %3d              %3d / %3d / %3d\n', K.S(k), ...
            sum(D.census(:, k)), sum(D.census(D.keptE, k)), sum(D.census(D.keptC, k)), sum(D.enc_rule(:, k)), ...
            sum(D.enc_rule(D.keptE, k)), sum(D.enc_rule(D.keptC, k)), sum(D.pass(:, k)), sum(D.pass(D.keptE, k)), sum(D.pass(D.keptC, k)));
    end
    fprintf('  joint pass at 0, 0.5 and 1: kept by the encounter rule %d, kept by the census %d; own-copy cooperation at s = 1 equals stored: %s\n', ...
        sum(all(D.pass(D.keptE, K.dev), 2)), sum(all(D.pass(D.keptC, K.dev), 2)), pf(D.okSelf));
end


function write_timing_run_losses(K, Z)
% The trial run's per-network loss shares under the three rules, written to timing_run_losses.txt only.
    LR = squeeze(mean(~(Z.rM > 0), 2)); LE = squeeze(mean(Z.eMut < 0.5, 2)); LC = squeeze(mean(~(Z.cMut > 0), 2));
    f = fullfile(K.rd, 'timing_run_losses.txt'); h = fopen(f, 'w'); c = onCleanup(@() fclose(h));
    fprintf(h, 'Timing run, raw values. network | from-rest L at s = %s | encounter-rule L at 0 0.5 1 | census L at 0 0.5 1\n', mat2str(K.S));
    for i = 1:numel(Z.focal)
        fprintf(h, '%d | %s| %s| %s\n', Z.focal(i), sprintf('%.4f ', LR(i, :)), sprintf('%.4f ', LE(i, K.dev)), sprintf('%.4f ', LC(i, K.dev)));
    end
    fprintf('timing run loss values written to %s (raw; not summarised here)\n', f);
end


function do_full(K, H, nW)
% All 500 networks into chunks/; refuses unless checks 1 to 4 and R0 passed, and unless the checks and
% the trial run ran under the same K.version.
    C = load(fullfile(K.rd, 'checks_a2.mat')); T = load(fullfile(K.rd, 'timing_run_manifest.mat'), 'T'); T = T.T;
    assert(C.check1 && C.check2 && C.check3 && C.check4 && C.checkR0, 'a check (1-4, R0) failed; not scoring');
    assert(strcmp(C.version, K.version) && strcmp(T.version, K.version), 'checks or timing run ran under another RD_VERSION');
    fprintf('timing run projection %.2f h on %d workers\n', T.h, T.nW);
    [Z, wall, ws, nW] = score(K, H, (1:500)', fullfile(K.rd, 'chunks'), nW);
    fprintf('  500 networks, wall %.1f min on %d workers, worker-seconds %.0f; non-finite rM %d\n', ...
        wall / 60, nW, ws, sum(~isfinite(Z.rM(:))));
end
