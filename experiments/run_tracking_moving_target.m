function out = run_tracking_moving_target(cfg)
% RUN_TRACKING_MOVING_TARGET  With the target moving, does lowering responsiveness alpha or ||x*||
% in an evolving population of design 1 change the load it carries? Lowers one factor to 0.70 at
% the split (generation 350) and evolves four arms for 150 generations.
%
% Produces : Table S14 (the rows with the target moving); Fig. S2; SI S5.1; Fig. 6 (read by
%            fig6_tracking_load and figS2_moving_target)
% Inputs   : publication/tracking_runs/tracking_runs_N1000.mat (the generation-350
%              populations, written by run_tracking_task)
%            publication/tracking_single_genotypes/lowering_settings.mat (written by run_tracking_single_genotypes;
%              the run checks its own (b, delta) against it)
% Outputs  : publication/tracking_moving_target/moving_target.mat and moving_target.txt
% Run      : run_tracking_moving_target()   (about 54 min on 2 workers); the tolerance (see
%              check_c1c) is that printed by a short run, run_tracking_moving_target(struct('pilot',
%              true)), and stored in moving_target_timing_run.mat, which the main run reads
% Seed     : replicate r continues with rng(20261001 + r, 'twister') in the control and both
%            lowered arms (paired by seed), 20261101 + r in the second control; sampled individuals
%            and the draws that score the per-mutation load use streams 20261201 + r, + 7919 r and
%            + 104729 r; bootstrap 20261002
%
% The split, the joint solve for (b, delta) on its best genotype, the re-centered optimum and the
% four arms are those of run_tracking_held_still, but no generations without mutation follow.
% Read over the last 60 generations: the load carried (one minus the population's mean fitness
% over its best genotype's, across the scoring window), mean fitness, the best genotype's
% tracking fitness, and the per-generation mutation cost (one minus the offspring's mean fitness
% over that of their unmutated parents).

    cfg = defaults(cfg);
    setup_process();
    D = getfield(load(cfg.inFile), 'out');
    P = D.P; R = D.moving; f = cfg.designIndex;
    E = environment(D, R, f, P, cfg);
    Rs = struct('popTh', R.popTh(f, cfg.repList), 'popFreq', R.popFreq(f, cfg.repList), ...
        'ThBest', R.ThBest(f, cfg.repList), 'maint', num2cell(R.maint(f, cfg.repList)));
    clear D R
    S1 = lowering_settings(cfg);

    diary_file(cfg.txtFile, true);
    fprintf('\n=========== TRACKING POPULATION TEST (%s) ===========\n', ...
        ternary(cfg.pilot, 'TIMING RUN', 'MAIN RUN'));
    fprintf('source %s | a = %.2f moving | N = %d | replicates %s | %d generations after the split (gen %d)\n', ...
        cfg.inFile, E.a, E.N, mat2str(cfg.repList), cfg.nGen, E.G0);
    fprintf('mutation distribution p = %.2f, sd = %.3f, s = %.1f | knob level m = %.2f | load window: last %d generations\n', ...
        P.mutation_prob, P.mutation_sd, P.selection_strength, cfg.m, cfg.loadWin);

    out = struct('cfg', cfg, 'P', P, 'E', rmfield(E, 'uPath'));
    start_pool(cfg.nWorkers);
    t0 = tic;
    nR = numel(cfg.repList); splits = cell(nR, 1);
    parfor q = 1:nR
        splits{q} = prepare_split(Rs(q), cfg.repList(q), E, P, cfg);
    end
    out.splits = splits; out.splitMin = toc(t0)/60;
    fprintf('split preparation: %.1f min for %d replicates\n', out.splitMin, nR);
    out.c1 = report_c1(splits, S1, E, cfg);
    out.verif = report_verification(splits, cfg);
    save(cfg.outFile, 'out');

    out.jobs = run_continuations(splits, E, P, cfg, out);
    out.contMin = toc(t0)/60 - out.splitMin;
    if cfg.pilot
        out.pilot = timing_run_tolerance(out, cfg);
    else
        out.analysis = analyse(out, Rs, cfg);
    end
    out.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    save(cfg.outFile, 'out');
    fprintf('\nSaved %s (%.1f min total)\n', cfg.outFile, toc(t0)/60);
    diary_file(cfg.txtFile, false);
end


% Configuration
function cfg = defaults(cfg)
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    pub = fullfile(grn_data_file('publication'));
    pilot = isfield(cfg, 'pilot') && cfg.pilot;
    d = struct('inFile', fullfile(pub, 'tracking_runs', 'tracking_runs_N1000.mat'), ...
        'stage1File', fullfile(pub, 'tracking_single_genotypes', 'lowering_settings.mat'), ...
        'outDir', fullfile(pub, 'tracking_moving_target'), 'pilot', false, ...
        'designIndex', 1, 'repList', 1:ternary(pilot, 3, 20), 'nGen', ternary(pilot, 30, 150), ...
        'loadWin', ternary(pilot, 20, 60), 'm', 0.7, 'contLevels', [0.9 0.8 0.7], ...
        'nWorkers', 2, 'seed', 20261001, 'seed2', 20261101, 'readSeed', 20261201, ...
        'nVer', 200, 'nDraws', 100, 'geoEvery', 10, 'betaEvery', 25, 'nTrack', 20, 'nBeta', 4, ...
        'trajEvery', 30, 'loadEvery', 5, ...
        'nBoot', 2000, 'bootSeed', 20261002, 'targetTol', 0.20, 'otherTol', 0.20, ...
        'flipDist', 0.5, 'waveSize', 8, 'tolerance', NaN);
    fn = fieldnames(d);
    for k = 1:numel(fn)
        if ~isfield(cfg, fn{k}), cfg.(fn{k}) = d.(fn{k}); end
    end
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    if ~cfg.pilot && isnan(cfg.tolerance), cfg.tolerance = stored_tolerance(cfg.outDir); end
    stem = ternary(cfg.pilot, 'moving_target_timing_run', 'moving_target');
    if ~isfield(cfg, 'outFile'), cfg.outFile = fullfile(cfg.outDir, [stem '.mat']); end
    if ~isfield(cfg, 'txtFile'), cfg.txtFile = fullfile(cfg.outDir, [stem '.txt']); end
    cfg.arms = {'CONTROL', 'CONTROL2', 'ALPHA-DOWN', 'XNORM-DOWN'};
end


function tol = stored_tolerance(outDir)
%STORED_TOLERANCE  the control-check tolerance the timing run printed and stored.
    f = fullfile(outDir, 'moving_target_timing_run.mat');
    if ~exist(f, 'file')
        error('run_tracking_moving_target:tolerance', ...
            'No tolerance given and %s is missing: run run_tracking_moving_target(struct(''pilot'', true)) first', f);
    end
    o = getfield(load(f), 'out');
    tol = o.pilot.tol.tolerance;
end


function E = environment(D, R, f, P, cfg)
%ENVIRONMENT  the moving-target run's circle, continued past the split (gen G0 = 350).
    E.a = D.aList(f); E.N = D.nPop; E.G0 = P.n_generations; E.cc = R.cc{f};
    E.U = R.winU{f}; E.O = R.winO{f}; E.wt = R.winW{f}; E.omega = R.omega;
    Gt = E.G0 + cfg.nGen;
    [E.uPath, ~, E.opt] = tracking.circle(E.omega, E.cc, Gt, P);
    Th0 = tracking.build_theta(E.a, D.w0(f), P);
    cc = tracking.steady(Th0, zeros(P.total, 1), P);
    E.ccCheck = max(abs(cc([3 4]) - E.cc));
    % the stored window is the last 36 generations of the original run
    gi = (E.G0 - size(E.U, 2) + 1):E.G0;
    % stored as winU = opt - cc (run_tracking_task.m), so rebuild it that way;
    % the continuation itself feeds u_path = pos, exactly as the original loop did
    Uw = zeros(size(E.U)); Uw(1:2, :) = E.opt(:, gi) - E.cc;
    E.winCheck = max(max(abs(E.O - E.opt(:, gi)))) + max(max(abs(E.U - Uw)));
end


% The split
function s = prepare_split(Rr, r, E, P, cfg)
%PREPARE_SPLIT  rebuild the generation-350 population, solve (b, delta) on ThBest, verify.
    t = tic;
    s.r = r;
    [M, s.nDistinct, s.topFreq] = rebuild_population(Rr.popTh, Rr.popFreq, E.N, P);
    s.pop = M;
    s.bestCheck = best_check(M, Rr.ThBest, E, P);
    s.dom = Rr.ThBest;                                 % the best genotype at the split
    gB = struct('Th', s.dom, 'U', E.U, 'wt', E.wt);
    geo0 = geom(s.dom, E.U, E.wt, 0, 1, P, []);
    s.alpha0 = geo0.alpha; s.xn0 = geo0.xn; s.valid0 = geo0.valid;
    [s.knobA, s.okA] = solve_arm(gB, P, cfg, geo0, 'A');
    [s.knobX, s.okX] = solve_arm(gB, P, cfg, geo0, 'X');
    s.shiftA = phenotype_shift(s.dom, s.knobA, E, P);
    s.shiftX = phenotype_shift(s.dom, s.knobX, E, P);
    s.wDom = [w_window(s.dom, [0 0], [0; 0], E, P), ...
              w_window(s.dom, s.knobA, [0; 0], E, P), w_window(s.dom, s.knobA, s.shiftA, E, P), ...
              w_window(s.dom, s.knobX, [0; 0], E, P), w_window(s.dom, s.knobX, s.shiftX, E, P)];
    s.ver = verify_population(M, s, E, P, cfg);
    s.sec = toc(t);
end


function V = verify_population(M, s, E, P, cfg)
%VERIFY_POPULATION per-individual alpha and ||x*|| ratios under each arm's (b, delta),
%   on nVer sampled individuals.
    st = RandStream('twister', 'Seed', cfg.readSeed + s.r);
    idx = randperm(st, size(M, 3), min(cfg.nVer, size(M, 3)));
    V.ra = nan(numel(idx), 2); V.rx = V.ra; V.valid = false(numel(idx), 3);
    knobs = {s.knobA, s.knobX};
    for q = 1:numel(idx)
        Th = M(:, :, idx(q));
        g0 = geom(Th, E.U, E.wt, 0, 1, P, []); V.valid(q, 1) = g0.valid;
        for k = 1:2
            if any(isnan(knobs{k})), continue; end
            g1 = geom(with_knob(Th, knobs{k}), E.U, E.wt, knobs{k}(1), 1, P, []);
            V.valid(q, k+1) = g1.valid;
            V.ra(q, k) = g1.alpha / g0.alpha; V.rx(q, k) = g1.xn / g0.xn;
        end
    end
end


% The (b, delta) solve, copied from run_tracking_single_genotypes
function [z, okAll] = solve_arm(g, P, cfg, geo0, arm)
%SOLVE_ARM  run_tracking_single_genotypes.m's continuation over cfg.contLevels (ending at cfg.m),
%   with its homotopy fallback; a failed intermediate level does not stop the march, as
%   there. Returns the last level's (b, delta), or NaN if Newton or the check from x0
%   failed.
    zz = [0; 0]; mPrev = 1; ok = false;
    for m = cfg.contLevels
        [tA, tX] = arm_targets(arm, m, geo0);
        [zz2, ok] = newton_bd(g, P, tA, tX, zz);
        if ~ok        % homotopy fallback: march from the last solved level in 8 sub-steps
            zs = zz; ok = true;
            for ms = mPrev + (m - mPrev)*(1:8)/8
                [ta, tx] = arm_targets(arm, ms, geo0);
                [zs, oks] = newton_bd(g, P, ta, tx, zs);
                if ~oks, ok = false; break; end
            end
            zz2 = zs;
        end
        if ok
            mPrev = m;
            e = geom(with_knob(g.Th, zz2), g.U, g.wt, zz2(1), 1, P, []);   % check from x0
            ok = e.valid && abs(e.alpha/tA - 1) <= 1e-6 && abs(e.xn/tX - 1) <= 1e-6;
            zz = zz2;
        end
    end
    okAll = ok;
    z = [zz(1), zz(2)];
    if ~okAll, z = [NaN NaN]; end
end


function [z, ok] = newton_bd(g, P, tA, tX, z)
%NEWTON_BD  2-D damped Newton on (b, delta) for log(alpha/tA) = log(xn/tX) = 0.
    ok = false;
    [r, geo] = bd_resid(g, P, z, tA, tX, []);
    if ~geo.valid, return; end
    for it = 1:40
        if max(abs(r)) < 1e-9, ok = true; return; end
        J = zeros(2); hs = 1e-6*max(1, abs(z));
        for k = 1:2
            zp = z; zp(k) = zp(k) + hs(k);
            J(:, k) = (bd_resid(g, P, zp, tA, tX, geo.X) - r) / hs(k);
        end
        s = -J \ r; lam = 1; moved = false;
        for h = 1:30
            zn = z + lam*s;
            [rn, gn] = bd_resid(g, P, zn, tA, tX, geo.X);
            if gn.valid && norm(rn) < norm(r), moved = true; break; end
            lam = lam/2;
        end
        if ~moved, return; end
        z = zn; r = rn; geo = gn;
    end
    ok = max(abs(r)) < 1e-9;
end


function [r, geo] = bd_resid(g, P, z, tA, tX, Xinit)
    Th2 = g.Th; Th2(3,3) = Th2(3,3) + z(2); Th2(4,4) = Th2(4,4) + z(2);
    geo = geom(Th2, g.U, g.wt, z(1), 1, P, Xinit);
    r = [log(max(geo.alpha, realmin)/tA); log(max(geo.xn, realmin)/tX)];
    if ~geo.valid, r = [Inf; Inf]; end
end


% Geometry, copied from run_tracking_single_genotypes
function geo = geom(Th, U, wt, b, gain, P, Xinit)
%GEOM  window-weighted alpha, ||x*|| at the fixed point reached from x0.
    n = P.total; C = [0 0 1 0; 0 0 0 1]; nE = size(U, 2);
    Ue = gain*U + [0; 0; b; b];
    if isempty(Xinit), X = integrate_fp(Th, Ue, P); else, X = Xinit; end
    [X, conv] = polish_fp(Th, Ue, P, X);
    a = 0; xn = 0; mr = -Inf;
    for e = 1:nE
        x = X(:, e); h = Th*x + Ue(:, e); act = h > 0;
        fp = zeros(n, 1); fp(act) = P.K ./ (P.K + h(act)).^2; Df = diag(fp);
        CP = C * ((P.mu*eye(n) - Df*Th) \ Df);
        a = a + wt(e)*norm(CP, 'fro'); xn = xn + wt(e)*norm(x);
        mr = max(mr, max(real(eig(Df*Th - P.mu*eye(n)))));
    end
    on = min(X(3, :)) > 0.5 && min(X(4, :)) > 0.5;
    geo = struct('alpha', a, 'xn', xn, 'maxRe', mr, 'X', X, ...
        'valid', conv && mr < 0 && on && all(isfinite(X(:))));
end


% Continuation after the split
function jobs = run_continuations(splits, E, P, cfg, out)
%RUN_CONTINUATIONS  every (replicate, arm) job, in waves; saved after each wave.
    nR = numel(splits); nA = numel(cfg.arms);
    [ja, jr] = ndgrid(1:nA, 1:nR); jr = jr(:); ja = ja(:);      % waves finish whole replicates
    jobs = cell(numel(jr), 1); t0 = tic;
    for w0 = 1:cfg.waveSize:numel(jr)
        wv = w0:min(w0 + cfg.waveSize - 1, numel(jr)); res = cell(numel(wv), 1);
        spw = splits(jr(wv)); aw = ja(wv);
        parfor q = 1:numel(wv)
            res{q} = continue_arm(spw{q}, aw(q), E, P, cfg);
        end
        jobs(wv) = res;
        out.jobs = jobs; save(cfg.outFile, 'out');
        fprintf('jobs %d-%d of %d done, %.1f min, saved\n', wv(1), wv(end), numel(jr), toc(t0)/60);
    end
end


function z = continue_arm(s, arm, E, P, cfg)
%CONTINUE_ARM  mutation and selection from the split population under one arm's (b, delta).
%   Per generation: mean fitness, the largest fitness w_max, and the mean fitness, in this
%   generation's environment, of this generation's unmutated parents. At the window points
%   (is_window_point) that generation's population is scored over the whole scoring window:
%   its mean window fitness and its best genotype's (for the load carried).
    t = tic;
    [knob, shift, seed] = arm_setup(s, arm, cfg);
    G = cfg.nGen;
    z = struct('r', s.r, 'arm', arm, 'knob', knob, 'shift', shift, 'seed', seed, ...
        'void', any(isnan(knob)) || any(isnan(shift)), 'meanFit', nan(1, G), 'wMax', nan(1, G), ...
        'wPar', nan(1, G), 'xnMean', nan(1, G), 'track', [], 'sec', 0);
    if z.void, return; end
    M = s.pop; aF = E.a + knob(2);
    M(3,3,:) = aF; M(4,4,:) = aF;
    N = size(M, 3); bb = [0; 0; knob(1); knob(1)];
    st = RandStream('twister', 'Seed', cfg.readSeed + 7919*s.r);    % readouts' own stream
    DW = draw_dw(P, cfg.nDraws, cfg.readSeed + 104729*s.r);
    z.track = track_point([], 0, M, knob, shift, E, P, cfg, DW, st);
    rng(seed, 'twister'); Mpar = [];
    for k = 1:G
        gen = E.G0 + k; u = E.uPath(:, gen) + bb; o = E.opt(:, gen) + shift;
        Xe = steady_batch(M, u, P); pe = Xe([3 4], :);
        fit = exp(-P.selection_strength * sum((pe - o).^2, 1))';
        z.meanFit(k) = mean(fit); z.wMax(k) = max(fit);
        z.xnMean(k) = mean(vecnorm(Xe, 2, 1));
        if ~isempty(Mpar)
            Xp = steady_batch(Mpar, u, P);
            z.wPar(k) = mean(exp(-P.selection_strength * sum((Xp([3 4], :) - o).^2, 1)));
        end
        if any(fit > 0); wsel = fit / sum(fit); else; wsel = ones(N,1)/N; end
        ix = randsample(1:N, N, true, wsel);
        M = M(:,:,ix); Mpar = M;
        mask = rand(size(M)) < P.mutation_prob;
        M(mask) = M(mask) + randn(nnz(mask),1) * P.mutation_sd;
        M(3,3,:) = aF; M(4,4,:) = aF;
        if mod(k, cfg.geoEvery) == 0 || is_window_point(k, cfg)
            z.track = track_point(z.track, k, M, knob, shift, E, P, cfg, DW, st);
        end
    end
    z.sec = toc(t);
end


function T = track_point(T, k, M, knob, shift, E, P, cfg, DW, st)
%TRACK_POINT  (i) population medians over nTrack sampled individuals of alpha and ||x*||,
%   and of beta_rms and L0 over the first nBeta of them (every betaEvery); (ii) at window
%   points, the population scored over the whole cycle: mean window fitness wBar, the best
%   genotype's window fitness wBest, and that genotype's alpha and ||x*||. Sampling uses
%   the readouts' own stream st; DW is common across arms and time points.
    t1 = tic;
    idx = randperm(st, size(M, 3), cfg.nTrack);
    a = nan(1, cfg.nTrack); x = a; v = false(1, cfg.nTrack);
    for q = 1:cfg.nTrack
        g = geom(M(:, :, idx(q)), E.U, E.wt, knob(1), 1, P, []);   % diagonal carries a + delta
        a(q) = g.alpha; x(q) = g.xn; v(q) = g.valid;
    end
    e = struct('k', k, 'alpha', median(a), 'xn', median(x), 'validShare', mean(v), ...
        'beta_rms', NaN, 'L0', NaN, 'flip', NaN, 'wBar', NaN, 'wBest', NaN, ...
        'alphaBest', NaN, 'xnBest', NaN, 'secGeo', toc(t1), 'secBeta', 0, 'secWin', 0);
    if mod(k, cfg.betaEvery) == 0
        t1 = tic;
        Ue = E.U + [0; 0; knob(1); knob(1)]; b = nan(3, cfg.nBeta);
        for q = 1:cfg.nBeta
            c = M(:, :, idx(q));
            m = load_of(c, c + DW, Ue, E.O + shift, E.wt, P, cfg.flipDist);
            b(:, q) = [m.beta_rms; m.L0; m.flip];
        end
        mb = median(b, 2); e.beta_rms = mb(1); e.L0 = mb(2); e.flip = mb(3);
        e.secBeta = toc(t1);
    end
    if is_window_point(k, cfg)
        t1 = tic;
        [e.wBar, e.wBest, ib] = window_eval(M, knob, shift, E, P);
        g = geom(M(:, :, ib), E.U, E.wt, knob(1), 1, P, []);
        e.alphaBest = g.alpha; e.xnBest = g.xn;
        e.secWin = toc(t1);
    end
    if isempty(T), T = e; else, T(end+1) = e; end
end


function DW = draw_dw(P, M, seed)
%DRAW_DW  M draws from the task's mutation distribution (W33 and W44 not mutated), as in
%   run_tracking_single_genotypes, on their own stream.
    st = RandStream('twister', 'Seed', seed); n = P.total; DW = zeros(n, n, M);
    for k = 1:M
        dW = zeros(n); mask = rand(st, n) < P.mutation_prob;
        dW(mask) = randn(st, nnz(mask), 1) * P.mutation_sd;
        dW(3,3) = 0; dW(4,4) = 0;
        DW(:, :, k) = dW;
    end
end


function m = load_of(Th, TA, Uu, Ou, wt, P, flipDist)
%LOAD_OF  run_tracking_single_genotypes's load_of: the per-mutation load L0, the shift per mutation
%   beta_rms, and the share of draws that move an output by more than flipDist.
    M = size(TA,3); s = P.selection_strength; nE = size(Uu,2);
    accW0 = zeros(M,1); accD2 = zeros(M,1); maxD = zeros(M,1); accP = 0;
    for g = 1:nE
        xp = tracking.steady(Th, Uu(:,g), P); pp = xp([3 4]);
        accP = accP + wt(g) * exp(-s * sum((pp - Ou(:,g)).^2));
        X = steady_batch(TA, Uu(:,g), P);
        dPar = X([3 4],:) - pp;
        accW0 = accW0 + wt(g) * exp(-s * sum(dPar.^2,1))';
        accD2 = accD2 + wt(g) * sum(dPar.^2,1)';
        maxD  = max(maxD, sqrt(sum(dPar.^2,1))');
    end
    m.w_par = accP; m.L0 = 1 - mean(accW0);
    m.beta_rms = sqrt(mean(accD2)); m.flip = mean(maxD > flipDist);
end


% Checks and reports
function S1 = lowering_settings(cfg)
%LOWERING_SETTINGS  The single-genotype run's (b, delta) at m = cfg.m for the moving-target
%   genotypes of this design.
    o = getfield(load(cfg.stage1File), 'out');
    k = find([o.G.ai] == 1 & [o.G.f] == cfg.designIndex);
    lv = find(abs(o.cfg.levels - cfg.m) < 1e-12);
    S1 = struct('r', [o.G(k).r], 'A', nan(numel(k), 2), 'X', nan(numel(k), 2));
    for q = 1:numel(k)
        z = o.res{k(q)};
        S1.A(q, :) = [z.A(lv).b, z.A(lv).delta]; S1.X(q, :) = [z.X(lv).b, z.X(lv).delta];
    end
end


function c1 = report_c1(splits, S1, E, cfg)
%REPORT_C1  The check against the stored run: environment; population rebuild -> ThBest; knob solver -> lowering_settings.mat.
    fprintf('\n---- CHECK AGAINST THE STORED RUN ----\n');
    sp = [splits{:}];
    fprintf('founder rest state vs stored cc max|diff| %.3g; stored window vs continued circle %.3g (both must be 0)\n', ...
        E.ccCheck, E.winCheck);
    dB = max([sp.bestCheck]);
    dK = 0; nOK = 0;
    for q = 1:numel(sp)
        i = find(S1.r == sp(q).r);
        if sp(q).okA && sp(q).okX
            nOK = nOK + 1;
            dK = max([dK, abs(sp(q).knobA ./ S1.A(i, :) - 1), abs(sp(q).knobX ./ S1.X(i, :) - 1)]);
        end
    end
    c1.bestCheck = dB; c1.knobRel = dK; c1.nOK = nOK;
    c1.pass = dB == 0 && nOK == numel(sp) && dK <= 1e-6 && E.ccCheck == 0 && E.winCheck == 0;
    fprintf('best-individual search on the rebuilt population -> ThBest: max|diff| = %.3g over %d replicates (must be 0)\n', dB, numel(sp));
    fprintf('knob solver on ThBest -> the (b, delta) of lowering_settings.mat at m = %.1f: %d/%d solved, max rel diff %.3g (<= 1e-6)\n', ...
        cfg.m, nOK, numel(sp), dK);
    fprintf('population structure at the split: distinct wirings median %d of 1000, largest clone median %.4f\n', ...
        round(median([sp.nDistinct])), median([sp.topFreq]));
    fprintf('CHECK AGAINST THE STORED RUN %s\n', ternary(c1.pass, 'PASS', 'FAIL'));
end


function V = report_verification(splits, cfg)
%REPORT_VERIFICATION the per-individual verification (a check that can fail), pooled and per replicate.
    sp = [splits{:}]; nm = {'ALPHA-DOWN', 'XNORM-DOWN'};
    fprintf('\n---- KNOB VERIFICATION AT THE SPLIT ----\n');
    fprintf('dominant genotype (ThBest) valid: %d/%d | knobs solved A %d/%d, X %d/%d\n', ...
        sum([sp.valid0]), numel(sp), sum([sp.okA]), numel(sp), sum([sp.okX]), numel(sp));
    V.rep = false(numel(sp), 2);
    for k = 1:2
        ra = []; rx = []; vf = nan(1, numel(sp));
        for q = 1:numel(sp)
            v = sp(q).ver; ok = v.valid(:, 1) & v.valid(:, k+1);
            a = v.ra(ok, k); x = v.rx(ok, k); ra = [ra; a]; rx = [rx; x]; %#ok<AGROW>
            vf(q) = mean(ok);
            V.rep(q, k) = ~isempty(a) && rep_pass(median(a), median(x), k, cfg);
        end
        V.pass(k) = ~isempty(ra) && rep_pass(median(ra), median(rx), k, cfg);
        V.medA(k) = median(ra); V.medX(k) = median(rx); V.validShare(k) = median(vf);
        [tgt, oth] = target_other(V.medA(k), V.medX(k), k);
        fprintf('%-10s pooled over %d sampled individuals: median alpha ratio %.4f, ||x*|| ratio %.4f | target factor %.4f (target %.2f within 20%%), other factor move %.4f (< 0.20) | valid share median %.3f | replicates failing %d/%d | %s\n', ...
            nm{k}, numel(ra), V.medA(k), V.medX(k), tgt, cfg.m, abs(oth - 1), V.validShare(k), ...
            sum(~V.rep(:, k)), numel(sp), ternary(V.pass(k), 'KEPT', 'NOT KEPT'));
        fprintf('           spread of per-individual ratios: alpha IQR [%.4f, %.4f], ||x*|| IQR [%.4f, %.4f]\n', ...
            prctile(ra, 25), prctile(ra, 75), prctile(rx, 25), prctile(rx, 75));
    end
    kb = reshape([sp.knobA], 2, [])'; kx = reshape([sp.knobX], 2, [])';
    wd = reshape([sp.wDom], 5, [])';
    fprintf('knobs on ThBest: A b = %.3f [%.3f, %.3f] delta = %.3f | X b = %.3f [%.3f, %.3f] delta = %.3f (median [range])\n', ...
        median(kb(:,1)), min(kb(:,1)), max(kb(:,1)), median(kb(:,2)), median(kx(:,1)), min(kx(:,1)), max(kx(:,1)), median(kx(:,2)));
    fprintf('ThBest window tracking fitness (median): control %.4f | A unshifted %.4f, re-centred %.4f | X unshifted %.4f, re-centred %.4f\n', ...
        median(wd(:,1)), median(wd(:,2)), median(wd(:,3)), median(wd(:,4)), median(wd(:,5)));
    V.knobA = kb; V.knobX = kx; V.wDom = wd;
end


function Pl = timing_run_tolerance(out, cfg)
%TIMING_RUN_TOLERANCE  the reproduction tolerance of the control check from the timing run's two controls.
    J = [out.jobs{:}]; J = J(~[J.void]);
    Pl.tol = seed_tolerance(J, cfg);
end


function T = seed_tolerance(J, cfg)
%SEED_TOLERANCE  CONTROL vs CONTROL2 per replicate over the last loadWin generations.
    rs = unique([J.r]); d = nan(numel(rs), 2);
    for q = 1:numel(rs)
        c1 = J([J.r] == rs(q) & [J.arm] == 1); c2 = J([J.r] == rs(q) & [J.arm] == 2);
        [l1, m1] = window_stats(c1, cfg); [l2, m2] = window_stats(c2, cfg);
        d(q, :) = [m2 - m1, l2 - l1];
    end
    T.diff = d; T.sdMaint = std(d(:, 1)); T.sdLoad = std(d(:, 2));
    T.tolerance = 3 * sqrt(mean(d(:, 1).^2));
    fprintf('\n---- seed-vs-seed (CONTROL2 - CONTROL), last %d generations ----\n', cfg.loadWin);
    fprintf('per replicate mean-fitness diff %s | load diff %s\n', mat2str(d(:,1)', 4), mat2str(d(:,2)', 4));
    fprintf('reproduction tolerance of the control check = 3 x RMS(mean-fitness diff) = %.4f\n', T.tolerance);
end


function [load, maint, benefit, mutCost] = window_stats(z, cfg)
%WINDOW_STATS  per replicate-arm, over the last loadWin generations:
%   load     = 1 - sum(wBar) / sum(wBest) over the window points in it (the load carried);
%   maint    = mean per-generation mean fitness;
%   benefit  = mean wBest (the best genotype's own tracking fitness);
%   mutCost  = 1 - sum(mean w) / sum(mean w of the unmutated parents), per generation.
    G = numel(z.meanFit); idx = G - cfg.loadWin + 1 : G;
    tk = z.track([z.track.k] > G - cfg.loadWin & isfinite([z.track.wBest]));
    load = 1 - sum([tk.wBar]) / sum([tk.wBest]);
    maint = mean(z.meanFit(idx)); benefit = mean([tk.wBest]);
    mutCost = 1 - sum(z.meanFit(idx)) / sum(z.wPar(idx));
end


% Analysis of the main run
function A = analyse(out, Rs, cfg)
    J = [out.jobs{:}]; rs = cfg.repList; nR = numel(rs);
    Y = nan(nR, 4, 4);                     % replicate x arm x [load maint benefit mutCost]
    for q = 1:nR
        for a = 1:4
            z = J([J.r] == rs(q) & [J.arm] == a);
            if isempty(z) || z.void, continue; end
            [Y(q,a,1), Y(q,a,2), Y(q,a,3), Y(q,a,4)] = window_stats(z, cfg);
        end
    end
    A.Y = Y; A.c1c = check_c1c(Y, Rs, cfg);
    keep = out.verif.rep;                  % per-replicate verification pass, arms A and X
    A.table = arm_table(Y, keep, cfg);
    A.push = report_added_input(J, rs, cfg);
    A.stationarity = report_stationarity(J, rs, cfg);
end


function c = check_c1c(Y, Rs, cfg)
% The check that the control's mean fitness over the last loadWin generations reproduces the
% stored mean fitness of the same replicates before the split (generations 291-350), within
% cfg.tolerance: 3 times the RMS of the second-control-minus-control differences in mean
% fitness, as printed by a short run with `cfg.pilot` = true (seed_tolerance).
    d = Y(:, 1, 2)' - [Rs.maint];
    c.meanDiff = mean(d, 'omitnan'); c.tol = cfg.tolerance;
    c.pass = abs(c.meanDiff) <= c.tol;
    fprintf('\nCONTROL CHECK: CONTROL mean fitness (last %d gens) - stored maint (gens 291-350): mean %.4f over %d replicates, tolerance %.4f -> %s\n', ...
        cfg.loadWin, c.meanDiff, sum(isfinite(d)), c.tol, ternary(c.pass, 'PASS', 'FAIL'));
end


function T = arm_table(Y, keep, cfg)
%ARM_TABLE  paired differences vs CONTROL, replicate bootstrap, per outcome.
    nm = {'LOAD', 'maintenance', 'benefit (best genotype window w)', 'mutation cost per generation'};
    st = RandStream('twister', 'Seed', cfg.bootSeed);
    nR = size(Y, 1); B = randi(st, nR, nR, cfg.nBoot);
    fprintf('\n---- ARM TABLE: level (median over replicates) and paired difference vs CONTROL, mean [95%% bootstrap CI], sign test ----\n');
    for o = 1:4
        fprintf('%s: CONTROL %.4f', nm{o}, median(Y(:,1,o), 'omitnan'));
        for a = 2:4
            ok = isfinite(Y(:,a,o)) & isfinite(Y(:,1,o));
            if a >= 3, ok = ok & keep(:, a-2); end
            d = Y(:,a,o) - Y(:,1,o); d(~ok) = NaN;
            bm = mean(d(B), 1, 'omitnan');
            ci = prctile(bm, [2.5 97.5]); p = sign_p(d(ok));
            T(o, a) = struct('level', median(Y(ok,a,o)), 'd', mean(d, 'omitnan'), 'ci', ci, ...
                'n', sum(ok), 'p', p, 'boot', bm); %#ok<AGROW>
            fprintf(' | %s %.4f, d %+.4f [%+.4f, %+.4f] n %d sign-p %.3g', cfg.arms{a}, ...
                T(o,a).level, T(o,a).d, ci(1), ci(2), sum(ok), p);
        end
        fprintf('\n');
    end
    r = T(1,3).boot ./ T(1,4).boot;
    fprintf('ratio of load effects ALPHA-DOWN / XNORM-DOWN: %.3f [%.3f, %.3f]\n', ...
        T(1,3).d / T(1,4).d, prctile(r, 2.5), prctile(r, 97.5));
    sd0 = std(Y(:,2,1) - Y(:,1,1), 'omitnan');
    fprintf('null arm CONTROL2 - CONTROL load SD %.4f -> minimum detectable effect (80%% power, 5%%) ~ %.4f\n', ...
        sd0, 2.8 * sd0 / sqrt(nR));
end


function S = report_stationarity(J, rs, cfg)
    fprintf('\n---- stationarity: the load at the first and second half of the load-window points (median over replicates) ----\n');
    G = numel(J(1).meanFit); h = G - cfg.loadWin/2;
    S = nan(4, 2);
    for a = 1:4
        v = nan(numel(rs), 2);
        for q = 1:numel(rs)
            z = J([J.r] == rs(q) & [J.arm] == a);
            if isempty(z) || z.void, continue; end
            tk = z.track(isfinite([z.track.wBest]) & [z.track.k] > G - cfg.loadWin);
            k = [tk.k]; w1 = k <= h; w2 = k > h;
            v(q, :) = [1 - sum([tk(w1).wBar])/sum([tk(w1).wBest]), 1 - sum([tk(w2).wBar])/sum([tk(w2).wBest])];
        end
        S(a, :) = median(v, 1, 'omitnan');
        fprintf('%-10s %.4f -> %.4f\n', cfg.arms{a}, S(a, 1), S(a, 2));
    end
end


% Utilities
function start_pool(n)
    pl = gcp('nocreate');
    if ~isempty(pl) && pl.NumWorkers == n, return; end
    if ~isempty(pl), delete(pl); end
    parpool('Processes', n);
    pids = fetchOutputs(parfevalOnAll(@() feature('getpid'), 1));
    if ispc
        for k = 1:numel(pids)      % this run's own workers only, by PID
            system(sprintf(['powershell -NoProfile -Command "(Get-Process -Id %d).PriorityClass' ...
                ' = ''BelowNormal''"'], pids(k)));
        end
        fprintf('pool of %d workers, PIDs %s, each set to BelowNormal\n', n, mat2str(pids'));
    else
        fprintf('pool of %d workers, PIDs %s, priority not changed on this system\n', n, mat2str(pids'));
    end
    wait(parfevalOnAll(@() maxNumCompThreads(1), 0));
end

function diary_file(f, on)
    if on
        if exist(f, 'file'), delete(f); end
        diary(f); diary on;
    else
        diary off;
    end
end
