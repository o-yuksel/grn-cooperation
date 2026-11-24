function out = run_tracking_single_genotypes(cfg)
% RUN_TRACKING_SINGLE_GENOTYPES  On single genotypes of the tracking task, how does the
% per-mutation load change as responsiveness alpha or ||x*|| is lowered with the other held?
% Lowers each through four levels and reads the slope of log load on log factor.
%
% Produces : Methods, "Non-social tracking task" (the input-gain test); main text,
%            "Responsiveness carries a mutation load in a non-social trait"; SI S5.1; Table S14
%            (the single-genotype rows, printed in lowering_settings.txt); lowering_settings.mat is also read by
%            run_tracking_moving_target and run_tracking_held_still, which check their own
%            (b, delta) against it
% Inputs   : publication/tracking_runs/tracking_runs_N1000.mat (the generation-350
%              best genotypes, written by run_tracking_task)
%            publication/tracking_runs/tracking_runs_N1000_load.mat (their
%              per-mutation loads, written by tracking_load_vs_identity; for a reproduction check)
% Outputs  : publication/tracking_single_genotypes/lowering_settings.mat and lowering_settings.txt
% Run      : run_tracking_single_genotypes()   (about 20 min on 3 workers)
% Seed     : genotype (target t: 1 moving, 2 held still; design index f; replicate r) draws its
%            mutations once, from rng(20260928 + 1000 f + r + 100000 (t - 1), 'twister'), for every
%            level; the reproduction check uses tracking_load_vs_identity's 20260917 + 1000 f + r;
%            the bootstrap over genotypes its own stream, seed 20260929
%
% The 240 genotypes are the generation-350 best genotypes of 4 designs x 2 targets x 30
% replicates; nothing evolves. The output bias b (a constant input to genes 3 and 4) and delta
% (added to W33 = W44, which mutation does not touch) move the factors: at each level m (0.9 to
% 0.6) a damped Newton solve sets alpha to m times its value at fixed ||x*|| (arm A) or the
% reverse (arm X), and the per-mutation load L0 is measured on 400 draws. alpha is the response
% of outputs 3 and 4 to an input at any of the four genes, averaged over the scoring window.
% Arm G scales the input gain g of the cue channels (0.6 to 1.5) on the moving-target genotypes;
% report_gain prints whether this also moves alpha.

    cfg = defaults(cfg);
    setup_process();
    D = getfield(load(cfg.inFile), 'out');
    P = D.P;
    G = collect_genotypes(D, cfg);
    clear D
    ref = getfield(load(cfg.loadRefFile), 'out');

    diary_file(cfg.txtFile, true);
    fprintf('\n=========== TRACKING TASK, SINGLE GENOTYPES ===========\n');
    fprintf('source %s\n%d genotypes | %d draws | levels %s | gains %s\n', cfg.inFile, ...
        numel(G), cfg.nDraws, mat2str(cfg.levels), mat2str(cfg.gains));
    fprintf('mutation distribution p = %.2f, sd = %.3f, s = %.1f, self-activations (3,3), (4,4) not mutated\n', ...
        P.mutation_prob, P.mutation_sd, P.selection_strength);

    out = struct('cfg', cfg, 'P', P);
    out.routing = routing_check(G, P, ref, cfg);
    start_pool(cfg.nWorkers);
    res = cell(numel(G), 1);
    t0 = tic;
    for f = cfg.designList
        idx = find([G.f] == f);
        Gs = G(idx); rs = cell(numel(idx), 1);
        parfor q = 1:numel(idx)
            rs{q} = one_genotype(Gs(q), P, cfg);
        end
        res(idx) = rs;
        out.res = res; out.G = strip_genotypes(G);
        save(cfg.outFile, 'out');
        fprintf('[a = %.2f] %d genotypes done, %.1f min elapsed, saved.\n', ...
            G(idx(1)).a, numel(idx), toc(t0)/60);
    end
    out.elapsedMin = toc(t0)/60;
    out.analysis = analyse(res, G, cfg);
    out.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    save(cfg.outFile, 'out');
    fprintf('\nSaved %s (%.1f min)\n', cfg.outFile, out.elapsedMin);
    diary_file(cfg.txtFile, false);
end


% Configuration
function cfg = defaults(cfg)
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    pub = fullfile(grn_data_file('publication'));
    d = struct('inFile', fullfile(pub, 'tracking_runs', 'tracking_runs_N1000.mat'), ...
        'loadRefFile', fullfile(pub, 'tracking_runs', 'tracking_runs_N1000_load.mat'), ...
        'outDir', fullfile(pub, 'tracking_single_genotypes'), 'nDraws', 400, ...
        'levels', [0.9 0.8 0.7 0.6], 'gains', [0.6 0.8 1.25 1.5], 'nWorkers', 3, ...
        'designList', 1:4, 'repList', 1:30, 'seed', 20260928, 'loadSeed', 20260917, ...
        'routeDraws', 2000, 'nBoot', 10000, 'bootSeed', 20260929, ...
        'otherTol', 0.20, 'fdTol', 1e-4, 'flipDist', 0.5);
    fn = fieldnames(d);
    for k = 1:numel(fn)
        if ~isfield(cfg, fn{k}), cfg.(fn{k}) = d.(fn{k}); end
    end
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    if ~isfield(cfg, 'outFile'), cfg.outFile = fullfile(cfg.outDir, 'lowering_settings.mat'); end
    if ~isfield(cfg, 'txtFile'), cfg.txtFile = fullfile(cfg.outDir, 'lowering_settings.txt'); end
end


function G = collect_genotypes(D, cfg)
    arms = {'moving', 'fixed'}; G = struct([]);
    for ai = 1:2
        R = D.(arms{ai});
        for f = cfg.designList
            for r = cfg.repList
                k = numel(G) + 1;
                G(k).Th = R.ThBest{f, r}; G(k).U = R.winU{f}; G(k).O = R.winO{f};
                G(k).wt = R.winW{f}; G(k).ai = ai; G(k).f = f; G(k).r = r;
                G(k).a = D.aList(f);
            end
        end
    end
end


function S = strip_genotypes(G)
    S = struct('ai', {G.ai}, 'f', {G.f}, 'r', {G.r}, 'a', {G.a}, 'Th', {G.Th});
end


% One genotype
function z = one_genotype(g, P, cfg)
%ONE_GENOTYPE  baseline, arm A, arm X (and arm G with the target moving), common draws.
    t = tic;
    geo0 = geom(g.Th, g.U, g.wt, 0, 1, P, []);
    DW = draw_dw(P, cfg.nDraws, cfg.seed + 1000*g.f + g.r + 100000*(g.ai - 1));
    z.base = level_eval(g.Th, 0, 0, 1, g, P, DW, cfg);
    z.alpha0 = geo0.alpha; z.xn0 = geo0.xn;
    z.A = solve_arm(g, P, DW, cfg, geo0, 'A');
    z.X = solve_arm(g, P, DW, cfg, geo0, 'X');
    z.G = [];
    if g.ai == 1
        for k = 1:numel(cfg.gains)
            z.G = cat_struct(z.G, level_eval(g.Th, 0, 0, cfg.gains(k), g, P, DW, cfg));
        end
    end
    z.sec = toc(t);
end


function L = solve_arm(g, P, DW, cfg, geo0, arm)
%SOLVE_ARM  continuation over the levels, each warm-started at the previous (b, delta).
    zz = [0; 0]; L = []; mPrev = 1;
    for k = 1:numel(cfg.levels)
        m = cfg.levels(k);
        [tA, tX] = arm_targets(arm, m, geo0);
        [zz2, ok, it] = newton_bd(g, P, tA, tX, zz);
        if ~ok        % homotopy fallback: march from the last solved level in 8 sub-steps
            zs = zz; ok = true;
            for ms = mPrev + (m - mPrev)*(1:8)/8
                [ta, tx] = arm_targets(arm, ms, geo0);
                [zs, oks, it] = newton_bd(g, P, ta, tx, zs);
                if ~oks, ok = false; break; end
            end
            zz2 = zs;
        end
        if ok, mPrev = m; end
        if ok
            e = level_eval(g.Th, zz2(1), zz2(2), 1, g, P, DW, cfg);
            % check from x0: the state development reaches must be the solved one
            ok = abs(e.alpha/tA - 1) <= 1e-6 && abs(e.xn/tX - 1) <= 1e-6;
            zz = zz2;
        else
            e = level_eval(g.Th, NaN, NaN, 1, g, P, DW, cfg);
        end
        e.m = m; e.newtonOK = ok; e.newtonIt = it;
        e.tA = tA; e.tX = tX;
        L = cat_struct(L, e);
    end
end


function [z, ok, it] = newton_bd(g, P, tA, tX, z)
%NEWTON_BD  2-D damped Newton on (b, delta) for log(alpha/tA) = log(xn/tX) = 0.
    ok = false;
    [r, geo] = bd_resid(g, P, z, tA, tX, []);
    if ~geo.valid, it = 0; return; end
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


% One level: geometry and load
function e = level_eval(Th, b, delta, gain, g, P, DW, cfg)
%LEVEL_EVAL  verification geometry (from x0) and the load under the common draws.
    e = empty_level();
    e.b = b; e.delta = delta; e.gain = gain;
    if isnan(b), return; end
    Th2 = Th; Th2(3,3) = Th2(3,3) + delta; Th2(4,4) = Th2(4,4) + delta;
    geo = geom(Th2, g.U, g.wt, b, gain, P, []);
    e.alpha = geo.alpha; e.xn = geo.xn; e.alpha_cue = gain*geo.alpha_cue;
    e.maxRe = geo.maxRe; e.hmin = geo.hmin; e.valid = geo.valid;
    e.fdErr = fd_check(Th2, g.U(:,1)*gain + [0;0;b;b], P, geo.CP{1}, geo.X(:,1));
    if ~geo.valid, return; end
    Ueff = gain*g.U + [0; 0; b; b];
    TA = Th2 + DW;
    m = load_of(Th2, TA, Ueff, g.O, g.wt, P, cfg.flipDist);
    e.L = m.L; e.L0 = m.L0; e.w_par = m.w_par; e.beta_rms = m.beta_rms;
    e.flip = m.flip;
    ll = lin_load(TA, Th2, geo, g.wt, P);
    e.L_identLin = ll.L; e.beta_lin_rms = ll.beta_rms;
end


function e = empty_level()
    e = struct('b', NaN, 'delta', NaN, 'gain', NaN, 'alpha', NaN, 'xn', NaN, ...
        'alpha_cue', NaN, 'maxRe', NaN, 'hmin', NaN, 'valid', false, 'fdErr', NaN, ...
        'L', NaN, 'L0', NaN, 'w_par', NaN, 'beta_rms', NaN, 'flip', NaN, ...
        'L_identLin', NaN, 'beta_lin_rms', NaN, 'm', NaN, 'newtonOK', false, ...
        'newtonIt', NaN, 'tA', NaN, 'tX', NaN);
end


function DW = draw_dw(P, M, seed)
%DRAW_DW  the mutation distribution of tracking_load_vs_identity's draw_mutants, as
%   perturbations Tm - Th, with W33 and W44 not mutated.
    rng(seed, 'twister'); n = P.total; DW = zeros(n, n, M);
    for k = 1:M
        dW = zeros(n); mask = rand(n) < P.mutation_prob;
        dW(mask) = randn(nnz(mask), 1) * P.mutation_sd;
        dW(3,3) = 0; dW(4,4) = 0;
        DW(:, :, k) = dW;
    end
end


% Geometry
function geo = geom(Th, U, wt, b, gain, P, Xinit)
%GEOM  window-weighted alpha, alpha_cue, ||x*|| at the fixed point reached from x0
%   (or, inside Newton's Jacobian and line search, polished from Xinit).
    n = P.total; C = [0 0 1 0; 0 0 0 1]; nE = size(U, 2);
    Ue = gain*U + [0; 0; b; b];
    if isempty(Xinit), X = integrate_fp(Th, Ue, P); else, X = Xinit; end
    [X, conv] = polish_fp(Th, Ue, P, X);
    a = 0; ac = 0; xn = 0; mr = -Inf; hmin = Inf; CPs = cell(1, nE);
    for e = 1:nE
        x = X(:, e); h = Th*x + Ue(:, e); act = h > 0;
        fp = zeros(n, 1); fp(act) = P.K ./ (P.K + h(act)).^2; Df = diag(fp);
        CP = C * ((P.mu*eye(n) - Df*Th) \ Df);
        a = a + wt(e)*norm(CP, 'fro'); ac = ac + wt(e)*norm(CP(:, 1:2), 'fro');
        xn = xn + wt(e)*norm(x);
        mr = max(mr, max(real(eig(Df*Th - P.mu*eye(n)))));
        hmin = min(hmin, min(abs(h))); CPs{e} = CP;
    end
    on = min(X(3, :)) > 0.5 && min(X(4, :)) > 0.5;
    geo = struct('alpha', a, 'alpha_cue', ac, 'xn', xn, 'maxRe', mr, 'hmin', hmin, ...
        'X', X, 'valid', conv && mr < 0 && on && all(isfinite(X(:))));
    geo.CP = CPs;
end


function err = fd_check(Th, u, P, CP, x0)
%FD_CHECK  central difference of the readout on all four input channels vs C*Phi.
    d = 1e-6; FD = zeros(2, P.total);
    for i = 1:P.total
        up = u; up(i) = up(i) + d; um = u; um(i) = um(i) - d;
        xp = polish_fp(Th, up, P, x0); xm = polish_fp(Th, um, P, x0);
        FD(:, i) = (xp([3 4]) - xm([3 4])) / (2*d);
    end
    err = norm(CP - FD, 'fro') / max(norm(FD, 'fro'), eps);
end


% Load, copied from tracking_load_vs_identity
function m = load_of(Th, TA, Uu, Ou, wt, P, flipDist)
%LOAD_OF  tracking_load_vs_identity.m's load_of, verbatim, plus the flip count.
    M = size(TA,3); s = P.selection_strength; nE = size(Uu,2);
    accW = zeros(M,1); accW0 = zeros(M,1); accP = 0; accE = 0;
    accD2 = zeros(M,1); accD = zeros(M,1); maxD = zeros(M,1);
    for g = 1:nE
        xp = tracking.steady(Th, Uu(:,g), P); pp = xp([3 4]);
        ep = pp - Ou(:,g);
        accP = accP + wt(g) * exp(-s * sum(ep.^2));
        accE = accE + wt(g) * sum(ep.^2);
        X = steady_batch(TA, Uu(:,g), P);
        dOpt = X([3 4],:) - Ou(:,g);
        dPar = X([3 4],:) - pp;
        accW  = accW  + wt(g) * exp(-s * sum(dOpt.^2,1))';
        accW0 = accW0 + wt(g) * exp(-s * sum(dPar.^2,1))';
        accD2 = accD2 + wt(g) * sum(dPar.^2,1)';
        accD  = accD  + wt(g) * sqrt(sum(dPar.^2,1))';
        maxD  = max(maxD, sqrt(sum(dPar.^2,1))');
    end
    m.w_par = accP; m.err_rms = sqrt(accE);
    m.Ew = mean(accW); m.Ew0 = mean(accW0);
    m.L = 1 - m.Ew / m.w_par;
    m.L0 = 1 - m.Ew0;
    m.beta_rms = sqrt(mean(accD2));
    m.beta_mean = mean(accD);
    m.flip = mean(maxD > flipDist);
end


function r = lin_load(TA, Th, geo, wt, P)
%LIN_LOAD  copied from tracking_load_vs_identity (g.CP / g.x read from this file's
%   geometry): the load the identity's first-order shift predicts.
    M = size(TA,3); s = P.selection_strength; nE = numel(geo.CP);
    accW = zeros(M,1); accD2 = zeros(M,1);
    for e = 1:nE
        CP = geo.CP{e}; xs = geo.X(:, e);
        d2 = zeros(M,1);
        for k = 1:M
            dp = CP * ((TA(:,:,k) - Th) * xs);
            d2(k) = sum(dp.^2);
        end
        accW  = accW  + wt(e) * exp(-s * d2);
        accD2 = accD2 + wt(e) * d2;
    end
    r.L = 1 - mean(accW);
    r.beta_rms = sqrt(mean(accD2));
end


% Checks
function rt = routing_check(G, P, ref, cfg)
%ROUTING_CHECK  this file's load at b = delta = 0, g = 1 against the stored loads of
%   tracking_load_vs_identity.
    fprintf('\n---- CHECK AGAINST THE STORED RUN: %s ----\n', cfg.loadRefFile);
    arms = {'moving', 'fixed'}; dL = 0; dL0 = 0; ra = 0; rx = 0; va = 0; vx = 0;
    for k = 1:numel(G)
        g = G(k); S = ref.S.(arms{g.ai});
        geo = geom(g.Th, g.U, g.wt, 0, 1, P, []);
        gv = identity_beta(g.Th, g.U, g.wt, P);          % tracking_load_vs_identity's own path
        va = max(va, abs(gv.alpha - S.alpha_op(g.f, g.r)));
        vx = max(vx, abs(gv.xnorm - S.xnorm_op(g.f, g.r)));
        ra = max(ra, abs(geo.alpha/gv.alpha - 1));
        rx = max(rx, abs(geo.xn/gv.xnorm - 1));
        if g.r ~= cfg.repList(1), continue; end
        rng(cfg.loadSeed + 1000*g.f + g.r, 'twister');
        TA = draw_mutants(g.Th, P, cfg.routeDraws);
        m = load_of(g.Th, TA, g.U, g.O, g.wt, P, cfg.flipDist);
        dL = max(dL, abs(m.L - S.L(g.f, g.r))); dL0 = max(dL0, abs(m.L0 - S.L0(g.f, g.r)));
    end
    rt = struct('maxAbsDiffL', dL, 'maxAbsDiffL0', dL0, 'verbAlpha', va, 'verbXn', vx, ...
        'maxRelAlpha', ra, 'maxRelXn', rx);
    rt.pass = dL == 0 && dL0 == 0 && va == 0 && vx == 0 && ra <= 1e-6 && rx <= 1e-6;
    fprintf('L, L0 (replicate %d of each cell, %d draws, the load file''s seed): max|diff| = %.3g, %.3g (must be 0)\n', ...
        cfg.repList(1), cfg.routeDraws, dL, dL0);
    fprintf('verbatim identity_beta against alpha_op, xnorm_op, all %d genotypes: max|diff| = %.3g, %.3g (must be 0)\n', ...
        numel(G), va, vx);
    fprintf('this file''s Newton-polished geometry against the verbatim one: max rel diff = %.3g, %.3g (<= 1e-6)\n', ...
        ra, rx);
    fprintf('CHECK AGAINST THE STORED RUN %s\n', ternary(rt.pass, 'PASS', 'FAIL'));
end


function g = identity_beta(Th, Uu, wt, P)
%IDENTITY_BETA  copied from tracking_load_vs_identity (only alpha and xnorm are used).
    n = P.total; C = [0 0 1 0; 0 0 0 1];
    mutable = true(n); mutable(3,3) = false; mutable(4,4) = false;
    acc = 0; aacc = 0; xacc = 0;
    nE = size(Uu,2);
    for g = 1:nE
        x = steady_fp(Th, Uu(:,g), P);
        h = Th*x + Uu(:,g); act = h > 0;
        fp = zeros(n,1); fp(act) = P.K ./ (P.K + h(act)).^2;
        Df = diag(fp);
        Phi = (P.mu*eye(n) - Df*Th) \ Df;
        CP = C*Phi;
        ca = sum(CP.^2, 1)'; xb = x.^2;
        a2 = 0;
        for k = 1:n, a2 = a2 + ca(k) * sum(xb(mutable(k,:))); end
        acc  = acc  + wt(g) * P.mutation_prob * P.mutation_sd^2 * a2;
        aacc = aacc + wt(g) * norm(CP,'fro');
        xacc = xacc + wt(g) * norm(x);
    end
    g = struct('beta', sqrt(acc), 'alpha', aacc, 'xnorm', xacc);
end


% Analysis
function A = analyse(res, G, cfg)
    A = struct();
    A.verif = report_verification(res, G, cfg);
    [sA, sX, keepA, keepX] = genotype_slopes(res, G, cfg, 'L0');
    A.sA = sA; A.sX = sX;
    A.primary = pooled(sA, sX, G, cfg, 'L0', keepA, keepX);
    fprintf('\n---- also reported: the same slopes on other load readouts ----\n');
    for s = {'L_identLin', 'beta_rms', 'L'}
        [a2, x2] = genotype_slopes(res, G, cfg, s{1});
        A.secondary.(s{1}) = pooled(a2, x2, G, cfg, s{1}, keepA, keepX);
    end
    A.cells = report_cells(sA, sX, G);
    A.gain = report_gain(res, G, cfg);
end


function [sA, sX, keepA, keepX] = genotype_slopes(res, G, cfg, met)
%GENOTYPE_SLOPES  per genotype OLS slope of log(met) on log(target factor), valid levels.
    nG = numel(G); sA = nan(nG, 1); sX = sA;
    for k = 1:nG
        z = res{k};
        sA(k) = arm_slope(z.base, z.A, 'alpha', 'xn', z.xn0, met, cfg);
        sX(k) = arm_slope(z.base, z.X, 'xn', 'alpha', z.alpha0, met, cfg);
    end
    V = verification_table(res, cfg);
    keepA = V.keepA; keepX = V.keepX;
end


function s = arm_slope(base, lev, fac, other, other0, met, cfg)
    ok = [lev.newtonOK] & [lev.valid] & [lev.fdErr] <= cfg.fdTol & ...
         abs([lev.(other)]/other0 - 1) <= cfg.otherTol;
    x = [base.(fac), lev(ok).(fac)]; y = [base.(met), lev(ok).(met)];
    good = x > 0 & y > 0 & isfinite(x) & isfinite(y);
    s = NaN;
    if sum(good) >= 3 && base.valid
        c = polyfit(log(x(good)), log(y(good)), 1); s = c(1);
    end
end


function V = verification_table(res, cfg)
%VERIFICATION_TABLE  per manipulation (arm A, arm X): the genotype-levels that fail the checks
%   and whether the arm is kept.
    arms = {'A', 'X'}; fac = {'alpha', 'xn'}; oth = {'xn', 'alpha'}; o0 = {'xn0', 'alpha0'};
    for q = 1:2
        n = 0; nFail = 0; tgt = []; othMove = [];
        for k = 1:numel(res)
            z = res{k}; lev = z.(arms{q});
            for j = 1:numel(lev)
                e = lev(j); n = n + 1;
                om = abs(e.(oth{q})/z.(o0{q}) - 1);
                pass = e.newtonOK && e.valid && e.fdErr <= cfg.fdTol && om <= cfg.otherTol;
                nFail = nFail + ~pass; othMove(end+1) = om; %#ok<AGROW>
                if e.m == min(cfg.levels)
                    tgt(end+1) = abs(log(e.(fac{q})/z.([fac{q} '0']))); %#ok<AGROW>
                end
            end
        end
        keep = nFail <= n/2 && median(tgt, 'omitnan') >= log(1/0.9) && ...
               median(othMove, 'omitnan') <= cfg.otherTol;
        V.(['keep' arms{q}]) = keep;
        V.(['fail' arms{q}]) = [nFail n];
        V.(['tgt' arms{q}]) = median(tgt, 'omitnan');
        V.(['oth' arms{q}]) = [median(othMove, 'omitnan') max(othMove)];
    end
end


function V = report_verification(res, G, cfg)
    V = verification_table(res, cfg);
    fprintf('\n---- MANIPULATION VERIFICATION ----\n');
    for q = {'A', 'X'}
        a = q{1};
        fprintf('arm %s: genotype-levels failing a verification check %d of %d | median |log target move| at m = %.1f: %.3f (needs >= %.3f) | other factor move median %.2e, max %.2e | KEPT: %s\n', ...
            a, V.(['fail' a])(1), V.(['fail' a])(2), min(cfg.levels), V.(['tgt' a]), log(1/0.9), ...
            V.(['oth' a])(1), V.(['oth' a])(2), ternary(V.(['keep' a]), 'yes', 'NO'));
    end
    fprintf('\nper level, median over genotypes [alpha ratio, ||x*|| ratio, alpha_cue ratio, b, delta, fdErr max, flip, n valid]\n');
    for q = {'A', 'X'}
        for j = 1:numel(cfg.levels)
            M = level_matrix(res, q{1}, j);
            fprintf('arm %s m = %.1f: %.4f %.4f %.4f  b = %.3f  delta = %.3f  fd %.1e  flip %.4f  valid %d/%d\n', ...
                q{1}, cfg.levels(j), median(M(:,1),'omitnan'), median(M(:,2),'omitnan'), ...
                median(M(:,3),'omitnan'), median(M(:,4),'omitnan'), median(M(:,5),'omitnan'), ...
                max(M(:,6)), median(M(:,7),'omitnan'), sum(M(:,8)), size(M,1));
        end
    end
    B = cellfun(@(z) z.base.flip, res); F = cellfun(@(z) z.base.fdErr, res);
    fprintf('baseline: flip median %.4f, fdErr max %.1e, hmin min %.2e\n', median(B), max(F), ...
        min(cellfun(@(z) z.base.hmin, res)));
end


function M = level_matrix(res, arm, j)
    M = nan(numel(res), 8);
    for k = 1:numel(res)
        z = res{k}; e = z.(arm)(j);
        ok = e.newtonOK && e.valid;
        M(k, :) = [e.alpha/z.alpha0, e.xn/z.xn0, e.alpha_cue/max(z.base.alpha_cue, eps), ...
                   e.b, e.delta, e.fdErr, e.flip, ok];
    end
end


function R = pooled(sA, sX, G, cfg, met, keepA, keepX)
%POOLED  mean slopes, paired bootstrap stratified by (target x design) cell.
    ok = isfinite(sA) & isfinite(sX); cellId = [G.ai]'*10 + [G.f]';
    idx = find(ok); cid = cellId(idx); cells = unique(cid);
    st = RandStream('twister', 'Seed', cfg.bootSeed);
    bA = nan(cfg.nBoot, 1); bX = bA;
    for bb = 1:cfg.nBoot
        pick = [];
        for c = cells'
            ic = idx(cid == c); pick = [pick; ic(randi(st, numel(ic), numel(ic), 1))]; %#ok<AGROW>
        end
        bA(bb) = mean(sA(pick)); bX(bb) = mean(sX(pick));
    end
    R = struct('met', met, 'n', numel(idx), 'sA', mean(sA(idx)), 'sX', mean(sX(idx)), ...
        'ciA', prctile(bA, [2.5 97.5]), 'ciX', prctile(bX, [2.5 97.5]));
    R.ratio = R.sA / R.sX; R.ciRatio = prctile(bA ./ bX, [2.5 97.5]);
    fprintf('[%s] n = %d genotypes | s_A (log %s vs log alpha) = %.4f [%.4f, %.4f] | s_X (vs log ||x*||) = %.4f [%.4f, %.4f] | ratio s_A/s_X = %.3f [%.3f, %.3f] | arms kept A %d X %d\n', ...
        met, R.n, met, R.sA, R.ciA(1), R.ciA(2), R.sX, R.ciX(1), R.ciX(2), ...
        R.ratio, R.ciRatio(1), R.ciRatio(2), keepA, keepX);
end


function C = report_cells(sA, sX, G)
    fprintf('\n---- per cell (not decisive): median [IQR] of per-genotype slopes ----\n');
    arms = {'moving', 'fixed'}; C = struct([]);
    for ai = 1:2
        for f = unique([G.f])
            k = [G.ai] == ai & [G.f] == f; a = sA(k); x = sX(k);
            fprintf('%-7s a = %.2f: s_A %.3f [%.3f, %.3f] | s_X %.3f [%.3f, %.3f] | n %d/%d\n', ...
                arms{ai}, G(find(k,1)).a, median(a,'omitnan'), prctile(a,25), prctile(a,75), ...
                median(x,'omitnan'), prctile(x,25), prctile(x,75), sum(isfinite(a) & isfinite(x)), sum(k));
            C = [C, struct('arm', ai, 'f', f, 'sA', a, 'sX', x)]; %#ok<AGROW>
        end
    end
end


function Gr = report_gain(res, G, cfg)
%REPORT_GAIN  arm G: whether the input gain moves alpha, and the slope of log L0 on the log of
%   the cue channel's responsiveness alpha_cue.
    fprintf('\n---- arm G (input gain, target moving; reported) ----\n');
    k = find([G.ai] == 1); nL = numel(cfg.gains);
    ra = nan(numel(k), nL); rx = ra; rc = ra; s = nan(numel(k), 1);
    for q = 1:numel(k)
        z = res{k(q)}; lev = z.G;
        ra(q,:) = [lev.alpha]/z.alpha0; rx(q,:) = [lev.xn]/z.xn0;
        rc(q,:) = [lev.alpha_cue]/max(z.base.alpha_cue, eps);
        x = [z.base.alpha_cue, lev.alpha_cue]; y = [z.base.L0, lev.L0];
        ok = x > 0 & y > 0 & isfinite(y);
        if sum(ok) >= 3, c = polyfit(log(x(ok)), log(y(ok)), 1); s(q) = c(1); end
    end
    for j = 1:nL
        fprintf('g = %.2f: alpha ratio median %.4f | ||x*|| ratio %.4f | alpha_cue ratio %.4f\n', ...
            cfg.gains(j), median(ra(:,j),'omitnan'), median(rx(:,j),'omitnan'), median(rc(:,j),'omitnan'));
    end
    ext = [1 nL]; mvA = median(abs(ra(:, ext) - 1), 'all', 'omitnan');
    mvX = median(abs(rx(:, ext) - 1), 'all', 'omitnan');
    Gr = struct('alphaMove', mvA, 'xnMove', mvX, 'slopeCue', s, 'keptAsAlpha', mvA >= 0.10 && mvX <= 0.20);
    fprintf('median |alpha move| at the extreme gains %.3f (>= 0.10 to count as an alpha manipulation), ||x*|| %.3f -> %s\n', ...
        mvA, mvX, ternary(Gr.keptAsAlpha, 'ALPHA CANDIDATE PASSES', 'does NOT move alpha: reported as a cue-responsiveness arm'));
    fprintf('slope of log L0 on log alpha_cue (g*||C Phi(:,1:2)||): median %.4f [IQR %.4f, %.4f], n %d\n', ...
        median(s,'omitnan'), prctile(s,25), prctile(s,75), sum(isfinite(s)));
end


% Utilities
function s = cat_struct(s, e)
    if isempty(s), s = e; else, s(end+1) = e; end
end

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
