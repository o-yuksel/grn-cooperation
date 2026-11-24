function out = run_tracking_task(cfg)
% RUN_TRACKING_TASK  Evolves the tracking task in its four designs, with the target moving and
% held still, measures the mean shift per mutation of each run's best genotype, and stores the
% populations and genotypes that every later tracking script reads.
%
% Produces : the mean shift per mutation of the generation-350 best genotypes, Table S13 (printed
%            as 'beta', mean per design and target); the stored runs behind Fig. 6, SI S5 and
%            Methods, "Non-social tracking task" (read by tracking_responsiveness,
%            tracking_load_vs_identity, tracking_load_excess, run_tracking_single_genotypes,
%            run_tracking_cut_vs_shift, run_tracking_moving_target and run_tracking_held_still)
% Inputs   : cfg (optional): nPop (1000), nReps (30), aList (each design's output
%              self-activation a, [0.70 0.28 0.20 0.15]), nGen (350), mutSd (0.020), betaDraws
%              (300), nWorkers (5); the rest as set in the code
% Outputs  : publication/tracking_runs/tracking_runs_N<nPop>.mat, per target
%              ('moving'; 'fixed', held still): readouts, final best genotype ThBest, generation-1
%              best ThG1, distinct final wirings popTh and their shares popFreq, scoring window;
%              tracking_runs_N<nPop>_run.txt (the console text)
% Run      : run_tracking_task()   (about 1.5 h on 5 workers)
% Seed     : replicate r evolves from rng(r, 'twister') under either target (shared founders); its
%            mean shift per mutation uses rng(20260917 + 1000 f + r, 'twister'), at snapshot k
%            20260917 + 7919 k + 1000 f + r; f is the index in aList (design 1 is f = 1 and design 2 is f = 3)
%
% Genes 1 and 2 receive the target's signal; genes 3 and 4 are the outputs, and a design fixes
% their self-activation W33 = W44 = a, never mutated. The mean shift per mutation is the mean
% distance the outputs move under one application of the mutation distribution (each other entry
% mutated with probability 0.08 by N(0, mutSd^2)), over betaDraws draws, at the probe input
% u = [0, 0.25, 0, 0]; betaOp is the same at the target's final position.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'nPop'),     cfg.nPop   = 1000; end     % population size N
    if ~isfield(cfg,'nReps'),    cfg.nReps  = 30;   end     % replicates per design and arm
    if ~isfield(cfg,'aList'),    cfg.aList  = [0.70 0.28 0.20 0.15]; end
    if ~isfield(cfg,'nGen'),     cfg.nGen   = 350;  end
    if ~isfield(cfg,'mutSd'),    cfg.mutSd  = 0.020; end    % sd of a mutation (SI S5, Specification)
    if ~isfield(cfg,'horizon'),  cfg.horizon = 60;  end
    if ~isfield(cfg,'snapGens'), cfg.snapGens = [1 25 50 100 150 200 250 300 350]; end
    if ~isfield(cfg,'nSample'),  cfg.nSample = 8;   end
    if ~isfield(cfg,'betaDraws'),cfg.betaDraws = 300; end   % draws per mean shift per mutation
    if ~isfield(cfg,'arms'),     cfg.arms = {'moving','fixed'}; end
    if ~isfield(cfg,'rngSeed'),  cfg.rngSeed = 20260917; end
    if ~isfield(cfg,'nWorkers'), cfg.nWorkers = 5; end
    dd = fullfile(grn_data_file('publication'),'tracking_runs');
    if ~isfield(cfg,'outFile'),  cfg.outFile = fullfile(dd, sprintf('tracking_runs_N%d.mat', cfg.nPop)); end
    if ~isfield(cfg,'txtFile'),  cfg.txtFile = fullfile(dd, sprintf('tracking_runs_N%d_run.txt', cfg.nPop)); end

    if ~exist(dd,'dir'), mkdir(dd); end
    diary_file(cfg.txtFile, true);

    P = tracking_config();
    P.mutation_sd   = cfg.mutSd;
    P.n_pop         = cfg.nPop;
    P.n_generations = cfg.nGen;
    P.horizon       = cfg.horizon;
    snapGens = cfg.snapGens(cfg.snapGens <= P.n_generations);

    % Angular step of the target per generation: 2 pi / 36 in arm 'moving', 0 in arm 'fixed'
    % (the stored arm name for the target held still)
    omegaOf = struct('moving', P.omega_fast, 'fixed', 0.0);
    uProbe = [0; P.trajectory_radius; 0; 0];

    fprintf('\n=== TRACKING TASK RUN, N = %d, %d replicates ===\n', cfg.nPop, cfg.nReps);
    fprintf('a = %s | %d gens | mutation %.3f at p = %.2f | horizon %d | snap %s\n', ...
        mat2str(cfg.aList), P.n_generations, P.mutation_sd, P.mutation_prob, P.horizon, mat2str(snapGens));
    fprintf('arms: %s (omega_fast = %.5f; arm fixed: target held still, omega = 0) | beta at %d draws\n', ...
        strjoin(cfg.arms, ', '), P.omega_fast, cfg.betaDraws);
    fprintf('mutation regime: per-entry p * 16 entries * N = %.1f; expected mutations per generation (14 mutable entries): %.1f\n', ...
        16*P.mutation_prob*cfg.nPop, 14*P.mutation_prob*cfg.nPop);
    start_pool(cfg.nWorkers);

    nA = numel(cfg.aList); nS = numel(snapGens); nR = cfg.nReps;
    out = struct('cfg',cfg,'snapGens',snapGens,'aList',cfg.aList,'P',P,'nPop',cfg.nPop,'nReps',nR);

    % Founder input gain w0 of each design: the value with the best generation-1 tracking
    w0s = zeros(1, nA);
    for f = 1:nA, w0s(f) = tracking.pick_w0(cfg.aList(f), P); end
    out.w0 = w0s;
    fprintf('pick_w0 = %s\n', mat2str(w0s, 4));

    tAll = tic;
    for ai_arm = 1:numel(cfg.arms)
        arm = cfg.arms{ai_arm};
        omega = omegaOf.(arm);
        R = struct('omega',omega,'a',cfg.aList,'n_reps',nR,'w0',w0s);
        R.maint = zeros(nA,nR); R.xnorm = zeros(nA,nR);
        R.beta  = zeros(nA,nR); R.betaOp = zeros(nA,nR);
        R.gen_track = zeros(nA,nR); R.plas_track = zeros(nA,nR);
        R.betaBest = nan(nA,nR,nS); R.betaMean = nan(nA,nR,nS); R.fitSnap = nan(nA,nR,nS);
        R.fit_traj = zeros(nA,P.n_generations); R.xnorm_traj = zeros(nA,P.n_generations);
        R.fit_rep = nan(nA,nR,P.n_generations);
        R.ThBest = cell(nA,nR); R.ThG1 = cell(nA,nR);
        R.popTh = cell(nA,nR); R.popFreq = cell(nA,nR);
        R.nDistinct = nan(nA,nR); R.uEnd = cell(nA,1); R.cc = cell(nA,1);
        R.winU = cell(nA,1); R.winO = cell(nA,1); R.winW = cell(nA,1);

        jobs = cell(nA*nR,1); jf = zeros(nA*nR,1); jr = jf; q = 0;
        for f = 1:nA
            for r = 1:nR, q = q+1; jf(q) = f; jr(q) = r; end
        end
        aList = cfg.aList; rngSeed = cfg.rngSeed;
        nSample = cfg.nSample; betaDraws = cfg.betaDraws;
        t0 = tic;
        parfor q = 1:numel(jf)
            jobs{q} = one_replicate(aList(jf(q)), w0s(jf(q)), omega, P, jr(q), ...
                snapGens, uProbe, rngSeed, jf(q), nSample, betaDraws);
        end
        tArm = toc(t0);
        fprintf('[%s] %d replicates in %.1f min (%.1f s per replicate on %d workers)\n', ...
            arm, numel(jf), tArm/60, tArm/numel(jf), cfg.nWorkers);

        accF = zeros(nA,P.n_generations); accX = zeros(nA,P.n_generations);
        for q = 1:numel(jf)
            f = jf(q); r = jr(q); z = jobs{q};
            R.maint(f,r) = z.maint; R.xnorm(f,r) = z.xnorm;
            R.gen_track(f,r) = z.gen_track; R.plas_track(f,r) = z.plas_track;
            R.beta(f,r) = z.beta; R.betaOp(f,r) = z.betaOp;
            R.betaBest(f,r,:) = z.betaBest; R.betaMean(f,r,:) = z.betaMean;
            R.fitSnap(f,r,:) = z.fitSnap;
            R.fit_rep(f,r,:) = z.fit;
            R.ThBest{f,r} = z.ThBest; R.ThG1{f,r} = z.ThG1;
            R.popTh{f,r} = z.popTh; R.popFreq{f,r} = z.popFreq;
            R.nDistinct(f,r) = z.nDistinct;
            accF(f,:) = accF(f,:) + z.fit; accX(f,:) = accX(f,:) + z.xnorm_traj;
            % The final input, the founder's output and the scoring window depend only on the
            % design and the arm, so replicate 1's serve every replicate
            if r == 1
                R.uEnd{f} = z.uEnd; R.cc{f} = z.cc;
                R.winU{f} = z.winU; R.winO{f} = z.winO; R.winW{f} = z.winW;
            end
        end
        R.fit_traj = accF/nR; R.xnorm_traj = accX/nR;
        R.elapsedMin = tArm/60;

        for f = 1:nA
            fprintf('[%s] a=%.2f | maint=%.3f+-%.3f | beta=%.4f | betaOp=%.4f | ||x*||=%.2f | gen=%.2f plas=%.2f | distinct %.0f/%d\n', ...
                arm, cfg.aList(f), mean(R.maint(f,:)), std(R.maint(f,:)), mean(R.beta(f,:)), ...
                mean(R.betaOp(f,:)), mean(R.xnorm(f,:)), mean(R.gen_track(f,:)), ...
                mean(R.plas_track(f,:)), median(R.nDistinct(f,:)), cfg.nPop);
        end
        out.(arm) = R;
        save(cfg.outFile, 'out', '-v7.3');    % incremental: the first arm survives a crash in the second
        fprintf('[%s] saved.\n', arm);
    end

    out.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    out.elapsedHours = toc(tAll)/3600;
    save(cfg.outFile, 'out', '-v7.3');
    fprintf('\nSaved %s  (%.2f h)\n', cfg.outFile, out.elapsedHours);
    diary_file(cfg.txtFile, false);
end


% =====================================================================
function z = one_replicate(a, w0, omega, P, r, snapGens, uProbe, rngSeed, f, nSample, betaDraws)
%ONE_REPLICATE  One evolutionary run plus every readout, self-seeded so parfor order is irrelevant.
%   `beta` and `betaOp` are the mean shifts per mutation of the final best genotype at the probe input and
%   at the arm's final input; at each snapshot k, the population's mean fitness, the shift of its
%   best individual and the mean shift of nSample individuals drawn at random.
    H = evolve_snap_fast(a, w0, omega, P, r, snapGens);
    o = track_readouts_fast(H, P);
    z = struct();
    z.maint = o.maint; z.xnorm = o.xnorm;
    z.gen_track = o.gen_track; z.plas_track = o.plas_track;
    z.fit = H.fit; z.xnorm_traj = H.xnorm;
    z.uEnd = H.uEnd; z.cc = H.cc;

    rng(rngSeed + 1000*f + r, 'twister');
    z.beta   = fragility_fast(H.pop{o.best}, uProbe, P, betaDraws);
    z.betaOp = fragility_fast(H.pop{o.best}, H.uEnd, P, betaDraws);

    nS = numel(snapGens);
    z.betaBest = nan(1,nS); z.betaMean = nan(1,nS); z.fitSnap = nan(1,nS);
    for k = 1:nS
        pk = H.snapPop{k}; fk = H.snapFit{k};
        [~, ib] = max(fk);
        z.fitSnap(k) = mean(fk);
        rng(rngSeed + 7919*k + 1000*f + r, 'twister');
        z.betaBest(k) = fragility_fast(pk{ib}, uProbe, P, betaDraws);
        idx = randperm(numel(pk), min(nSample, numel(pk)));
        bs = zeros(1, numel(idx));
        for q = 1:numel(idx)
            bs(q) = fragility_fast(pk{idx(q)}, uProbe, P, betaDraws);
        end
        z.betaMean(k) = mean(bs);
    end

    % ---- the genotypes the other tracking files read: the final best, the generation-1 best,
    % and the distinct final wirings (Mall: one row per individual, its 4 x 4 wiring in column
    % order; unique rows) with the share of individuals carrying each ----
    z.ThBest = H.pop{o.best};
    p1 = H.snapPop{1}; f1 = H.snapFit{1}; [~, b1] = max(f1);
    z.ThG1 = p1{b1};
    np = numel(H.pop);
    Mall = reshape(cat(3, H.pop{:}), P.total*P.total, np)';
    [Uw, ~, iw] = unique(Mall, 'rows', 'stable');
    fr = accumarray(iw,1); fr = fr(:)'/np;
    z.popTh = Uw; z.popFreq = fr; z.nDistinct = size(Uw,1);

    % ---- the scoring window: the last min(36, horizon) generations, each with its input and
    % optimum; a repeated point is kept once, weighted by how often it occurs ----
    G = P.n_generations; n = min(36, P.horizon); gi = (G-n+1):G;
    pos = H.opt(:,gi) - H.cc;
    U = zeros(P.total, n); U(1:2,:) = pos; O = H.opt(:,gi);
    [~, iu, ic] = unique([U; O]', 'rows', 'stable');
    z.winU = U(:,iu); z.winO = O(:,iu);
    wt = accumarray(ic,1); z.winW = wt(:)'/sum(wt);
end


% =====================================================================
function H = evolve_snap_fast(a, w0, omega, P, seed, snapGens)
%EVOLVE_SNAP_FAST  One evolutionary run of the tracking task (SI S5, Specification), the settled
%   states of the whole population solved together by steady_batch. Founders: the design's
%   founder wiring plus N(0, 0.01^2) on every entry, W33 = W44 reset to a. Each generation:
%   fitness exp(-gamma d^2) against the optimum at the target's position, N offspring drawn with
%   replacement in proportion to fitness (uniformly if every fitness is zero), then each entry
%   mutated with probability P.mutation_prob by N(0, P.mutation_sd^2) and W33 = W44 reset to a.
    rng(seed, 'twister');
    tot = P.total; N = P.n_pop;
    Th0 = tracking.build_theta(a, w0, P);
    cc  = tracking.steady(Th0, zeros(tot,1), P); cc = cc([3 4]);
    [u_path, ~, opt] = tracking.circle(omega, cc, P.n_generations, P);

    M = zeros(tot, tot, N);
    for i = 1:N
        T = Th0 + randn(tot) * 0.01;
        T(3,3) = a; T(4,4) = a;
        M(:,:,i) = T;
    end

    G = P.n_generations;
    H.fit = zeros(1,G); H.xnorm = zeros(1,G);
    H.exp = zeros(2,G); H.ref = zeros(2,G);
    H.snapPop = cell(1, numel(snapGens)); H.snapFit = cell(1, numel(snapGens));
    uref = zeros(tot,1);

    for gen = 1:G
        % Scored at the target's input u (Xe) and, for the tracking readouts, without input (Xr)
        u = u_path(:,gen); o = opt(:,gen);
        Xe = steady_batch(M, u,    P);
        Xr = steady_batch(M, uref, P);
        pe = Xe([3 4],:); pr = Xr([3 4],:);
        fit = exp(-P.selection_strength * sum((pe - o).^2, 1))';
        xn  = vecnorm(Xe, 2, 1)';

        H.fit(gen) = mean(fit); H.xnorm(gen) = mean(xn);
        H.exp(:,gen) = mean(pe,2); H.ref(:,gen) = mean(pr,2);

        k = find(snapGens == gen, 1);
        if ~isempty(k)
            H.snapPop{k} = squeeze(num2cell(M, [1 2]));
            H.snapFit{k} = fit;
        end

        if any(fit > 0); wsel = fit / sum(fit); else; wsel = ones(N,1)/N; end
        ix = randsample(1:N, N, true, wsel);
        M = M(:,:,ix);
        mask = rand(size(M)) < P.mutation_prob;
        M(mask) = M(mask) + randn(nnz(mask),1) * P.mutation_sd;
        M(3,3,:) = a; M(4,4,:) = a;
    end
    H.pop = squeeze(num2cell(M, [1 2]));
    H.cc = cc; H.opt = opt; H.uEnd = u_path(:,G);
end


% =====================================================================
function out = track_readouts_fast(H, P)
%TRACK_READOUTS_FAST  Readouts over the last P.horizon generations: maint, the population's mean
%   fitness; xnorm, its mean ||x*||; gen_track, one minus the mean distance of the output without
%   input from the optimum, over the radius; plas_track, how much the input brings the output
%   closer, on the same scale. best: the individual with the highest mean fitness over the
%   scoring window (the last min(36, horizon) generations).
    G = P.n_generations; idx = (G - P.horizon + 1):G; R = P.trajectory_radius;
    out.maint = mean(H.fit(idx));
    out.xnorm = mean(H.xnorm(idx));
    gerr = mean( sqrt(sum((H.ref(:,idx) - H.opt(:,idx)).^2, 1)) ) / R;
    eerr = mean( sqrt(sum((H.exp(:,idx) - H.opt(:,idx)).^2, 1)) ) / R;
    out.gen_track  = max(0, 1 - gerr);
    out.plas_track = max(0, gerr - eerr);

    n = min(36, P.horizon); gi = (G - n + 1):G;
    pos = H.opt(:,gi) - H.cc; u = zeros(P.total, n); u(1:2,:) = pos;
    np = numel(H.pop);
    M = reshape(cat(3, H.pop{:}), P.total, P.total, np);
    d2 = zeros(np, n);
    for g = 1:n
        X = steady_batch(M, u(:,g), P);
        d2(:,g) = sum((X([3 4],:) - H.opt(:,gi(g))).^2, 1)';
    end
    f = mean(exp(-P.selection_strength * d2), 2);
    best = -inf; out.best = 1;
    for i = 1:np
        if f(i) > best; best = f(i); out.best = i; end
    end
end


% =====================================================================
function start_pool(n)
    try
        pl = gcp('nocreate');
        if ~isempty(pl) && pl.NumWorkers == n, return; end
        if ~isempty(pl), delete(pl); end
        parpool('Processes', n);
    catch ME
        fprintf('parpool failed (%s) -- running serial\n', ME.message);
    end
end

function diary_file(f, on)
    try
        if on
            d = fileparts(f); if ~exist(d,'dir'), mkdir(d); end
            if exist(f,'file'), delete(f); end
            diary(f); diary on;
        else
            diary off;
        end
    catch
    end
end
