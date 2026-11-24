function out = membership_test_vs_encounter_rule(cfg)
% MEMBERSHIP_TEST_VS_ENCOUNTER_RULE  How often do the membership test and the runs' encounter rule
% agree on the cooperation condition (C1), and how often does a member's encounter with its clone
% have more than one outcome under the runs' rule? Read on generation-500 wirings of nine settings:
% the four signaling architectures and the five spatial settings (populations of 50, well mixed
% or on rings).
%
% Produces : SI S4.5: agreement on (C1) and Cohen's kappa (printed block [1]); member pairings
%            with the clone that have more than one outcome (block [3])
% Inputs   : cfg (optional struct; defaults in fill_defaults: K 20, nNeigh 20, nScreen 2,
%              nBoot 10000, nWorkers 2 (at most 2), nGen1 10000, stage 'arch')
%            publication/<architecture>_100reps.mat   (run_signaling_architectures)
%            publication/sweep_spatial_N50_<setting>_100reps_*.mat   (run_spatial_n50; the
%              newest matching file)
% Outputs  : publication/membership_vs_encounter_9settings/membership_vs_encounter_<stage>.mat (one per stage)
%            and membership_vs_encounter_agreement.txt there (the report, from the stage files present)
% Run      : membership_test_vs_encounter_rule(struct('stage', s)) for s = 'arch' (the four
%            architectures), 'struct' (the spatial settings), 'gen1' (random generation-1 wirings,
%            block [5], not read in SI S4.5) (not timed); then struct('stage', 'report') (about
%            8 min).
% Seed     : sample of 30 individuals per signaling replicate under rng(42) (the session's
%            generator); job j of a stage under rng(20260927 + offset + j, 'twister'), offset 0,
%            1e6, 2e6 for 'arch', 'struct', 'gen1'; other seeds beside their blocks; membership_test
%            seeds each wiring itself
%
% Runs' rule here: each player starts from U(0, 0.01), 300 steps; a wiring meets (C1) when it
% cooperates in at least half of the 2K actions of K = 20 encounters with its clone.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
    if ~strcmp(cfg.stage, 'report')
        st = run_stage(cfg);
        save(fullfile(cfg.outDir, ['membership_vs_encounter_' cfg.stage '.mat']), 'st', '-v7');
        out = []; return;
    end
    S = {};
    for s = {'arch', 'struct', 'gen1'}
        f = fullfile(cfg.outDir, ['membership_vs_encounter_' s{1} '.mat']);
        if isfile(f), L = load(f, 'st'); S{end+1} = L.st; end %#ok<AGROW>
    end
    out = assemble_stages(S, cfg);
    txt = fullfile(cfg.outDir, 'membership_vs_encounter_agreement.txt');
    if exist(txt, 'file'), delete(txt); end
    diary(txt); print_report(out); diary off;
end


% =========================================================================
function cfg = fill_defaults(cfg)
    pub = fullfile(grn_data_file('publication'));
    def = struct('dataDir', pub, 'outDir', fullfile(pub, 'membership_vs_encounter_9settings'), ...
        'K', 20, 'nNeigh', 20, 'nScreen', 2, 'rngSeed', 20260927, 'bootSeed', 20260928, ...
        'nBoot', 10000, 'nWorkers', 2, 'nGen1', 10000, 'stage', 'arch');
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    assert(cfg.nWorkers <= 2, 'use at most 2 workers (cfg.nWorkers)');
end


% =========================================================================
function st = run_stage(cfg)
% One stage: build its jobs, run them on <= 2 workers, return the checkpoint.
    t0 = tic; C = setting_table();
    switch cfg.stage
        case 'arch',   cond = 1:4;  off = 0;
        case 'struct', cond = 5:9;  off = 1e6;
        case 'gen1',   cond = [];   off = 2e6;
    end
    [samp, jobs, xjobs] = draw_samples(C, cond, cfg, off);
    g1 = struct([]);
    if strcmp(cfg.stage, 'gen1'), [g1, jobs] = gen1_jobs(cfg, jobs, off); end
    fprintf('stage %s: %d distinct-wiring jobs, %d cross-check jobs\n', cfg.stage, numel(jobs), numel(xjobs));
    p = gcp('nocreate');
    if ~isempty(p) && p.NumWorkers > cfg.nWorkers, delete(p); p = []; end
    if isempty(p), parpool('Processes', cfg.nWorkers); end
    P = struct('K', cfg.K, 'nNeigh', cfg.nNeigh, 'nScreen', cfg.nScreen); D = grid_deltas();
    res = cell(numel(jobs), 1);
    parfor j = 1:numel(jobs)
        res{j} = analyse_wiring(jobs{j}, P, D);
    end
    fprintf('wiring analysis %.1f min\n', toc(t0) / 60);
    xres = cell(numel(xjobs), 1);
    parfor j = 1:numel(xjobs)
        xres{j} = crosscheck_rep(xjobs{j}, P.K);
    end
    fprintf('stage %s total %.1f min\n', cfg.stage, toc(t0) / 60);
    st = struct('stage', cfg.stage, 'cond', cond, 'samp', samp, 'g1', g1, 'res', {res}, ...
        'xres', {xres}, 'cfg', cfg, 'runMin', toc(t0) / 60, 'date', ...
        char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm')));
end


% =========================================================================
function out = assemble_stages(S, cfg)
% Concatenates stage checkpoints; job indices are offset into one res list.
    C = setting_table(); out = struct('C', C([]), 'samp', struct([]), 'g1', struct([]), ...
        'res', {{}}, 'xres', {{}}, 'K', cfg.K, 'cfg', cfg, 'runMin', 0, 'date', '', 'stages', '');
    for k = 1:numel(S)
        st = S{k}; o = numel(out.res);
        for i = 1:numel(st.samp), st.samp(i).firstJob = st.samp(i).firstJob + o; end
        for a = 1:numel(st.g1), st.g1(a).firstJob = st.g1(a).firstJob + o; end
        out.C = [out.C, C(st.cond)];
        if isempty(out.samp), out.samp = st.samp; else, out.samp = [out.samp, st.samp]; end
        if isempty(out.g1), out.g1 = st.g1(:); elseif ~isempty(st.g1), out.g1 = [out.g1; st.g1(:)]; end
        out.res = [out.res; st.res]; out.xres = [out.xres; st.xres];
        out.runMin = out.runMin + st.runMin; out.date = st.date;
        out.stages = strtrim([out.stages ' ' st.stage]); out.K = st.cfg.K; out.cfg = st.cfg;
    end
end


% =========================================================================
function C = setting_table()
% The nine signaling and spatial settings, in the order of cooperation_condition_shares_gen500. nsub:
% individuals sampled per replicate (Inf = all).
    add = @(f,l,p,ns,a) struct('family',f,'label',l,'pat',p,'nsub',ns,'arch',a);
    C = add('signaling','graded separate','graded_separate_100reps.mat',30,1);
    C(2)  = add('signaling','graded contrast','graded_contrast_100reps.mat',30,2);
    C(3)  = add('signaling','saturating separate','saturating_separate_100reps.mat',30,3);
    C(4)  = add('signaling','saturating contrast','saturating_contrast_100reps.mat',30,4);
    C(5)  = add('spatial','well-mixed','sweep_spatial_N50_wellmixed_100reps_*.mat',Inf,0);
    C(6)  = add('spatial','ring k=1','sweep_spatial_N50_ring_k1_100reps_*.mat',Inf,0);
    C(7)  = add('spatial','ring k=2','sweep_spatial_N50_ring_k2_100reps_*.mat',Inf,0);
    C(8)  = add('spatial','ring k=5','sweep_spatial_N50_ring_k5_100reps_*.mat',Inf,0);
    C(9)  = add('spatial','ring k=10','sweep_spatial_N50_ring_k10_100reps_*.mat',Inf,0);
end


function s = success_pct(A)
% Success in % per setting: the share of the sampled replicates with final cooperation >= 0.95.
    s = arrayfun(@(a) 100 * mean(a.succ), A)';
end


% =========================================================================
function [samp, jobs, xjobs] = draw_samples(C, cond, cfg, off)
% Draws the sample of the settings `cond` before any membership test (signaling: rng(42) once
% before the first setting, 30 individuals per replicate; spatial: every individual, no draw)
% and makes one job per distinct wiring. One setting in memory at a time. For the four architectures the
% cross-check jobs (control (c)) are cut from the same load.
    jobs = {}; xjobs = {}; samp = struct('rep', {}, 'ic', {}, 'nRep', {}, 'succ', {}, 'firstJob', {});
    rng(42);
    for i = cond
        [reps, sp] = load_setting(cfg.dataDir, C(i));
        n = numel(reps);
        pick = cell(numel(reps), 1);
        for k = 1:numel(reps)
            nAg = size(reps{k}.finalNetworkB, 3);
            if isinf(C(i).nsub), pick{k} = 1:nAg; else, pick{k} = randperm(nAg, min(C(i).nsub, nAg)); end
        end
        Wl = []; rl = [];
        for k = 1:n
            M = double(reps{k}.finalNetworkM(:,:,pick{k})); B = double(reps{k}.finalNetworkB(:,:,pick{k}));
            Wl = [Wl; reshape(M .* B, 16, [])']; rl = [rl; k * ones(numel(pick{k}), 1)]; %#ok<AGROW>
        end
        [Wu, ~, ic] = unique(Wl, 'rows');
        s = struct('rep', rl, 'ic', ic, 'nRep', n, 'firstJob', numel(jobs) + 1, ...
            'succ', cellfun(@(r) double(r.historyCooperation(end)) >= 0.95, reps(1:n)));
        samp(end+1) = s; %#ok<AGROW>
        for u = 1:size(Wu, 1)
            jobs{end+1} = make_job(reshape(Wu(u,:), 4, 4), sp, cfg, off + numel(jobs) + 1, i, u); %#ok<AGROW>
        end
        if C(i).arch > 0, xjobs = [xjobs, crosscheck_jobs(reps(1:n), sp, C(i).arch, cfg)]; end %#ok<AGROW>
        fprintf('%-20s %d sampled, %d distinct\n', C(i).label, numel(rl), size(Wu, 1));
        clear reps
    end
end


% =========================================================================
function [reps, sp] = load_setting(dataDir, c)
    f = dir(fullfile(dataDir, c.pat));
    assert(~isempty(f), 'No file matching %s', c.pat);
    [~, mx] = max([f.datenum]);
    d = load(fullfile(dataDir, f(mx).name));
    if isfield(d, 'results'), reps = d.results; cf = d.config; else, reps = d.res; cf = d.cfg; end
    sp = sim_params(cf.sensingMode, cf.hillCoefficient, cf.discreteOutput, cf.michaelisConstant, ...
        cf.degradationRate, cf.timeStep, cf.sensorIntegrationRate);
end


% =========================================================================
function sp = sim_params(sm, hc, disc, km, mu, dt, ir)
    sp = struct('sensingMode', sm, 'hillCoefficient', hc, 'discreteOutput', disc, ...
        'michaelisConstant', km, 'degradationRate', mu, 'timeStep', dt, ...
        'developmentSteps', 600, 'payoffDC', 0, 'sensorIntegrationRate', ir);
end


% =========================================================================
function job = make_job(W, sp, cfg, j, cond, u)
% cp = the membership test's parameters (membership_test, simulate_pair, 600 steps).
% j already carries the stage offset.
    cp = struct('sensingMode', sp.sensingMode, 'hillCoefficient', sp.hillCoefficient, ...
        'discreteOutput', sp.discreteOutput, 'michaelisConstant', sp.michaelisConstant, ...
        'degradationRate', sp.degradationRate, 'timeStep', sp.timeStep, ...
        'developmentSteps', 600, 'backend', 'simulate_pair');
    job = struct('W', W, 'sp', sp, 'cp', cp, 'seed', cfg.rngSeed + j, 'cond', cond, 'u', u);
end


% =========================================================================
function [g1, jobs] = gen1_jobs(cfg, jobs, off)
% cooperation_condition_shares.m's generation-1 sample: rng(2026), randn(4) * 0.1, the same
% wirings for every architecture (membership_test restores the caller's stream).
    A = {'graded separate', 'separate', Inf; 'graded contrast', 'contrast', Inf;
         'saturating separate', 'separate', 1; 'saturating contrast', 'contrast', 1};
    rng(2026); Wg = zeros(4, 4, cfg.nGen1);
    for i = 1:cfg.nGen1, Wg(:,:,i) = randn(4) * 0.1; end
    g1 = struct('name', A(:,1), 'firstJob', 0, 'n', cfg.nGen1);
    for a = 1:4
        sp = sim_params(A{a,2}, A{a,3}, false, 1.0, 0.1, 1.0, 1.0);
        g1(a).firstJob = numel(jobs) + 1;
        for i = 1:cfg.nGen1
            jobs{end+1} = make_job(Wg(:,:,i), sp, cfg, off + numel(jobs) + 1, -a, i); %#ok<AGROW>
        end
    end
end


% =========================================================================
function D = grid_deltas()
% The membership test's grid: linspace(-3*0.2, 3*0.2, 11) (sigma_mut = 0.2) with the exact
% zero skipped, 16 x 10 = 160 single-weight mutants.
    d = linspace(-0.6, 0.6, 11); d = d(d ~= 0);
    [dd, jj, ii] = ndgrid(d, 1:4, 1:4);     % membership_test loop order: ii, jj, delta
    D = [ii(:), jj(:), dd(:)];
end


% =========================================================================
function r = analyse_wiring(job, P, D)
% Everything measured on one distinct wiring: the membership test (membership_test) and its batch
% emulation (control (b)); K encounters with its clone from random starts at 300 steps (the
% runs' rule) and at 600 steps, and one from the fixed start at 300 steps, all sharing the K
% start pairs. (C2) under the runs' rule is tested on up to P.nNeigh grid mutants drawn from
% the pool of those that defect against the wiring under the membership test's dynamics or in
% one of P.nScreen screening encounters.
    W = job.W; sp = job.sp; K = P.K;
    cs = membership_test(W, job.cp);
    r.is_weak = cs.is_weak; r.is_css = cs.is_css;
    X0 = 0.005 * ones(4, 1);
    F = grid_pairs(W, sp, 600, X0, X0, D, 'simpair');           % membership test, emulated
    r.em_c1 = F.c1 && F.c2;
    r.em_css = r.em_c1 && ~any(F.aR & ~F.aM);
    r.clsDefect = ~F.aM;                                         % grid mutant defects (membership test)
    rng(job.seed, 'twister');
    X1 = rand(4, K) * 0.01; X2 = rand(4, K) * 0.01;
    [c3a, c3b] = encounter_batch(repmat(W, 1, 1, K), repmat(W, 1, 1, K), X1, X2, sp, 300, 'encounter_rule');
    [c6a, c6b] = encounter_batch(repmat(W, 1, 1, K), repmat(W, 1, 1, K), X1, X2, sp, 600, 'encounter_rule');
    [f3a, f3b] = encounter_batch(W, W, X0, X0, sp, 300, 'encounter_rule');
    r.pC1 = mean([c3a c3b]); r.pCC = mean(c3a & c3b); r.multi300 = numel(unique(2*c3a + c3b)) > 1;
    r.pC1_600 = mean([c6a c6b]); r.multi600 = numel(unique(2*c6a + c6b)) > 1;
    r.c1_f300 = f3a && f3b;
    r.g300 = []; r.g600 = []; r.gf300 = []; r.nbr = []; r.pool = 0;
    if r.pC1 >= 0.5 || r.is_weak
        S = grid_pairs(W, sp, 300, X1(:, 1:P.nScreen), X2(:, 1:P.nScreen), D, 'encounter_rule');
        pool = find(r.clsDefect | S.pDM > 0);
        r.pool = numel(pool);
        r.nbr = sort(pool(randperm(numel(pool), min(P.nNeigh, numel(pool)))));
        r.g300 = grid_pairs(W, sp, 300, X1, X2, D(r.nbr, :), 'encounter_rule');
    end
    if r.is_css
        r.g600 = grid_pairs(W, sp, 600, X1, X2, D(r.nbr, :), 'encounter_rule');
        r.gf300 = grid_pairs(W, sp, 300, X0, X0, D(r.nbr, :), 'encounter_rule');
    end
end


% =========================================================================
function g = grid_pairs(W, sp, steps, X1, X2, D, snapMode)
% The wiring (slot 1) against its clone and against each grid mutant in D (slot 2), with the
% same K start pairs for every pairing. Per mutant: pDM, the share of encounters in which it
% defects; retal, the share of those in which the wiring defects too; d, the mutant's payoff
% minus the clone's in the same draw (T, R, P, S = 5, 3, 1, 0); nT, encounters where it earns T.
    K = size(X1, 2); nN = size(D, 1);
    Wn = repmat(W, 1, 1, nN);
    for q = 1:nN, Wn(D(q,1), D(q,2), q) = Wn(D(q,1), D(q,2), q) + D(q,3); end
    W2 = cat(3, repmat(W, 1, 1, K), Wn(:,:,repelem(1:nN, K)));
    [a1, a2] = encounter_batch(repmat(W, 1, 1, K * (nN + 1)), W2, repmat(X1, 1, nN + 1), ...
        repmat(X2, 1, nN + 1), sp, steps, snapMode);
    a1 = reshape(a1, K, nN + 1); a2 = reshape(a2, K, nN + 1);
    aR = a1(:, 2:end); aM = a2(:, 2:end);
    pay = @(r, m) 1 * (~r & ~m) + 0 * (~r & m) + 5 * (r & ~m) + 3 * (r & m);
    d = pay(aR, aM) - pay(a1(:,1), a2(:,1));
    nD = sum(~aM, 1);
    g.c1 = all(a1(:,1)); g.c2 = all(a2(:,1));                    % used when K = 1
    g.aR = aR; g.aM = aM;
    if K > 1, g = rmfield(g, {'aR', 'aM'}); end
    g.pC1clone = single(mean([a1(:,1); a2(:,1)]));
    g.pDM = single(mean(~aM, 1)); g.pDR = single(mean(~aR, 1));
    g.retal = single(sum(~aR & ~aM, 1) ./ max(nD, 1)); g.retal(nD == 0) = NaN;
    g.dbar = single(mean(d, 1)); g.sd = single(std(d, 0, 1));
    g.nT = uint8(sum(aR & ~aM, 1));
    code = 2 * aR + aM; g.multi = any(code ~= code(1,:), 1);
end


% =========================================================================
function xjobs = crosscheck_jobs(reps, sp, a, cfg)
% Control (c): the gen-500 snapshot (historyM/B, stored BEFORE evaluation) of every
% replicate of architecture a, with its recorded historyCooperation(500). Seed of replicate r:
% rng(rngSeed + 5e6 + 1000 a + r, 'twister'), rngSeed 20260927 by default.
    xjobs = cell(1, numel(reps));
    for r = 1:numel(reps)
        R = reps{r};
        assert(R.historyGenerations(end) == 500);
        xjobs{r} = struct('arch', a, 'rep', r, 'M', R.historyM(:,:,:,end), ...
            'B', R.historyB(:,:,:,end), 'sp', sp, 'obs', double(R.historyCooperation(500)), ...
            'succ', double(R.historyCooperation(end)) >= 0.95, ...
            'seed', cfg.rngSeed + 5e6 + 1000 * a + r);
    end
end


% =========================================================================
function x = crosscheck_rep(job, K)
% Predicted focal cooperation under the runs' uniform pairing, integrated exactly
% over distinct wirings: P(focal u meets v) = (c_v - [u == v]) / (N - 1).
    N = size(job.M, 3);
    [Wu, ~, ic] = unique(reshape(double(job.M) .* double(job.B), 16, N)', 'rows');
    c = accumarray(ic, 1); nu = numel(c);
    Pm = (repmat(c', nu, 1) - eye(nu)) / (N - 1);
    Wr = reshape(Wu', 4, 4, nu);
    [iu, iv] = ndgrid(1:nu, 1:nu); iu = iu(:); iv = iv(:);
    rng(job.seed, 'twister');
    ridx = repelem(1:numel(iu), K);
    X1 = rand(4, numel(ridx)) * 0.01; X2 = rand(4, numel(ridx)) * 0.01;
    a1 = encounter_batch(Wr(:,:,iu(ridx)), Wr(:,:,iv(ridx)), X1, X2, job.sp, 300, 'encounter_rule');
    Aeng = reshape(mean(reshape(a1, K, []), 1), nu, nu);
    X0 = 0.005 * ones(4, numel(iu));
    a0 = encounter_batch(Wr(:,:,iu), Wr(:,:,iv), X0, X0, job.sp, 600, 'simpair');
    Acls = reshape(double(a0), nu, nu);
    x = struct('arch', job.arch, 'rep', job.rep, 'obs', job.obs, 'succ', job.succ, 'nu', nu, ...
        'predEng', sum(c .* sum(Pm .* Aeng, 2)) / N, 'predCls', sum(c .* sum(Pm .* Acls, 2)) / N);
end


% =========================================================================
%  AGGREGATION
% =========================================================================
function V = variants()
% Membership under the runs' rule, variants {name, theta, payoff rule}; row 1 is the reference variant.
% A wiring is a member when it meets (C1) (p_C1 >= theta) and (C2): it defects back in at least
% theta of the encounters with every sampled mutant that defects in at least half, and no
% sampled mutant exceeds (below). Payoff rules: 'test', mean d > 0 and above 1.6449 sd / sqrt(K)
% (one-sided 0.05, paired); 'strict', mean d > 0; 'Troute', as 'test' and the mutant earned T
% at least once. theta = 1.0 is a control that must move the share of members (control (d)).
    V = {'primary 0.5/test', 0.5, 'test'; 'sens 0.9/test', 0.9, 'test'; ...
         'broken 1.0/test', 1.0, 'test'; 'strict 0.5', 0.5, 'strict'; ...
         'Troute 0.5', 0.5, 'Troute'; 'strict 0.9', 0.9, 'strict'};
end


function m = enc_rule_member(r, theta, rule, K)
    m = false;
    if r.pC1 < theta || isempty(r.g300), return; end
    m = enc_rule_c2(r.g300, theta, rule, K);
end


function c2 = enc_rule_c2(g, theta, rule, K)
    def = g.pDM >= 0.5;
    retOK = all(g.retal(def) >= theta);
    c2 = retOK && ~any(exceeds(g, rule, K));
end


function ex = exceeds(g, rule, K)
    db = double(g.dbar); sd = double(g.sd);
    ex = db > 0;
    if ~strcmp(rule, 'strict'), ex = ex & db > 1.6449 * sd / sqrt(K); end
    if strcmp(rule, 'Troute'), ex = ex & g.nT > 0; end
end


function U = unit_table(res, idx, K)
% One row per distinct wiring (jobs idx): the scalars every item reads.
    V = variants(); n = numel(idx);
    U = struct('weak', false(n,1), 'css', false(n,1), 'emc1', false(n,1), 'emcss', false(n,1), ...
        'pC1', nan(n,1), 'pCC', nan(n,1), 'pC1_600', nan(n,1), 'c1f300', false(n,1), ...
        'multi300', false(n,1), 'multi600', false(n,1), 'enc_rule', false(n, size(V,1)));
    for q = 1:n
        r = res{idx(q)};
        U.weak(q) = r.is_weak; U.css(q) = r.is_css; U.emc1(q) = r.em_c1; U.emcss(q) = r.em_css;
        U.pC1(q) = r.pC1; U.pCC(q) = r.pCC; U.pC1_600(q) = r.pC1_600; U.c1f300(q) = r.c1_f300;
        U.multi300(q) = r.multi300; U.multi600(q) = r.multi600;
        for v = 1:size(V,1), U.enc_rule(q,v) = enc_rule_member(r, V{v,2}, V{v,3}, K); end
    end
end


function A = aggregate(out)
% Per-condition tables: U (distinct wirings) and per-sample index into U.
    A = struct([]);
    for i = 1:numel(out.C)
        s = out.samp(i); nU = max(s.ic);
        idx = s.firstJob + (0:nU-1);
        A(i).U = unit_table(out.res, idx, out.K); A(i).idx = idx;
        A(i).ic = s.ic; A(i).rep = s.rep; A(i).succ = s.succ(:); A(i).nRep = s.nRep;
        A(i).nTot = accumarray(s.rep, 1, [s.nRep 1]);
        A(i).nCls = accumarray(s.rep, double(A(i).U.css(s.ic)), [s.nRep 1]);
        A(i).nEncRule = zeros(s.nRep, size(A(i).U.enc_rule, 2));
        for v = 1:size(A(i).U.enc_rule, 2)
            A(i).nEncRule(:,v) = accumarray(s.rep, double(A(i).U.enc_rule(s.ic, v)), [s.nRep 1]);
        end
    end
end


% =========================================================================
%  REPORT
% =========================================================================
function print_report(out)
    fprintf('\nMEMBERSHIP TEST vs THE RUNS'' ENCOUNTER RULE: AGREEMENT ON MEMBERSHIP  (run %s UTC, %.1f min, K = %d draws, seed %d)\n', ...
        out.date, out.runMin, out.K, out.cfg.rngSeed);
    fprintf('Membership test (test): x0 = 0.005 both players, 600 steps, grid of 160 mutants (membership_test).\n');
    fprintf('Runs'' encounter rule (rule): x0 ~ U(0, 0.01) per player per draw, 300 steps, the runs'' snapshot rule; K draws, common starts per sampled wiring.\n');
    fprintf('Workers: %d; C2 on <= %d sampled defecting grid mutants per member. Stages present: %s (%d of 9 settings).\n', ...
        out.cfg.nWorkers, out.cfg.nNeigh, out.stages, numel(out.C));
    if isempty(out.C), print_item5(out); return; end
    A = aggregate(out);
    print_reproduction(out, A);
    if ~isempty(out.xres), print_crosscheck(out); end
    print_item1(out, A);
    print_item2(out, A);
    print_item3(out, A);
    print_item4(out, A);
    if ~isempty(out.g1), print_item5(out); end
    if numel(out.C) == 9, print_extra_checks(out, A); end
end


% =========================================================================
%  BLOCK [P]: descriptive diagnostics, printed after the main report
% =========================================================================
function print_extra_checks(out, A)
    fprintf('\n[P] DIAGNOSTICS (descriptive)\n');
    fprintf('  (P1) membership confusion, member by the test vs member by the rule (primary) [both, test only, rule only, neither]; rule-only by the test''s class\n');
    T = zeros(1, 4); E = zeros(1, 2);
    for i = 1:numel(A)
        U = A(i).U; c = U.css(A(i).ic); e = U.enc_rule(A(i).ic, 1); w = U.weak(A(i).ic);
        t = [sum(c & e), sum(c & ~e), sum(~c & e), sum(~c & ~e)]; T = T + t;
        E = E + [sum(~c & e & w), sum(~c & e & ~w)];
        fprintf('    %-20s [%5d %5d %5d %5d]  agreement %.3f\n', out.C(i).label, t, (t(1) + t(4)) / sum(t));
    end
    fprintf('    %-20s [%5d %5d %5d %5d]  agreement %.3f; rule-only members that the test calls cooperate_only %d, C1-fail %d\n', ...
        'pooled', T, (T(1) + T(4)) / sum(T), E);
    suc = success_pct(A); nt = arrayfun(@(a) sum(a.nTot), A)';
    c1C = arrayfun(@(a) sum(a.U.weak(a.ic)), A)' ./ nt; c1E = arrayfun(@(a) sum(a.U.pC1(a.ic) >= 0.5), A)' ./ nt;
    fprintf('  (P2) what r can discriminate: r(success) of the C1 fraction by the test %+.3f; by the rule %+.3f\n', ...
        corr(c1C, suc), corr(c1E, suc));
    A2 = later_defector_coverage(out, A);
    share = arrayfun(@(a) sum(a.nEncRule(:,1)), A2)' ./ nt;
    B = boot_r(A2, out.cfg);
    fprintf('  (P3) primary membership by the rule with every grid mutant that invades under the test also tested (re-run on the same starts):\n');
    for i = 1:numel(A2)
        fprintf('    %-20s share of members by the rule %6.2f%% -> %6.2f%%  (%d distinct members re-tested, %d lost)\n', out.C(i).label, ...
            100 * sum(A(i).nEncRule(:,1)) / nt(i), 100 * share(i), A2(i).nRetest, A2(i).nLost);
    end
    fprintf('    r = %+.4f [%+.3f, %+.3f] (replicate-cluster bootstrap); architecture order kept %s\n', corr(share, suc), ...
        B.ciRE(1,:), yn(issorted(share(1:4)) && numel(unique(share(1:4))) == 4));
end


function A = later_defector_coverage(out, A)
% Re-draws each stage's jobs (same order, seeds and W as the run), finds the grid mutants that
% invade each primary runs'-rule member under the membership test's dynamics (fixed start, 600
% steps: the member cooperates, the mutant defects), and tests the ones the sample of up to 20
% missed on the member's own K start pairs. (C2) is a conjunction over mutants, so the missed
% ones can be tested separately.
    cfg = out.cfg; C = setting_table(); D = grid_deltas(); K = out.K; X0 = 0.005 * ones(4, 1);
    stages = {1:4, 0; 5:9, 1e6};
    for s = 1:2
        [~, jobs] = draw_samples(C, stages{s,1}, cfg, stages{s,2});
        for i = stages{s,1}
            a = A(i); a.nRetest = 0; a.nLost = 0; U = a.U;
            for q = find(U.enc_rule(:,1))'
                j = a.idx(q); r = out.res{j}; job = jobs{j - A(stages{s,1}(1)).idx(1) + 1};
                assert(isequal(job.cond, i) && job.u == q);
                F = grid_pairs(job.W, job.sp, 600, X0, X0, D, 'simpair');
                miss = setdiff(find(F.aR & ~F.aM), r.nbr);
                if isempty(miss), continue; end
                a.nRetest = a.nRetest + 1;
                rng(job.seed, 'twister'); X1 = rand(4, K) * 0.01; X2 = rand(4, K) * 0.01;
                if ~enc_rule_c2(grid_pairs(job.W, job.sp, 300, X1, X2, D(miss, :), 'encounter_rule'), 0.5, 'test', K)
                    U.enc_rule(q, 1) = false; a.nLost = a.nLost + 1;
                end
            end
            a.U = U; a.nEncRule(:,1) = accumarray(a.rep, double(U.enc_rule(a.ic, 1)), [a.nRep 1]);
            A2(i) = a; %#ok<AGROW>
        end
    end
    A = A2;
end


% =========================================================================
function print_reproduction(out, A)
    fprintf('\n[0] CONTROLS (a) the share of members by the membership test and (b) batch emulation of the membership test\n');
    fprintf('  %-20s %8s %8s | %8s | %6s %6s | %s\n', 'setting', 'n_m', 'n', ...
        'share', 'succ', 'd.wir', 'emulation mismatches C1 / member (distinct)');
    mm = 0; nu = 0;
    for i = 1:numel(out.C)
        c = out.C(i); U = A(i).U;
        nv = sum(A(i).nCls); nt = sum(A(i).nTot);
        e1 = sum(U.emc1 ~= U.weak); e2 = sum(U.emcss ~= U.css); mm = mm + e1 + e2; nu = nu + numel(U.css);
        sr = 100 * mean(A(i).succ);
        fprintf('  %-20s %8d %8d | %7.2f%% | %5.0f%% %6d | %d / %d\n', c.label, nv, nt, ...
            100 * nv / nt, sr, numel(U.css), e1, e2);
    end
    share = arrayfun(@(a) sum(a.nCls) / sum(a.nTot), A)'; suc = success_pct(A);
    fprintf('  r(share of members, success) over these %d settings = %+.4f\n', numel(share), corr(share, suc));
    okB = mm / nu <= 0.001;
    fprintf('  (b) %d mismatches over %d distinct wirings: %s\n', mm, nu, pf(okB));
end


% =========================================================================
function print_crosscheck(out)
    fprintf('\n[C] CONTROL (c): predicted focal cooperation vs recorded historyCooperation(500), gen-500 snapshot, exact pairing integral\n');
    fprintf('  %-20s %4s | %-40s | %-40s\n', 'architecture', 'reps', 'runs'' rule: mean|obs-pred|  med|z|  %|z|>3', ...
        'membership test dynamics: same');
    X = [out.xres{:}]; names = {'graded separate', 'graded contrast', 'saturating separate', 'saturating contrast'};
    for a = 1:4
        x = X([X.arch] == a); obs = 1 - [x.obs];
        [e1, z1, s1] = zstats(obs, 1 - [x.predEng]); [e2, z2, s2] = zstats(obs, 1 - [x.predCls]);
        fprintf('  %-20s %4d | %8.2f pts %8.2f %8.1f%%              | %8.2f pts %8.2f %8.1f%%\n', names{a}, numel(x), ...
            e1, z1, s1, e2, z2, s2);
    end
    fprintf('  (z = (obs - pred) / sqrt(max(pred(1-pred), 1e-3) / N), N = 1000: one realized pairing per generation.\n');
    fprintf('   Sampling noise alone gives a median |z| of about 0.67.)\n');
end


function [e, z, s] = zstats(obs, pred)
    zz = (obs - pred) ./ sqrt(max(pred .* (1 - pred), 1e-3) / 1000);
    e = 100 * mean(abs(obs - pred)); z = median(abs(zz)); s = 100 * mean(abs(zz) > 3);
end


% =========================================================================
function [P, Cc, Cl, D, Dc, PCC] = pooled_c1(A)
% Per sampled individual: p_C1 and the membership test's class (1 member, 2 cooperates but does
% not retaliate, 3 fails (C1)), pooled over settings; Cl, the cluster (setting x replicate).
    cat = @(U) 1 * (U.css) + 2 * (U.weak & ~U.css) + 3 * (~U.weak);
    P = []; Cc = []; Cl = []; D = []; Dc = []; PCC = [];
    for i = 1:numel(A)
        U = A(i).U; cu = cat(U);
        P = [P; U.pC1(A(i).ic)]; PCC = [PCC; U.pCC(A(i).ic)]; Cc = [Cc; cu(A(i).ic)]; %#ok<AGROW>
        Cl = [Cl; i * 1000 + A(i).rep]; D = [D; U.pC1]; Dc = [Dc; cu]; %#ok<AGROW>
    end
end


function print_item1(out, A)
    fprintf('\n[1] C1 AGREEMENT: p_C1 under the runs'' rule = P(cooperate with own clone), %d draws x 2 players\n', out.K);
    fprintf('  Per sampled wiring (copy-weighted, as the share of members is), 9 settings pooled; "distinct" = each wiring once.\n');
    [P, Cc, Cl, D, Dc, PCC] = pooled_c1(A);
    names = {'membership test member', 'cooperate_only', 'C1 fail'};
    fprintf('  %-18s %6s %6s | %6s %15s %15s | %6s %6s %6s | %s\n', 'class', 'n', 'dist.', 'median', ...
        'IQR', 'range', '>=0.5', '>=0.9', '=1', 'share >= 0.5 [cluster CI]  distinct median');
    for k = 1:3
        p = P(Cc == k); d = D(Dc == k);
        ci = cluster_ci(double(p >= 0.5), Cl(Cc == k), 2000, 11);
        fprintf('  %-18s %6d %6d | %6.3f [%5.3f, %5.3f] [%5.3f, %5.3f] | %6.3f %6.3f %6.3f | %.3f [%.3f, %.3f]  %.3f\n', ...
            names{k}, numel(p), numel(d), median(p), quantile(p, [.25 .75]), min(p), max(p), mean(p >= 0.5), ...
            mean(p >= 0.9), mean(p == 1), mean(p >= 0.5), ci, median(d));
    end
    fprintf('  members: P(CC) median %.3f; mean p_C1 %.3f\n', median(PCC(Cc == 1)), mean(P(Cc == 1)));
    for th = [0.5 0.9]
        cls = Cc <= 2; enc_rule = P >= th;
        fprintf('  confusion, rule threshold %.1f [test C1 & rule C1, test C1 & not, not & rule C1, neither] = [%d %d %d %d]  agreement %.4f, kappa %.3f\n', ...
            th, sum(cls & enc_rule), sum(cls & ~enc_rule), sum(~cls & enc_rule), sum(~cls & ~enc_rule), mean(cls == enc_rule), kappa(cls, enc_rule));
    end
    fprintf('  per setting: members median p_C1 | members share >= 0.5 | cooperate_only median | C1-fail share >= 0.5 | n members\n');
    for i = 1:numel(A)
        [p, cc] = pooled_c1(A(i));
        fprintf('    %-20s %6.3f | %6.3f | %6.3f | %6.3f | %d\n', out.C(i).label, median(p(cc == 1)), ...
            mean(p(cc == 1) >= 0.5), median(p(cc == 2)), mean(p(cc == 3) >= 0.5), sum(cc == 1));
    end
end


function k = kappa(a, b)
    po = mean(a == b); pe = mean(a) * mean(b) + mean(~a) * mean(~b);
    k = (po - pe) / (1 - pe);
end


function ci = cluster_ci(v, cl, nB, seed)
% Percentile CI of mean(v), resampling clusters (condition x replicate).
    st = rng(); rng(seed);
    [~, ~, g] = unique(cl); s = accumarray(g, v(:)); n = accumarray(g, 1); G = numel(s);
    b = randi(G, G, nB); m = sum(s(b), 1) ./ sum(n(b), 1);
    ci = prctile(m, [2.5 97.5]); rng(st);
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end


% =========================================================================
function m = member_c2(r, K)
% A member (membership test) against the sampled grid mutants (up to 20, drawn from the
% defecting pool) that defect against it under the membership test's dynamics.
    g = r.g300; dc = r.clsDefect(r.nbr);
    rt = double(g.retal(dc)); rt = rt(~isnan(rt));
    m = struct('nDef', sum(dc), 'nClsDef', sum(r.clsDefect), 'pool', r.pool, 'nSamp', numel(r.nbr), ...
        'meanRetal', mean(rt), 'minRetal', min([rt, NaN]), ...
        'allRetal05', all(rt >= 0.5), 'stillDef', mean(g.pDM(dc) >= 0.5), ...
        'pDM', mean(double(g.pDM(dc))), 'pDR', mean(double(g.pDR(dc))), ...
        'exTest', any(exceeds(g, 'test', K)), 'exStrict', any(exceeds(g, 'strict', K)), ...
        'exT', any(exceeds(g, 'Troute', K)), 'exDefCls', any(exceeds(g, 'test', K) & dc), ...
        'c2', enc_rule_c2(g, 0.5, 'test', K), 'c2_09', enc_rule_c2(g, 0.9, 'test', K), ...
        'multi', mean(g.multi(dc)), 'retalAll', rt);
    if isempty(rt), m.meanRetal = NaN; m.allRetal05 = NaN; end
end


function [M, cl, R] = member_table(res, jobIdx, K)
% member_c2 for every job in jobIdx (one row per SAMPLE, copies repeated).
    [uj, ~, back] = unique(jobIdx);
    mu = arrayfun(@(j) member_c2(res{j}, K), uj, 'uni', 0); mu = [mu{:}];
    M = mu(back); cl = []; R = [];
    for q = 1:numel(uj), R = [R, mu(q).retalAll]; end %#ok<AGROW>
end


function print_member_c2(M, cl, R, label)
    f = @(x) x(~isnan(x));
    mr = f([M.meanRetal]); mn = f([M.minRetal]);
    fprintf('  %s: %d sampled members; %d with no grid mutant that defects under the test\n', label, numel(M), sum([M.nDef] == 0));
    fprintf('    grid mutants per member (of 160) that defect under the test: median %g [%g, %g]; defecting pool median %g; sampled median %g, of them defecting under the test median %g\n', ...
        median([M.nClsDef]), quantile([M.nClsDef], [.25 .75]), median([M.pool]), median([M.nSamp]), median([M.nDef]));
    fprintf('    rule: share of them still defecting in >= 0.5 of draws  mean %.3f\n', mean(f([M.stillDef])));
    fprintf('    rule: member defects back | mutant defected: member mean    median %.3f [%.3f, %.3f] range [%.3f, %.3f]\n', ...
        median(mr), quantile(mr, [.25 .75]), min(mr), max(mr));
    fprintf('            worst-case mutant (min)                  median %.3f [%.3f, %.3f]; share of members with every one >= 0.5: %.3f\n', ...
        median(mn), quantile(mn, [.25 .75]), mean(f(double([M.allRetal05]))));
    fprintf('            unconditional P(member defects) vs them mean %.3f; P(mutant defects) mean %.3f\n', ...
        mean([M.pDR], 'omitnan'), mean([M.pDM], 'omitnan'));
    fprintf('    pooled over (distinct member, mutant) pairs: retaliation median %.3f [%.3f, %.3f], share >= 0.5 %.3f, n = %d\n', ...
        median(R), quantile(R, [.25 .75]), mean(R >= 0.5), numel(R));
    fprintf('    some mutant out-earns the clone (rule): test %.3f  strict %.3f  T-route %.3f  (one that defects under the test: %.3f)\n', ...
        mean([M.exTest]), mean([M.exStrict]), mean([M.exT]), mean([M.exDefCls]));
    c2 = double([M.c2]);
    if ~isempty(cl), ci = cluster_ci(c2, cl, 2000, 12); else, ci = [NaN NaN]; end
    fprintf('    C2 PASS UNDER THE RULE among members by the test: %.3f [%.3f, %.3f] at 0.5; %.3f at 0.9\n', mean(c2), ci, mean([M.c2_09]));
end


function print_item2(out, A)
    fprintf('\n[2] C2 AGREEMENT: members by the test against the grid mutants that defect on them under the test (runs'' rule, K = %d)\n', out.K);
    J = []; cl = [];
    for i = 1:numel(A)
        mem = A(i).U.css(A(i).ic);
        J = [J; A(i).idx(A(i).ic(mem))']; cl = [cl; i * 1000 + A(i).rep(mem)]; %#ok<AGROW>
    end
    [M, ~, R] = member_table(out.res, J, out.K);
    print_member_c2(M, cl, R, sprintf('%d settings pooled, copy-weighted', numel(A)));
    fprintf('  per setting: members | mean retaliation | min retaliation median | C2 pass under the rule\n');
    for i = 1:numel(A)
        mem = A(i).U.css(A(i).ic); Mi = member_table(out.res, A(i).idx(A(i).ic(mem))', out.K);
        if isempty(Mi), fprintf('    %-20s 0\n', out.C(i).label); continue; end
        fprintf('    %-20s %5d | %.3f | %.3f | %.3f\n', out.C(i).label, numel(Mi), mean([Mi.meanRetal], 'omitnan'), ...
            median([Mi.minRetal], 'omitnan'), mean([Mi.c2]));
    end
end


% =========================================================================
function print_item3(out, A)
    fprintf('\n[3] CAUSE: start rule vs step count (same K start pairs for the two random-start rows)\n');
    [P, Cc] = pooled_c1(A); P6 = []; F3 = []; M3 = []; M6 = [];
    for i = 1:numel(A)
        U = A(i).U; P6 = [P6; U.pC1_600(A(i).ic)]; F3 = [F3; U.c1f300(A(i).ic)]; %#ok<AGROW>
        M3 = [M3; U.multi300(A(i).ic)]; M6 = [M6; U.multi600(A(i).ic)]; %#ok<AGROW>
    end
    cls = Cc <= 2;
    fprintf('  C1 outcome vs the membership test (fixed start, 600 steps), all sampled wirings (n = %d):\n', numel(P));
    rows = {'fixed start, 300 steps', F3; 'random start, 600 steps', P6 >= 0.5; 'random start, 300 steps (the runs'' encounter rule)', P >= 0.5};
    for k = 1:3
        e = rows{k,2};
        fprintf('    %-34s agreement %.4f  kappa %.3f  | members C1 %.3f  cooperate_only C1 %.3f  C1-fail C1 %.3f\n', rows{k,1}, ...
            mean(e == cls), kappa(cls, e), mean(e(Cc == 1)), mean(e(Cc == 2)), mean(e(Cc == 3)));
    end
    fprintf('    members median p_C1: random 600 %.3f, random 300 %.3f\n', median(P6(Cc == 1)), median(P(Cc == 1)));
    fprintf('  Multistability (share of clone pairings with more than one outcome over K draws):\n');
    nm = {'member', 'cooperate_only', 'C1 fail'};
    for k = 1:3
        fprintf('    %-8s random 300 %.3f   random 600 %.3f\n', nm{k}, mean(M3(Cc == k)), mean(M6(Cc == k)));
    end
    print_item3_c2(out, A);
end


function print_item3_c2(out, A)
% Members only: the same (C2) statistics under the three dynamics.
    J = [];
    for i = 1:numel(A), J = [J; A(i).idx(A(i).ic(A(i).U.css(A(i).ic)))']; end %#ok<AGROW>
    [uj, ~, back] = unique(J); K = out.K;
    fld = {'gf300', 'fixed start, 300 steps', 1; 'g600', 'random start, 600 steps', K; 'g300', 'random start, 300 steps', K};
    fprintf('  C2 for members by the test (copy-weighted, n = %d): retaliation vs grid mutants that defect under the test | C2 pass under the rule 0.5 | pairings multistable\n', numel(J));
    for f = 1:3
        v = nan(numel(uj), 3);
        for q = 1:numel(uj)
            r = out.res{uj(q)}; g = r.(fld{f,1}); dc = r.clsDefect(r.nbr);
            v(q,:) = [mean(double(g.retal(dc)), 'omitnan'), enc_rule_c2(g, 0.5, 'test', fld{f,3}), mean(g.multi(dc))];
        end
        v = v(back, :);
        fprintf('    %-26s mean %.3f | %.3f | %.3f\n', fld{f,2}, mean(v(:,1), 'omitnan'), mean(v(:,2)), mean(v(:,3), 'omitnan'));
    end
end


% =========================================================================
function print_item4(out, A)
    V = variants(); nV = size(V, 1); n = numel(A);
    suc = success_pct(A); nt = arrayfun(@(a) sum(a.nTot), A)';
    shareC = arrayfun(@(a) sum(a.nCls), A)' ./ nt;
    shareE = cell2mat(arrayfun(@(a) sum(a.nEncRule, 1), A', 'uni', 0)) ./ nt;
    B = boot_r(A, out.cfg);
    fprintf('\n[4] SHARE OF MEMBERS UNDER THE RUNS'' RULE. Share %% [95%% cluster-bootstrap CI over replicates]\n');
    fprintf('  %-20s %5s | %-22s | %-22s %s\n', 'setting', 'succ', 'membership test (paper)', 'runs'' rule (primary)', ...
        strjoin(V(2:end,1)', ' | '));
    for i = 1:n
        fprintf('  %-20s %4.0f%% | %6.2f [%5.2f, %5.2f]   | %6.2f [%5.2f, %5.2f]  ', out.C(i).label, suc(i), ...
            100 * shareC(i), 100 * B.ciShareC(i,:), 100 * shareE(i,1), 100 * B.ciShareE(i,:,1));
        fprintf(' %6.2f', 100 * shareE(i, 2:end)); fprintf('\n');
    end
    rC = corr(shareC, suc);
    fprintf('\n  r(share of members, success) over these %d settings, CI: replicate-cluster bootstrap | setting bootstrap | Fisher z\n', n);
    fprintf('  %-20s %+.4f [%+.3f, %+.3f] | [%+.3f, %+.3f] | [%+.3f, %+.3f]\n', 'membership test', rC, B.ciRC, ...
        cond_boot(shareC, suc), fisher(rC, n));
    for v = 1:nV
        rE = corr(shareE(:,v), suc);
        fprintf('  %-20s %+.4f [%+.3f, %+.3f] | [%+.3f, %+.3f] | [%+.3f, %+.3f]  d(r) %+.4f [%+.3f, %+.3f]  arch order kept %s (boot %.3f)\n', ...
            V{v,1}, rE, B.ciRE(v,:), cond_boot(shareE(:,v), suc), fisher(rE, n), rE - rC, B.ciD(v,:), ...
            yn(issorted(shareE(1:4, v)) && numel(unique(shareE(1:4,v))) == 4), B.pOrder(v));
    end
    print_family_r(out, shareC, shareE(:,1), suc);
    fprintf('\n  CONTROL (d): the broken threshold 1.0 moves the share of members by %+.2f points on average (primary to broken), max %.2f\n', ...
        100 * mean(shareE(:,3) - shareE(:,1)), 100 * max(abs(shareE(:,3) - shareE(:,1))));
    fprintf('  rule / test ratio of the share of members (primary): median %.3f, range [%.3f, %.3f]\n', ...
        median(shareE(:,1) ./ shareC), min(shareE(:,1) ./ shareC), max(shareE(:,1) ./ shareC));
end


function B = boot_r(A, cfg)
% Replicate-cluster bootstrap within every setting; the membership test and every runs'-rule
% variant from the same resamples, so the difference in r is paired. Resamples from
% rng(bootSeed) (20260928 by default, the session's generator); the caller's stream is restored.
    st = rng(); rng(cfg.bootSeed); nB = cfg.nBoot; n = numel(A); nV = size(A(1).nEncRule, 2);
    oC = zeros(n, nB); oE = zeros(n, nB, nV); sU = zeros(n, nB);
    for i = 1:n
        b = randi(A(i).nRep, A(i).nRep, nB); T = sum(A(i).nTot(b), 1);
        oC(i,:) = sum(A(i).nCls(b), 1) ./ T; sU(i,:) = mean(A(i).succ(b), 1);
        for v = 1:nV, e = A(i).nEncRule(:,v); oE(i,:,v) = sum(e(b), 1) ./ T; end
    end
    rc = colcorr(oC, sU); B.ciRC = prctile(rc, [2.5 97.5]);
    B.ciShareC = prctile(oC, [2.5 97.5], 2);
    for v = 1:nV
        re = colcorr(oE(:,:,v), sU);
        B.ciRE(v,:) = prctile(re, [2.5 97.5]); B.ciD(v,:) = prctile(re - rc, [2.5 97.5]);
        B.pOrder(v) = mean(all(diff(oE(1:4,:,v), 1, 1) > 0, 1));
        B.ciShareE(:,:,v) = prctile(oE(:,:,v), [2.5 97.5], 2);
    end
    rng(st);
end


function r = colcorr(X, Y)
    X = X - mean(X, 1); Y = Y - mean(Y, 1);
    r = sum(X .* Y, 1) ./ sqrt(sum(X.^2, 1) .* sum(Y.^2, 1));
end


function ci = cond_boot(x, y)
    st = rng(); rng(31415); nB = 10000; n = numel(x); b = randi(n, n, nB);
    ok = arrayfun(@(k) numel(unique(b(:,k))) >= 3, 1:nB);
    r = colcorr(x(b(:,ok)), y(b(:,ok))); ci = prctile(r, [2.5 97.5]); rng(st);
end


function ci = fisher(r, n)
    if n > 3, ci = tanh(atanh(r) + [-1 1] * 1.96 / sqrt(n - 3)); else, ci = [NaN NaN]; end
end
function s = yn(t), if t, s = 'yes'; else, s = 'NO'; end, end


function print_family_r(out, shareC, shareE, suc)
    fam = {'signaling', 'spatial'};
    fprintf('  within-family r (descriptive): ');
    for f = 1:2
        m = strcmp({out.C.family}', fam{f});
        if sum(m) < 3, continue; end
        fprintf('%s test %+.3f rule %+.3f;  ', fam{f}, corr(shareC(m), suc(m)), corr(shareE(m), suc(m)));
    end
    fprintf('\n');
end


% =========================================================================
function print_item5(out)
    fprintf('\n[5] GENERATION 1: %d random wirings per architecture (cooperation_condition_shares.m sample, rng(2026), randn(4) * 0.1)\n', out.g1(1).n);
    V = variants(); K = out.K;
    fprintf('  %-20s | %-13s %-13s | %6s | %-22s %-22s %s\n', 'architecture', 'test C1 %', 'test mem %', 'emul.', ...
        'rule C1 >= 0.5 %', 'rule mem % [Wilson]', strjoin(V(2:end,1)', ' | '));
    for a = 1:4
        idx = out.g1(a).firstJob + (0:out.g1(a).n - 1);
        U = unit_table(out.res, idx, K); n = numel(idx);
        [lo, hi] = wilson(sum(U.enc_rule(:,1)), n);
        fprintf('  %-20s | %6.2f        %5.2f         | %6d | %6.2f                 %5.2f [%4.2f, %4.2f]  ', out.g1(a).name, ...
            100 * mean(U.weak), 100 * mean(U.css), sum(U.emc1 ~= U.weak) + sum(U.emcss ~= U.css), ...
            100 * mean(U.pC1 >= 0.5), 100 * mean(U.enc_rule(:,1)), 100 * lo, 100 * hi);
        fprintf(' %5.2f', 100 * mean(U.enc_rule(:, 2:end), 1)); fprintf('\n');
    end
    for a = 1:4
        idx = out.g1(a).firstJob + (0:out.g1(a).n - 1);
        U = unit_table(out.res, idx, K);
        cc = 1 * U.css + 2 * (U.weak & ~U.css) + 3 * (~U.weak);
        fprintf('  %-20s p_C1 median [IQR] member %.3f [%.3f, %.3f] (n %d) | cooperate_only %.3f | C1 fail share >= 0.5 %.3f | members C1 >= 0.5 %.3f\n', ...
            out.g1(a).name, median(U.pC1(cc == 1)), quantile(U.pC1(cc == 1), [.25 .75]), sum(cc == 1), ...
            median(U.pC1(cc == 2)), mean(U.pC1(cc == 3) >= 0.5), mean(U.pC1(cc == 1) >= 0.5));
        cls = U.weak; enc_rule = U.pC1 >= 0.5;
        fprintf('  %-20s C1 confusion [both, test only, rule only, neither] = [%d %d %d %d]; clone multistable: member %.3f cooperate_only %.3f fail %.3f\n', '', ...
            sum(cls & enc_rule), sum(cls & ~enc_rule), sum(~cls & enc_rule), sum(~cls & ~enc_rule), mean(U.multi300(cc == 1)), ...
            mean(U.multi300(cc == 2)), mean(U.multi300(cc == 3)));
        if any(U.css)
            [M, ~, R] = member_table(out.res, idx(U.css)', K);
            print_member_c2(M, [], R, sprintf('%s gen-1 members', out.g1(a).name));
        end
    end
end


function [lo, hi] = wilson(k, n)
    p = k / n; z = 1.96; d = 1 + z^2 / n;
    c = (p + z^2 / (2*n)) / d; h = z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d;
    lo = c - h; hi = c + h;
end
