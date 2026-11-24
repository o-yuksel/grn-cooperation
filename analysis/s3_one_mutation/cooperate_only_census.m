function out = cooperate_only_census(cfg)
% COOPERATE_ONLY_CENSUS  Runs the per-mutation census on the wirings that cooperate with their
% clone but are not members (cooperators that do not retaliate; group 'cooperate_only'), or again on
% members (group 'member', for the checks).
%
% Produces : SI S3.3, "Cooperators that do not retaliate lose cooperation less often than
%            members" (the rows cooperate_only_losses reads)
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures)
%            cfg: dataDir (the publication folder; required), outFile, `modes` (architectures 1
%              to 4), reps (one replicate list per architecture, [] for all 100), group,
%              `ejectDraws`, and the census settings K (30), nSnaps (5), steadyStateFrom (100),
%              pToggle (0.2), sdWeight (0.2), rngSeed, useParallel, minXnorm (0.5), maxXnorm
%              (14), minR (1e-6)
% Outputs  : cfg.outFile (publication/cooperate_only_census/full_cooperate_only.mat in the paper's run):
%              W, one row per scored wiring with the census's fields (readout r, ||x*||, alpha,
%              beta, the K draws and whether each ends cooperation, the wiring W), cfg, snapshot
%              generations, timings
% Run      : through run_cooperate_only_census('full')   (6.8 min on 2 workers)
% Seed     : rng(rngSeed + 1000 m + r + 500000 g) for architecture m and replicate r, g = 1 for
%            'cooperate_only' and 0 for 'member'; rngSeed 20261016 unless set; generator not named (in a
%            pool, the worker's default, Threefry)
%
% The code is per_mutation_census's (snapshots, mutation distribution, measures), with the same
% rule for kept rows (ok): readout against the clone at least minR, ||x*|| in [minXnorm,
% maxXnorm]. With `ejectDraws` false (the default for 'cooperate_only') no draw is tested for membership
% (`ejected` is NaN); membership_test restores the random stream, so the draws and the
% cooperation losses are the same either way.

    cfg = fill_defaults(cfg);
    sensing_modes = { ...
        struct('name','graded_separate',     'sensingMode','separate',     'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','graded_contrast',     'sensingMode','contrast', 'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','saturating_separate', 'sensingMode','separate',     'hillCoefficient',1,  'discreteOutput',false), ...
        struct('name','saturating_contrast', 'sensingMode','contrast', 'hillCoefficient',1,  'discreteOutput',false)};

    jobs = {};          % one job per (architecture, replicate)
    snapGens = cell(4,1);
    for mi = 1:numel(cfg.modes)
        m = cfg.modes(mi);
        mn = sensing_modes{m}.name;
        d = load(fullfile(cfg.dataDir, sprintf('%s_100reps.mat', mn)), 'results');
        reps = cfg.reps{mi};
        if isempty(reps), reps = 1:numel(d.results); end
        snaps = extract_snapshots(d.results, reps, cfg);
        clear d;
        snapGens{m} = unique(cell2mat(cellfun(@(s) [s.S.gen], snaps, 'UniformOutput', false)));
        fprintf('%s: %d reps, snapshot generations %s\n', mn, numel(reps), mat2str(snapGens{m}));
        [params, sim_params] = mode_params(sensing_modes{m});
        for k = 1:numel(snaps)
            jobs{end+1} = struct('m', m, 'mode', mn, 'snap', snaps{k}, ...
                'params', params, 'sim_params', sim_params); %#ok<AGROW>
        end
    end

    nJ = numel(jobs);
    rows = cell(nJ,1); tJob = nan(nJ,1);
    t0 = tic;
    if cfg.useParallel
        parfor j = 1:nJ
            tj = tic;
            rows{j} = process_rep(jobs{j}, cfg);
            tJob(j) = toc(tj);
        end
    else
        for j = 1:nJ
            tj = tic;
            rows{j} = process_rep(jobs{j}, cfg);
            tJob(j) = toc(tj);
        end
    end
    W = vertcat(rows{:});
    out = struct('W', W, 'cfg', cfg, 'snapGens', {snapGens}, ...
        'wallSec', toc(t0), 'jobSec', tJob, ...
        'jobMode', cellfun(@(j) j.m, jobs)', 'jobRep', cellfun(@(j) j.snap.rep, jobs)');
    if ~isempty(cfg.outFile)
        save(cfg.outFile, '-struct', 'out', '-v7.3');
        fprintf('Saved %s (%d rows, %.1f min wall)\n', cfg.outFile, numel(W), out.wallSec/60);
    end
end


% =========================================================================
function cfg = fill_defaults(cfg)
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    dflt = struct('dataDir', '', 'outFile', '', 'modes', 1:4, 'reps', [], ...
        'group', 'cooperate_only', 'ejectDraws', [], 'K', 30, 'nSnaps', 5, ...
        'steadyStateFrom', 100, 'pToggle', 0.2, 'sdWeight', 0.2, ...
        'rngSeed', 20261016, 'useParallel', true, ...
        'minXnorm', 0.5, 'maxXnorm', 14, 'minR', 1e-6);
    f = fieldnames(dflt);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = dflt.(f{i}); end
    end
    assert(~isempty(cfg.dataDir), 'cfg.dataDir required (the checkout holding the .mat data)');
    assert(any(strcmp(cfg.group, {'member','cooperate_only'})), 'group must be member or cooperate_only');
    if isempty(cfg.ejectDraws), cfg.ejectDraws = strcmp(cfg.group, 'member'); end
    if isempty(cfg.reps), cfg.reps = repmat({[]}, 1, numel(cfg.modes)); end
    if ~iscell(cfg.reps), cfg.reps = repmat({cfg.reps}, 1, numel(cfg.modes)); end
end


function [params, sim_params] = mode_params(sm)
% Identical to the census's two structs.
    params = struct('sensingMode', sm.sensingMode, ...
        'hillCoefficient', sm.hillCoefficient, ...
        'discreteOutput',  sm.discreteOutput, ...
        'michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0, ...
        'developmentSteps',600,'backend','simulate_pair');
    sim_params = struct('sensingMode', sm.sensingMode, ...
        'hillCoefficient', sm.hillCoefficient, ...
        'discreteOutput',  sm.discreteOutput, ...
        'michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0, ...
        'developmentSteps',600,'payoffDC',0,'sensorIntegrationRate',1.0, ...
        'externalSignal1',zeros(4,1),'externalSignal2',zeros(4,1));
end


% =========================================================================
function rows = process_rep(job, cfg)
% One architecture and replicate: test each distinct wiring with the membership test and
% score the members ('member') or the non-members ('cooperate_only'), with clone counts as weights
    isMemberRun = strcmp(cfg.group, 'member');
    seed = cfg.rngSeed + 1000*job.m + job.snap.rep + 500000*(~isMemberRun);
    rng(seed);
    rows = [];
    S = job.snap.S; params = job.params; sim_params = job.sim_params;
    for s = 1:numel(S)
        Wu = S(s).Wu; counts = S(s).counts; nU = size(Wu,1); nG = 4;
        cls = cell(nU,1); isMem = false(nU,1);
        for u = 1:nU
            cls{u} = membership_test(reshape(Wu(u,:),nG,nG), params);
            isMem(u) = cls{u}.is_css;
        end
        if isMemberRun, pick = find(isMem)'; else, pick = find(~isMem)'; end
        for u = pick
            Mi = S(s).Mfirst(:,:,u); Bi = S(s).Bfirst(:,:,u);
            rec = wiring_metrics(Mi, Bi, params, sim_params, cfg);
            if isempty(rec), continue; end
            rec.mode = job.m; rec.rep = job.snap.rep; rec.gen = S(s).gen;
            rec.count = counts(u);  rec.N = S(s).N;
            rec.is_css = cls{u}.is_css; rec.is_weak = cls{u}.is_weak;
            rec.alpha_css = cls{u}.alpha; rec.alpha_p = cls{u}.alpha_p;
            rec.beta_css = cls{u}.beta;
            rec.nDistinct = nU;
            if isempty(rows), rows = rec; else, rows(end+1) = rec; end %#ok<AGROW>
        end
    end
    rows = rows(:);
end


% =========================================================================
function rec = wiring_metrics(Mi, Bi, params, sim_params, cfg)
% The census's wiring_metrics without its slopes (the census's rho, which this group's reading does
% not use); only the per-draw membership_test call is switchable.
    rec = [];
    nG = 4; Msq = ones(nG); Kv = [];
    W = Mi .* Bi;

    % The clone pair's settled state x* and readout r = x_C - x_D; a wiring that does not
    % cooperate with its clone (r <= 0) is not scored
    [~, fs, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sim_params);
    x = fs(1:nG);
    r = x(1) - x(2);
    if ~isfinite(r) || r <= 0, return; end
    xn = norm(x);

    % alpha and beta of the readout with the partner's signal held fixed (alpha_and_beta)
    sPart = partner_signal(x, sim_params);
    th = alpha_and_beta(W, x, sPart, sim_params.degradationRate, ...
                         sim_params.michaelisConstant, [1;-1;0;0]);
    if th.valid
        alpha_SI = th.alpha; beta_SI = th.beta;
    else
        alpha_SI = NaN; beta_SI = NaN;
    end

    % K draws of the mutation distribution (Methods, "Game and evolution"): connection (g1, g2)
    % uniform; with probability pToggle the connection is toggled, else its weight gains
    % N(0, sdWeight^2). A draw ends cooperation when the mutant's readout with its own clone is
    % not positive.
    K = cfg.K;
    dWv = nan(K,1);
    slot = nan(K,1); isTog = false(K,1); ejected = nan(K,1); lostR = false(K,1);
    for t = 1:K
        Mm = Mi; Bm = Bi;
        g1 = randi(nG); g2 = randi(nG);
        tog = rand() < cfg.pToggle;
        if tog
            Mm(g1,g2) = ~Mm(g1,g2);
        else
            Bm(g1,g2) = Bm(g1,g2) + randn()*cfg.sdWeight;
        end
        Wm = Mm .* Bm;
        q = sub2ind([nG nG], g1, g2);
        dl = Wm(g1,g2) - W(g1,g2);
        dWv(t) = dl; slot(t) = q; isTog(t) = tog;
        if cfg.ejectDraws
            ejected(t) = ~membership_test(Wm, params).is_css;
        end
        [~, fm, ~] = simulate_pair(Msq, Wm, Kv, Msq, Wm, Kv, sim_params);
        lostR(t) = ~((fm(1) - fm(2)) > 0);
    end

    ok = (xn >= cfg.minXnorm) && (xn <= cfg.maxXnorm) && (r >= cfg.minR);
    rec = struct('mode',0,'rep',0,'gen',0,'count',0,'N',0, 'ok', ok, ...
        'r', r, 'xn', xn, 'alpha_SI', alpha_SI, 'beta_SI', beta_SI, ...
        'theorem_valid', th.valid, ...
        'pExit', mean(ejected), 'pLostReadout', mean(lostR), ...
        'dW', dWv, 'slot', slot, ...
        'isToggle', isTog, 'ejected', ejected, 'lostReadout', lostR, ...
        'is_css', false, 'is_weak', false, 'alpha_css', NaN, 'alpha_p', NaN, ...
        'beta_css', NaN, 'nDistinct', 0, ...
        'W', W);   % the wiring W = M .* B
end


% =========================================================================
function s = partner_signal(x, sim_params)
% The census's partner_signal, verbatim.
    nG = numel(x);
    intRate = 1.0;
    if isfield(sim_params,'sensorIntegrationRate') && ~isempty(sim_params.sensorIntegrationRate)
        intRate = sim_params.sensorIntegrationRate;
    end
    gC = 1; gD = 2; gS1 = 3; gS2 = 4;
    if isfield(sim_params,'discreteOutput') && sim_params.discreteOutput
        oC = double(x(gC) > x(gD)); oD = 1 - oC;
    elseif isfinite(sim_params.hillCoefficient)
        n = sim_params.hillCoefficient;
        cn = x(gC)^n; dn = x(gD)^n; den = cn + dn + 1e-12;
        oC = cn/den; oD = dn/den;
    else
        oC = x(gC); oD = x(gD);
    end
    s = zeros(nG,1);
    if strcmp(stored_value('sensingMode', sim_params.sensingMode), 'separate')
        s(gS1) = intRate*oC; s(gS2) = intRate*oD;
    else
        sig = intRate*(oC - oD);
        s(gS1) = sig; s(gS2) = -sig;
    end
end
