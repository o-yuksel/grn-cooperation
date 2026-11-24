function out = run_saturating_separate_to_gen1000(cfg)
% RUN_SATURATING_SEPARATE_TO_GEN1000  Does the share of members under saturating separate keep
% drifting after generation 500? Continues its 100 populations to generation 1,000 and reads the
% share of members, its slope from generation 700 on, and success.
%
% Produces : SI S3.2, "After generation 300, the share of members shows no resolved drift in
%            three architectures" (saturating separate continued to generation 1,000)
% Inputs   : publication/saturating_separate_100reps.mat   (run_signaling_architectures)
%            publication/member_share_by_generation_satsep100.mat   (the share of members to generation
%              500; share_of_members_by_generation); 'control' reads more (see control below)
% Outputs  : publication/saturating_separate_to_gen1000/main.mat   (checkpoint per replicate; a rerun
%              resumes from it), summary.txt ('summary'), control.txt ('control')
% Run      : run_saturating_separate_to_gen1000(struct('mode', 'main')), then
%            run_saturating_separate_to_gen1000(struct('mode', 'summary'))   (2 workers; not
%            timed: 3 replicates on 3 workers take about 1 min). 'control' (the default) checks
%            the stored inputs; struct('mode', 'main', 'reps', 1:3) runs replicates 1 to 3.
% Seed     : replicate r runs under rng(31000000 + 1000 r + 1, 'twister'); the seeds of
%            'summary' and 'control' are given beside those functions

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    switch cfg.mode
        case 'control', out = control(cfg);
        case 'main',    out = run_main(cfg, cfg.reps);
        case 'summary', out = summary(cfg);
        otherwise, error('unknown mode %s', cfg.mode);
    end
end


function cfg = fill_defaults(cfg)
    pub = fullfile(grn_data_file('publication'));
    def = struct('mode', 'control', 'pub', pub, 'outDir', fullfile(pub, 'saturating_separate_to_gen1000'), ...
        'file', 'saturating_separate_100reps.mat', 'nGen', 500, 'snapEvery', 20, ...
        'seedBase', 31000000, 'nWorkers', 2, 'reps', 1:100, ...
        'from', 700, 'nBoot', 10000, 'bootSeed', 20260922, 'nSub', 30);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if ~isfield(cfg, 'minFreeGB'), cfg.minFreeGB = 1.2 * cfg.nWorkers; end
    cfg.params = struct('michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, ...
        'developmentSteps', 600, 'backend', 'simulate_pair', 'sensingMode', 'separate', ...
        'hillCoefficient', 1, 'discreteOutput', false);
    cfg.ckFile = fullfile(cfg.outDir, 'main.mat');
end


% =========================================================================
function out = control(cfg)
% Checks of the inputs before the continuation: the generation-500 share of members on the
% 30-per-replicate sample against the stored count (PASS or FAIL), the stored success count
% (printed), the stored generation-500 population tested again for membership against the trace file (PASS or
% FAIL), and the stored 300-500 mean slope (printed). Reads also publication/graded_separate_100reps.mat
% and graded_contrast_100reps.mat (run_signaling_architectures) and css_predicts_sensing_gen500.mat
% (gen500_share_of_members); replays rng(42) (generator not named).
    diary_on(fullfile(cfg.outDir, 'control.txt'));
    fprintf('INPUT CHECKS (%s)\n', datestr(now, 'yyyy-mm-dd HH:MM'));
    % Check 1a (not the cooperation condition): the Table S3 generation-500 share
    % of members, replaying the rng(42) stream of gen500_share_of_members: the first two
    % architectures' draws are consumed, not tested.
    pats = {'graded_separate_100reps.mat', 'graded_contrast_100reps.mat', cfg.file};
    rng(42);
    for m = 1:3
        d = load(fullfile(cfg.pub, pats{m}), 'results');
        Ws = zeros(0, 16);
        for r = 1:numel(d.results)
            Bf = d.results{r}.finalNetworkB; Mf = d.results{r}.finalNetworkM;
            pick = randperm(size(Bf, 3), min(cfg.nSub, size(Bf, 3)));
            if m == 3
                for a = pick
                    W = double(Mf(:,:,a)) .* double(Bf(:,:,a)); Ws(end+1, :) = W(:)'; %#ok<AGROW>
                end
            end
        end
    end
    nPass = 0;
    for k = 1:size(Ws, 1), nPass = nPass + membership_test(reshape(Ws(k,:), 4, 4), cfg.params).is_css; end
    S = load(fullfile(cfg.pub, 'css_predicts_sensing_gen500.mat'));
    out.c1a = [nPass, size(Ws, 1)];
    out.c1aPass = nPass == S.n_in_Vm(3) && size(Ws, 1) == S.n_subsample(3);
    fprintf('Check 1a Table S3 generation-500 share of members (rng(42), 30 per replicate): %d/%d = %.4f; stored %d/%d = %.4f -> %s\n', ...
        nPass, size(Ws, 1), nPass / size(Ws, 1), S.n_in_Vm(3), S.n_subsample(3), S.reach_occupancy(3), pf(out.c1aPass));
    % Check 1b: stored success
    hcEnd = cellfun(@(x) x.historyCooperation(end), d.results);
    out.c1b = mean(hcEnd >= 0.95);
    fprintf('Check 1b stored success (cooperation >= 0.95 at generation 500): %d/100 = %.2f\n', ...
        sum(hcEnd >= 0.95), out.c1b);
    % Check 1c: the stored generation-500 population, every individual, against the trace file
    V = load(fullfile(cfg.pub, 'member_share_by_generation_satsep100.mat')); o = V.out;
    share500 = nan(100, 1);
    for r = 1:100
        share500(r) = share_full(d.results{r}.historyM(:,:,:,end), d.results{r}.historyB(:,:,:,end), cfg.params);
    end
    out.c1cMaxDiff = max(abs(share500 - o.occ(:, end)));
    fprintf('Check 1c stored generation-500 snapshot reclassified (all individuals, 100 replicates): mean %.4f; trace file %.4f; max |diff| %.3g -> %s\n', ...
        mean(share500), mean(o.occ(:, end)), out.c1cMaxDiff, pf(out.c1cMaxDiff == 0));
    sl = rep_slopes(double(o.gens(:))', o.occ, 300, 1000);
    out.storedSlope = mean(sl);
    fprintf('Check 1d stored mean slope over generations 300-500: %+.4f per 100 generations\n', ...
        out.storedSlope);
    out.pass = out.c1aPass && out.c1cMaxDiff == 0;
    fprintf('INPUT CHECKS OVERALL: %s\n', pf(out.pass));
    diary off;
end


% =========================================================================
function out = run_main(cfg, reps)
% One evolve_population run per replicate r, from its stored final population (finalNetworkM/B,
% the population entering generation 501) under its stored config, for 500 generations, storing
% every 20th generation and the last (global generation = 500 + local). rngSeedBase is
% 31000000 + 1000 r; evolve_population adds the replicate index 1. Each stored generation's share
% of members tests every individual, once per distinct wiring (share_full).
    d = load(fullfile(cfg.pub, cfg.file), 'results', 'config');
    jobs = cell(numel(reps), 1);
    for i = 1:numel(reps)
        r = reps(i);
        c = d.config;
        c.numGenerations = cfg.nGen; c.numReplicates = 1;
        c.historySampleInterval = cfg.snapEvery; c.alwaysSaveGenerations = [1 cfg.nGen];
        c.initPopulationM = d.results{r}.finalNetworkM; c.initPopulationB = d.results{r}.finalNetworkB;
        c.rngSeedBase = cfg.seedBase + 1000 * r;
        jobs{i} = struct('key', sprintf('r%03d', r), 'r', r, 'config', c, 'params', cfg.params);
    end
    clear d
    open_pool(cfg);
    res = run_pool(jobs, cfg.ckFile);
    out = [res{:}];
    fprintf('mean worker time per replicate %.1f s (evolution %.1f s, membership test %.1f s)\n', ...
        mean([out.sec]), mean([out.secEvo]), mean([out.sec]) - mean([out.secEvo]));
end


function r = run_job(job)
% seedUsed records whether the run's first stored population is the supplied start
    t = tic;
    R = evolve_population(job.config); x = R{1};
    r.secEvo = toc(t);
    r.key = job.key; r.r = job.r;
    r.hc = x.historyCooperation(:)';
    r.gens = 500 + double(x.historyGenerations(:)');      % global generation
    r.seedUsed = isequal(x.historyM(:,:,:,1), logical(job.config.initPopulationM)) && ...
        isequal(x.historyB(:,:,:,1), single(job.config.initPopulationB));
    nK = numel(r.gens); r.occ = nan(1, nK);
    for k = 1:nK
        r.occ(k) = share_full(x.historyM(:,:,:,k), x.historyB(:,:,:,k), job.params);
    end
    r.finalM = x.finalNetworkM; r.finalB = x.finalNetworkB;
    r.sec = toc(t);
end


% =========================================================================
function out = summary(cfg)
% The summary of the continued runs: the replicate-mean share of members at each stored
% generation from 300 on, the mean slope from cfg.from (700) to 1,000 with a cluster bootstrap,
% the share at 1,000 and its mean over 700-1,000, the 30-per-replicate sample of Table S3 on the
% final populations, success at 1,000 and the fate transitions, and the autocorrelation time.
% The last line says whether the slope interval includes zero. Seeds: the
% bootstrap draws from RandStream('mt19937ar', 'Seed', 20260922); the 30 wirings per replicate
% are sampled after rng(42) (generator not named).
    C = load(cfg.ckFile, 'ckKey', 'ckRes');
    res = [C.ckRes{:}];
    [~, o] = sort([res.r]); res = res(o);
    assert(numel(res) == 100 && isequal([res.r], 1:100), 'main run incomplete: %d replicates', numel(res));
    assert(all([res.seedUsed]), 'a seeded start was not used');
    diary_on(fullfile(cfg.outDir, 'summary.txt'));
    fprintf('SUMMARY (%s): saturating separate continued from generation 500 to 1,000, 100 replicates\n', datestr(now, 'yyyy-mm-dd HH:MM'));
    g = res(1).gens; share = vertcat(res.occ);
    V = load(fullfile(cfg.pub, 'member_share_by_generation_satsep100.mat')); vo = V.out;
    st = RandStream('mt19937ar', 'Seed', cfg.bootSeed);
    fprintf('\nreplicate-mean share of members (all individuals) at snapshot generations:\n');
    show = [vo.gens(vo.gens >= 300), g];
    allShare = [vo.occ(:, vo.gens >= 300), share];
    for k = 1:numel(show), fprintf('  gen %4d  %.4f\n', show(k), mean(allShare(:, k))); end
    % Mean of the per-replicate slopes over 700-1000, per 100 generations; whole replicates
    % resampled
    sl = rep_slopes(g, share, cfg.from, 1000);
    [m, ci] = boot_mean(sl, cfg.nBoot, st);
    out.slope = m; out.slopeCI = ci; out.slopes = sl;
    fprintf('\nmean slope over generations %d-1000 (%d snapshots): %+.4f [%+.4f, %+.4f] per 100 generations (cluster bootstrap, %d, seed %d); IQR of replicate slopes [%+.4f, %+.4f]\n', ...
        cfg.from, sum(g >= cfg.from), m, ci, cfg.nBoot, cfg.bootSeed, prctile(sl, [25 75]));
    [m2, ci2] = boot_mean(rep_slopes(g, share, 500, 700), cfg.nBoot, st);
    fprintf('context: mean slope over generations 501-700: %+.4f [%+.4f, %+.4f]\n', m2, ci2);
    % The share of members at generation 1,000, and its mean over the stored generations 700-1,000
    [o1000, ciO] = boot_mean(share(:, end), cfg.nBoot, st);
    out.share1000 = o1000; out.share1000CI = ciO;
    fprintf('share of members at generation 1000 (snapshot, all individuals): %.4f [%.4f, %.4f]; at 500 (stored trace) %.4f\n', ...
        o1000, ciO, mean(vo.occ(:, end)));
    [o7, ci7] = boot_mean(mean(share(:, g >= cfg.from), 2), cfg.nBoot, st);
    fprintf('mean share of members over the 700-1000 snapshots: %.4f [%.4f, %.4f]\n', o7, ci7);
    % The Table S3 sample on the final populations: 30 wirings per replicate after rng(42); the
    % printed generation-500 value beside it is fixed text
    rng(42); nP = 0; nS = 0; perRep = nan(100, 1);
    for r = 1:100
        pick = randperm(size(res(r).finalB, 3), cfg.nSub); k0 = 0;
        for a = pick
            W = double(res(r).finalM(:,:,a)) .* double(res(r).finalB(:,:,a));
            k0 = k0 + membership_test(W, cfg.params).is_css;
        end
        perRep(r) = k0 / cfg.nSub; nP = nP + k0; nS = nS + cfg.nSub;
    end
    [~, ciT] = boot_mean(perRep, cfg.nBoot, st);
    fprintf('Table S3 method on the final populations (rng(42), 30 per replicate): %d/%d = %.4f [%.4f, %.4f] (generation 500: 0.4077)\n', ...
        nP, nS, nP / nS, ciT);
    % success at 1000
    hcEnd = arrayfun(@(x) x.hc(end), res);
    nSucc = sum(hcEnd >= 0.95); [~, pci] = binofit(nSucc, 100);
    out.success = nSucc / 100;
    fprintf('success at generation 1000 (cooperation >= 0.95): %d/100 = %.2f [%.3f, %.3f] (Clopper-Pearson); at generation 500: 58/100\n', ...
        nSucc, nSucc / 100, pci);
    d0 = load(fullfile(cfg.pub, cfg.file), 'results');
    s500 = cellfun(@(x) x.historyCooperation(end) >= 0.95, d0.results(:));
    s1000 = hcEnd(:) >= 0.95;
    fprintf('fate transitions 500 -> 1000: S->S %d, S->F %d, F->S %d, F->F %d\n', ...
        sum(s500 & s1000), sum(s500 & ~s1000), sum(~s500 & s1000), sum(~s500 & ~s1000));
    % Autocorrelation time over 700-981 (spacing 20), as in share_of_members_drift
    Y = share(:, g >= cfg.from & g < 1000); Y = Y - mean(Y, 2);
    a = Y(:, 1:end-1); b = Y(:, 2:end); rho1 = sum(a(:).*b(:)) / sqrt(sum(a(:).^2) * sum(b(:).^2));
    tau = -20 / log(max(rho1, eps));
    fprintf('lag-1 autocorrelation of within-replicate deviation (700-981, spacing 20): %.3f, tau %.1f generations\n', rho1, tau);
    fprintf('\nslope interval %s 0\n', ternary(ci(1) <= 0 && ci(2) >= 0, 'includes', 'excludes'));
    diary off;
end


% =========================================================================
function o = share_full(M, B, params)
%SHARE_FULL  The share of members over every individual, each distinct wiring tested once
%   (share_of_members_by_generation's rule).
    N = size(M, 3);
    Wall = reshape(double(M) .* double(B), 16, N)';
    [Wu, ~, ic] = unique(Wall, 'rows');
    mem = false(size(Wu, 1), 1);
    for u = 1:size(Wu, 1), mem(u) = membership_test(reshape(Wu(u, :), 4, 4), params).is_css; end
    o = mean(mem(ic));
end


function sl = rep_slopes(g, share, from, to)
% Least-squares slope of each replicate's share over generations from..to, per 100 generations
    w = g >= from & g <= to; gw = g(w);
    sl = nan(size(share, 1), 1);
    for i = 1:size(share, 1)
        p = polyfit(gw(:), share(i, w)', 1); sl(i) = 100 * p(1);
    end
end


function [m, ci] = boot_mean(x, nBoot, stream)
    x = x(~isnan(x)); m = mean(x); n = numel(x);
    bm = nan(nBoot, 1);
    for b = 1:nBoot, bm(b) = mean(x(randi(stream, n, n, 1))); end
    ci = prctile(bm, [2.5 97.5]);
end


function open_pool(cfg)
% cfg.nWorkers processes, one thread each; on Windows at below-normal priority (PowerShell), and
% stops if less than cfg.minFreeGB of memory is free (both skipped on other systems).
    if ~isempty(gcp('nocreate')), return; end
    if ispc
        [~, sys] = memory;
        freeGB = sys.PhysicalMemory.Available / 2^30;
        fprintf('available RAM before the pool: %.2f GB (need %.1f)\n', freeGB, cfg.minFreeGB);
        assert(freeGB >= cfg.minFreeGB, 'not enough free RAM for %d workers', cfg.nWorkers);
    end
    p = parpool('Processes', cfg.nWorkers);
    wait(parfevalOnAll(p, @maxNumCompThreads, 0, 1));
    if ispc
        f = parfevalOnAll(p, @() feature('getpid'), 1);
        pids = [feature('getpid'); fetchOutputs(f)];
        cmd = sprintf(['powershell -NoProfile -Command "Get-Process -Id %s | ForEach-Object ', ...
            '{ $_.PriorityClass = ''BelowNormal'' }"'], strjoin(arrayfun(@num2str, pids', 'uni', 0), ','));
        st = system(cmd);
        fprintf('pool: %d workers, 1 thread each; BelowNormal on PIDs %s (status %d)\n', cfg.nWorkers, mat2str(pids'), st);
    else
        fprintf('pool: %d workers, 1 thread each; priority and free memory not changed or checked on this system\n', cfg.nWorkers);
    end
end


function res = run_pool(jobs, ck)
% parfeval over jobs, saving each finished job to the checkpoint file ck (keyed by job.key);
% jobs already in the checkpoint are not rerun.
    key = cellfun(@(j) j.key, jobs, 'uni', 0);
    res = cell(numel(jobs), 1);
    ckKey = {}; ckRes = {};
    if exist(ck, 'file')
        C = load(ck, 'ckKey', 'ckRes'); ckKey = C.ckKey; ckRes = C.ckRes;
        [tf, loc] = ismember(key, ckKey);
        res(tf) = ckRes(loc(tf));
        fprintf('resumed %d of %d jobs from %s\n', sum(tf), numel(jobs), ck);
    end
    todo = find(cellfun(@isempty, res));
    if isempty(todo), return; end
    p = gcp;
    for n = numel(todo):-1:1
        F(n) = parfeval(p, @run_job, 1, jobs{todo(n)}); %#ok<AGROW>
    end
    t0 = tic;
    for n = 1:numel(todo)
        [f, r] = fetchNext(F);
        res{todo(f)} = r;
        [tf, loc] = ismember(r.key, ckKey);
        if tf, ckRes{loc} = r; else, ckKey{end+1} = r.key; ckRes{end+1} = r; end %#ok<AGROW>
        save(ck, 'ckKey', 'ckRes');
        fprintf('job %d/%d done (%s, %.0f s) at %.1f min\n', n, numel(todo), r.key, r.sec, toc(t0) / 60);
    end
end


function diary_on(f)
    if exist(f, 'file'), delete(f); end
    diary(f);
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
