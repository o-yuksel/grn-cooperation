function out = run_repair_test(cfg)
% RUN_REPAIR_TEST  Does restoring retaliation keep failing populations from collapsing? Restores it
% in each failing population at its restart point (the repair), or moves the same individuals as
% far without restoring it (the sham), and restarts both with and without the defector that later
% took over.
%
% Produces : collapse in the repair and sham arms and the checks of the design, for Table S10,
%            Fig. S1, SI S4.3, main text "How cooperating populations lose cooperation" and Methods
%            "Failure route" (read by repair_test_summary, repair_test_without_defector,
%            repair_vs_sham_checks and figS1_repair_test)
% Inputs   : under the data root's publication folder: repair_restart_points/repair_restart_points_events.mat
%              (run_repair_restart_points); repair_test_observed/repair_test_observed.mat
%              (run_repair_test_observed); collapse_census/collapse_census.mat (census_before_failure);
%              <architecture>_100reps.mat (run_signaling_architectures; the earlier snapshots)
% Outputs  : repair_matched/repair_matched_build_targets.mat ('build': targets, arm populations, checks)
%            repair_matched/repair_matched_main.mat ('main': the runs of every arm)
%            repair_matched/repair_matched_main_ckpt.mat (finished jobs; a call with resume true
%              continues from it, and 'main' stops if it exists and resume is false)
% Run      : run_repair_test(struct('mode', 'build'))   (about 4 min on 3 workers), then
%            run_repair_test(struct('mode', 'main'))   (about 85 min on 3 workers; only after a build
%              that passed)
% Seed     : replicate r of event e under rng(1000000 + 1000 e + r, 'twister') (in evolve_population),
%            as in the observed arm; the other seeds are given at target_pool, start_coop and make_job.
%            The seed ranges do not overlap in the paper's run: the sham start's rngSeed + 300000 + e
%            and a restart-point genotype's rngSeed + 300000 + 1000 e + j meet only for an event
%            number e above 1,000 (the runs have 91), and rngSeed + 5000000 + 10000 e + row meets
%            the next event only for a pool above 10,000 genotypes (the largest has 1,139)
%
% The repair replaces every cooperator that does not retaliate by its nearest member (Frobenius
% distance between wirings); the sham moves the same individuals about as far, to other such
% cooperators. Targets come from the restart point, the previous 100 generations and the matched
% population that held cooperation, and must cooperate with their clone in at least 99% of encounters
% under the runs' encounter rule, because the membership test (fixed start, 600 steps) and the runs
% (random start, 300 steps) can disagree. Every arm uses the same ten seeds; 'build' stops on any
% failed check (listed at run_build).
    %% Parameters
    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    assert(any(strcmp(which('evolve_population'), ...
        fullfile(root, 'model', 'evolve_population.m'))), 'evolve_population is not the copy in this release folder');
    assert(cfg.nWorkers <= 3, 'use at most 3 workers (cfg.nWorkers)');
    open_pool(cfg);

    %% Inputs
    E = load(fullfile(cfg.dataDir, 'repair_restart_points', 'repair_restart_points_events.mat'), 'out');
    ev = E.out.events;
    ev = ev([ev.usable]);
    fprintf('events in repair_restart_points_events.mat: %d with a member at the restart point, %d main, %d run-time\n', ...
        numel(ev), sum([ev.inRun]), sum([ev.isPilot]));

    %% Compute
    switch cfg.mode
        case 'build', out = run_build(ev, cfg);
        case 'main',  out = run_arms(ev, load_build(cfg), {'Bp', 'Sp', 'B0p', 'S0p'}, cfg);
        otherwise, error('unknown mode %s', cfg.mode);
    end
    out.cfg = cfg;

    %% Save
    if strcmp(cfg.mode, 'main')                 % run_build saves to build_file(cfg) itself
        save(fullfile(cfg.outDir, 'repair_matched_main.mat'), 'out', '-v7.3');
    end
end


function cfg = fill_defaults(cfg)
    def = struct('dataDir', grn_data_file('publication'), 'outDir', '', 'rngSeed', 20260928, 'mode', 'build', 'k', 10, ...
        'm', 10, 'nGen', 100, 'lambda', 0.005, 'nWorkers', 3, 'tau', 0.99, 'maxExcluded', 20, ...
        'collapseFloor', 0.5, 'nCheck', 5, 'minFreeGB', [], 'seedOffset', 0, 'poolGens', 100, 'resume', false);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    % The target rule of the paper's runs, recorded in the build: the wide pool (restart point,
    % earlier snapshots, matched control) and a sham from the repair target's source
    cfg.pool = 'wide'; cfg.matchSource = true;
    if isempty(cfg.minFreeGB), cfg.minFreeGB = 1.2 * cfg.nWorkers; end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'repair_matched'); end
    cfg.injDir = fullfile(cfg.dataDir, 'repair_test_observed');
end


function open_pool(cfg)
% At most cfg.nWorkers processes, one thread each; on Windows opened only when enough memory is
% free, and the client and the workers run at BelowNormal priority (both skipped on other systems).
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
        fprintf('pool: %d workers, 1 thread each; BelowNormal on PIDs %s (status %d)\n', ...
            cfg.nWorkers, mat2str(pids'), st);
    else
        fprintf('pool: %d workers, 1 thread each; priority and free memory not changed or checked on this system\n', cfg.nWorkers);
    end
end


function out = run_build(ev, cfg)
% Targets, exclusions and the unfiltered-target check per event (parallel over events), then the
% distance check and the starting-cooperation check; the build is saved before either can stop.
% Checks, each an error when it fails: an event is left out when one of its cooperators that do not
% retaliate has no qualifying repair or sham target, and the build stops when more than 20 main
% events are left out; without the filter, targets from the restart point must reproduce the stored
% repair and sham of repair_restart_points_events.mat exactly; the sham-to-repair distance ratio per replaced
% individual must have a median within [0.95, 1.05]; over the included main events, the starting
% cooperation (entering generation 1 without the added defector; observed = the stored value) must
% give median(repair - observed) >= -0.02, median(sham - observed) >= -0.02 and
% repair >= observed - 0.05 in at least 90% of events.
    t0 = tic;
    ev = attach_pool(ev, cfg);
    T = cell(numel(ev), 1);
    parfor i = 1:numel(ev)
        T{i} = event_targets(ev(i), cfg);
    end
    T = [T{:}];
    fprintf('targets built for %d events in %.1f min (tau = %.2f, pool %s)\n', ...
        numel(T), toc(t0) / 60, cfg.tau, cfg.pool);
    out = struct('targets', T, 'buildMin', toc(t0) / 60);
    out.exclusion = report_exclusion(T, cfg);
    out.cdist = report_cdist(T);
    out.cstart = report_cstart(T);
    out.reported = report_targets(T);
    out.cfg = cfg;
    save(build_file(cfg), 'out', '-v7.3');
    if out.exclusion.stop
        error('STOP: %d of %d main events excluded (> %d); do not run', out.exclusion.nMain, sum([T.inRun]), cfg.maxExcluded);
    end
    assert(out.cdist.pass, 'C-dist FAIL');
    if ~out.cstart.pass, error('C-start FAIL: do not run'); end
end


function f = build_file(cfg)
    f = fullfile(cfg.outDir, 'repair_matched_build_targets.mat');
end


function ev = attach_pool(ev, cfg)
% Wide pool: the replicate's earlier stored snapshots with sg in [gPre - poolGens, gPre), and
% the snapshot at kPre of the matched population that held cooperation, where it exists. One
% architecture file at a time, cleared after use.
    Cn = load(fullfile(cfg.dataDir, 'collapse_census', 'collapse_census.mat'), 'out'); Cn = Cn.out;
    I = load(fullfile(cfg.injDir, 'repair_test_observed.mat'), 'out'); st = I.out.events;
    [ev.poolM] = deal({}); [ev.poolB] = deal({}); [ev.poolSrc] = deal([]); [ev.poolGen] = deal([]);
    for a = unique([ev.arch])
        d = load(fullfile(cfg.dataDir, Cn(a).file), 'results');
        sg = Cn(a).sg;
        for i = find([ev.arch] == a)
            ev(i) = pool_snapshots(ev(i), d.results, sg, st([st.num] == ev(i).num), cfg);
        end
        clear d
    end
    n = arrayfun(@(e) sum(e.poolSrc == 2), ev); nc = arrayfun(@(e) sum(e.poolSrc == 3), ev);
    fprintf('wide pool: earlier snapshots per event %d to %d (median %g); control snapshot in %d of %d events\n', ...
        min(n), max(n), median(n), sum(nc), numel(ev));
end


function e = pool_snapshots(e, results, sg, s, cfg)
    assert(strcmp(s.id, e.id), 'event bookkeeping (%s)', e.id);
    R = results{e.rep};
    assert(sg(e.kPre) == e.gPre && isequal(R.historyM(:,:,:,e.kPre), e.Mpre) && ...
        isequal(R.historyB(:,:,:,e.kPre), e.Bpre), '%s: restart-point snapshot does not match the architecture file', e.id);
    ks = find(sg >= e.gPre - cfg.poolGens & sg < e.gPre);
    for k = ks(:)'
        e.poolM{end+1} = R.historyM(:,:,:,k); e.poolB{end+1} = R.historyB(:,:,:,k);
        e.poolSrc(end+1) = 2; e.poolGen(end+1) = sg(k);
    end
    if isfinite(s.ctrlRep)
        e.poolM{end+1} = results{s.ctrlRep}.historyM(:,:,:,e.kPre);
        e.poolB{end+1} = results{s.ctrlRep}.historyB(:,:,:,e.kPre);
        e.poolSrc(end+1) = 3; e.poolGen(end+1) = sg(e.kPre);
    end
end


function t = event_targets(e, cfg)
% One event: classes, the unfiltered-target check, filtered targets from the pool, the arms'
% populations and their starting cooperation.
    tic0 = tic;
    [Gu, ig, cnt, Wg] = genotypes({e.Mpre}, {e.Bpre}, size(e.Mpre, 3));
    cls = test_membership(Wg, e.sp);
    mem = find(cls == 1); cooperate_only = find(cls == 2);
    ind = find(cls(ig) == 2);
    assert(isequal(ind, changed(e.Mpre, e.Bpre, e.Mrep, e.Brep)), '%s: cooperate_only set differs from the stored repair', e.id);
    assert(isequal(ind, changed(e.Mpre, e.Bpre, e.Msham, e.Bsham)), '%s: stored sham changed other individuals', e.id);
    % Unfiltered-target check: without the filter, targets from the restart point reproduce the
    % stored repair and sham populations exactly
    [tb0, db0, ts0] = choose(Wg, cnt, cooperate_only, mem, cooperate_only, @(j) true, []);
    [M, B] = place(e, ind, Gu, tb0(ig(ind)));
    okR = isequal(M, e.Mrep) && isequal(B, e.Brep);
    [M, B] = place(e, ind, Gu, ts0(ig(ind)));
    okS = isequal(M, e.Msham) && isequal(B, e.Bsham);
    assert(okR && okS, 'C-code FAIL: %s unfiltered targets differ from the stored repair (%d) or sham (%d)', e.id, okR, okS);
    % filtered targets from the pool; the sham from the repair target's source
    P = target_pool(e, cfg);
    qs = unique(P.ig(ind));
    [tb, db, ts, hasS] = choose(P.Wg, P.cnt, qs, find(P.cls == 1), find(P.cls == 2), @(j) pool_c(P, j, e, cfg) >= cfg.tau, P.src);
    noB = qs(tb(qs) == 0); noS = qs(~hasS(qs));
    pm = find(P.cls == 1 & P.src == 1);
    pm = pm(arrayfun(@(j) isKey(P.memo, j), pm));
    cm = arrayfun(@(j) P.memo(j), pm);
    t = struct('id', e.id, 'num', e.num, 'arch', e.arch, 'inRun', e.inRun, 'isPilot', e.isPilot, ...
        'nC1oGenotypes', numel(qs), 'nMember', numel(mem), 'nMemberPool', sum(P.cls == 1), ...
        'nC1oPool', sum(P.cls == 2), 'nPool', size(P.Gu, 1), 'nCcomputed', P.memo.Count, ...
        'cMemberMedian', median(cm), 'excluded', ~isempty(noB) || ~isempty(noS), ...
        'nNoB', numel(noB), 'nNoS', numel(noS), 'ind', ind, 'nReplaced', numel(ind), 'codeOK', okR && okS, ...
        'Mb', [], 'Bb', [], 'Ms', [], 'Bs', [], 'dB', [], 'dS', [], 'dBstored', db0(ig(ind)), 'ratio', [], ...
        'cB', [], 'cS', [], 'srcB', [], 'srcS', [], 'coop1A', e.coop1A, 'coop1Arepro', NaN, ...
        'coop1Bstored', e.coop1B, 'coop1Bp', NaN, 'coop1Sp', NaN, 'sec', NaN);
    if ~t.excluded
        t = fill_targets(t, e, P, tb, db, ts, cfg);
    end
    t.sec = toc(tic0);
    fprintf('  %s: pool %d genotypes (%d members, %d cooperate_only), c computed for %d, excluded %d, %.0f s\n', ...
        e.id, t.nPool, t.nMemberPool, t.nC1oPool, t.nCcomputed, t.excluded, t.sec);
end


function t = fill_targets(t, e, P, tb, db, ts, cfg)
    ind = t.ind; g = P.ig(ind);
    [t.Mb, t.Bb] = place(e, ind, P.Gu, tb(g));
    [t.Ms, t.Bs] = place(e, ind, P.Gu, ts(g));
    assert(isequal(changed(e.Mpre, e.Bpre, t.Mb, t.Bb), ind) && isequal(changed(e.Mpre, e.Bpre, t.Ms, t.Bs), ind), ...
        'C-dist FAIL: %s repair and sham replace different individuals', e.id);
    t.dB = wdist(e.Mpre, e.Bpre, t.Mb, t.Bb, ind); t.dS = wdist(e.Mpre, e.Bpre, t.Ms, t.Bs, ind);
    assert(max(abs(t.dB - db(g))) < 1e-9, '%s: repair distance bookkeeping', e.id);
    t.ratio = t.dS ./ t.dB;
    t.cB = arrayfun(@(j) pool_c(P, j, e, cfg), tb(g)); t.cS = arrayfun(@(j) pool_c(P, j, e, cfg), ts(g));
    t.srcB = P.src(tb(g)); t.srcS = P.src(ts(g));
    assert(isequal(t.srcB, t.srcS), '%s: sham source differs from repair source', e.id);
    [t.coop1Arepro, t.coop1Bp, t.coop1Sp] = start_coop(e, t, cfg);
end


function [Gu, ig, cnt, Wg] = genotypes(Ms, Bs, N)
% Distinct genotypes over the snapshots in Ms/Bs (counts summed); ig maps the first N
% individuals (the restart point, which comes first). As repair() of run_repair_test_observed for a
% single snapshot.
    G = cell(numel(Ms), 1);
    for k = 1:numel(Ms)
        n = size(Ms{k}, 3);
        G{k} = [reshape(double(Ms{k}), 16, n); reshape(double(Bs{k}), 16, n)]';
    end
    [Gu, ~, igAll] = unique(vertcat(G{:}), 'rows'); cnt = accumarray(igAll, 1);
    ig = igAll(1:N);
    Wg = Gu(:, 1:16) .* Gu(:, 17:32);
end


function cls = test_membership(Wg, sp)
% Class per distinct genotype by the membership test, as in repair() of run_repair_test_observed:
% 1 member, 2 cooperates but does not retaliate, 3 other.
    [Wu, ~, iw] = unique(Wg, 'rows');
    cw = zeros(size(Wu, 1), 1);
    for u = 1:size(Wu, 1)
        rr = membership_test(reshape(Wu(u, :), 4, 4), sp);
        cw(u) = 3 - rr.is_weak - rr.is_css;
    end
    cls = cw(iw);
end


function P = target_pool(e, cfg)
% The target pool: the restart point first (source 1), then the earlier snapshots (2) and the
% matched population's snapshot (3); abundance = counts summed over the pooled snapshots. A genotype's source is
% the lowest source holding it.
% Seeds: a
% restart-point genotype gets rngSeed + 300000 + 1000 num + j (j = its row among the
% restart-point genotypes), any other genotype rngSeed + 5000000 + 10000 num + row.
    N = size(e.Mpre, 3);
    Ms = [{e.Mpre}, e.poolM]; Bs = [{e.Bpre}, e.poolB]; src = [1, e.poolSrc];
    [P.Gu, P.ig, P.cnt, P.Wg] = genotypes(Ms, Bs, N);
    s = cell2mat(arrayfun(@(k) repmat(src(k), size(Ms{k}, 3), 1), (1:numel(Ms))', 'uni', 0));
    Gall = cell2mat(cellfun(@(M, B) [reshape(double(M), 16, []); reshape(double(B), 16, [])]', ...
        Ms(:), Bs(:), 'uni', 0));
    [~, ia] = ismember(Gall, P.Gu, 'rows');
    P.src = accumarray(ia, s, [size(P.Gu, 1), 1], @min);
    P.cls = test_membership(P.Wg, e.sp);
    Grestart = [reshape(double(e.Mpre), 16, N); reshape(double(e.Bpre), 16, N)]';
    [inRestart, jp] = ismember(P.Gu, unique(Grestart, 'rows'), 'rows');
    P.seed = cfg.rngSeed + 5000000 + 10000 * e.num + (1:size(P.Gu, 1))';
    P.seed(inRestart) = cfg.rngSeed + 300000 + 1000 * e.num + jp(inRestart);
    P.memo = containers.Map('KeyType', 'double', 'ValueType', 'double');
end


function c = pool_c(P, j, e, cfg)
% Clone cooperation c of pooled genotype j under the runs' encounter rule: N copies, mean of
% nCheck evaluations, under the genotype's fixed seed; memoised in P.memo (a handle), so computed
% once.
    if isKey(P.memo, j), c = P.memo(j); return; end
    N = size(e.Mpre, 3); g = P.Gu(j, :);
    M = repmat(reshape(logical(g(1:16)), 4, 4), 1, 1, N);
    B = repmat(reshape(single(g(17:32)), 4, 4), 1, 1, N);
    rng(P.seed(j), 'twister');
    c = coop_now(M, B, e.config, cfg.nCheck);
    P.memo(j) = c;
end


function [tb, db, ts, hasS] = choose(Wg, cnt, qs, memC, c1oC, ok, srcM)
% Per genotype q in qs (cooperates, does not retaliate): tb = the repair target, the nearest
% member among memC with ok() (ties to the most abundant, as repair()); ts = the sham target, a
% genotype among c1oC (cooperates, does not retaliate) with a different wiring and ok(), whose
% distance is closest to the tb distance (ties to the most abundant, as sham() of
% run_repair_restart_points); with srcM, only from the source of the repair target. 0 = no target.
% hasS = a qualifying sham candidate exists, decided even when tb = 0.
% Candidates are visited by level of the sort key and ok() is asked only until a level
% holds a qualifying one: identical to filtering first, since ok(j) is fixed per genotype.
    tb = zeros(size(Wg, 1), 1); db = nan(size(tb)); ts = zeros(size(tb)); hasS = false(size(tb));
    for q = qs(:)'
        if ~isempty(memC)
            [tb(q), db(q)] = first_level(memC, sqrt(sum((Wg(memC, :) - Wg(q, :)).^2, 2)), cnt, ok);
        end
        if isempty(c1oC), continue; end
        dc = sqrt(sum((Wg(c1oC, :) - Wg(q, :)).^2, 2));
        keep = dc > 0;
        if ~isempty(srcM) && tb(q) > 0, keep = keep & srcM(c1oC) == srcM(tb(q)); end
        cand = c1oC(keep); dc = dc(keep);
        % Without a repair target, only whether a sham candidate exists
        if isnan(db(q))
            hasS(q) = first_level(cand, dc, cnt, ok) > 0;
        else
            ts(q) = first_level(cand, abs(dc - db(q)), cnt, ok); hasS(q) = ts(q) > 0;
        end
    end
end


function [t, v] = first_level(rows, key, cnt, ok)
% Lowest key level holding a row with ok(); among those, the most abundant (first on ties).
    t = 0; v = NaN;
    for L = unique(key)'
        at = rows(key == L);
        at = at(arrayfun(ok, at));
        if isempty(at), continue; end
        [~, b] = max(cnt(at)); t = at(b); v = L; return;
    end
end


function [M, B] = place(e, ind, Gu, tgt)
% The restart-point population with individual ind(n) set to distinct genotype tgt(n) (as
% repair()).
    M = e.Mpre; B = e.Bpre;
    for n = 1:numel(ind)
        M(:,:,ind(n)) = reshape(logical(Gu(tgt(n), 1:16)), 4, 4);
        B(:,:,ind(n)) = reshape(single(Gu(tgt(n), 17:32)), 4, 4);
    end
end


function [a, b, s] = start_coop(e, t, cfg)
% Starting cooperation: observed (recomputed) and repaired under rng(rngSeed + 200000 + num), the
% seed of the stored repair's value in run_repair_test_observed; sham under
% rng(rngSeed + 300000 + num), the seed of the stored sham's value in run_repair_restart_points.
    rng(cfg.rngSeed + 200000 + e.num, 'twister');
    a = coop_now(e.Mpre, e.Bpre, e.config, cfg.nCheck);
    b = coop_now(t.Mb, t.Bb, e.config, cfg.nCheck);
    rng(cfg.rngSeed + 300000 + e.num, 'twister');
    s = coop_now(t.Ms, t.Bs, e.config, cfg.nCheck);
end


function idx = changed(M0, B0, M, B)
% Individuals whose genotype (M, B) differs from their genotype at the restart point.
    N = size(M0, 3);
    idx = find(any(reshape(M ~= M0, 16, N), 1) | any(reshape(B ~= B0, 16, N), 1));
    idx = idx(:);
end


function d = wdist(M0, B0, M, B, sub)
% Frobenius distance of W = M.*B from the restart-point wiring, per individual in sub.
    W0 = double(M0(:,:,sub)) .* double(B0(:,:,sub));
    W1 = double(M(:,:,sub)) .* double(B(:,:,sub));
    d = sqrt(squeeze(sum(sum((W1 - W0).^2, 1), 2)));
    d = d(:);
end


function x = report_exclusion(T, cfg)
    main = [T.inRun]; exc = [T.excluded];
    noB = [T.nNoB] > 0; noS = [T.nNoS] > 0;
    x = struct('nMain', sum(exc & main), 'idsMain', {{T(exc & main).id}}, ...
        'nPilot', sum(exc & ~main), 'idsPilot', {{T(exc & ~main).id}}, ...
        'nMainNoB', sum(noB & main), 'nMainNoS', sum(noS & main), 'nMainNoSonly', sum(noS & ~noB & main), ...
        'byArch', arrayfun(@(a) sum(exc & main & [T.arch] == a), 1:4), ...
        'nMainArch', arrayfun(@(a) sum(main & [T.arch] == a), 1:4));
    x.stop = x.nMain > cfg.maxExcluded;
    fprintf('\nEXCLUSION (tau = %.2f, pool %s): %d of %d main events excluded: %s\n', cfg.tau, cfg.pool, ...
        x.nMain, sum(main), strjoin(x.idsMain, ', '));
    fprintf('  by architecture (GS GC SS SC): %s of %s\n', mat2str(x.byArch), mat2str(x.nMainArch));
    fprintf('  main events lacking a repair target for some cooperate_only genotype: %d; lacking a sham target: %d (sham only: %d)\n', ...
        x.nMainNoB, x.nMainNoS, x.nMainNoSonly);
    fprintf('  of the five run-time events, excluded: %d (%s)\n', x.nPilot, strjoin(x.idsPilot, ', '));
    for t = T(exc)
        fprintf('    %s: of %d cooperate_only genotypes, %d lack a repair target, %d lack a sham target (pool: %d members, %d cooperate_only)\n', ...
            t.id, t.nC1oGenotypes, t.nNoB, t.nNoS, t.nMemberPool, t.nC1oPool);
    end
    fprintf('  stop rule (> %d main events excluded): %s\n', cfg.maxExcluded, yn(x.stop));
    fprintf('C-code (unfiltered targets reproduce stored repair and sham on all %d events): %s\n', numel(T), pf(all([T.codeOK])));
end


function x = report_cdist(T)
    T = T(~[T.excluded]);
    r = vertcat(T.ratio);
    em = arrayfun(@(t) median(t.ratio), T);
    x = struct('median', median(r), 'iqr', quantile(r, [0.25 0.75]), 'n', numel(r), 'sameInd', true);
    x.pass = x.median >= 0.95 && x.median <= 1.05;
    fprintf(['C-dist: sham/repair distance ratio median %.4f [IQR %.3f, %.3f] over %d individuals (events'' medians %.3f to %.3f); ', ...
        'same individuals: YES (asserted per event): %s\n'], x.median, x.iqr, x.n, min(em), max(em), pf(x.pass));
end


function x = report_cstart(T)
    T = T([T.inRun] & ~[T.excluded]);
    a = [T.coop1A]'; b = [T.coop1Bp]'; s = [T.coop1Sp]';
    x = struct('n', numel(T), 'medBA', median(b - a), 'medSA', median(s - a), ...
        'shareB', mean(b >= a - 0.05), 'nB', sum(b >= a - 0.05), ...
        'medA', median(a), 'medB', median(b), 'medS', median(s), ...
        'nReproA', sum([T.coop1Arepro]' == a), 'medBstored', median([T.coop1Bstored]));
    x.pass = x.medBA >= -0.02 && x.medSA >= -0.02 && x.shareB >= 0.90;
    fprintf('C-start (%d included main events): observed start recomputed = stored on %d of %d events\n', x.n, x.nReproA, x.n);
    fprintf('  medians: observed %.3f, repair %.3f, sham %.3f (stored repair %.3f)\n', x.medA, x.medB, x.medS, x.medBstored);
    fprintf('  median(repair - observed) %+.4f (>= -0.02: %s); median(sham - observed) %+.4f (>= -0.02: %s)\n', ...
        x.medBA, yn(x.medBA >= -0.02), x.medSA, yn(x.medSA >= -0.02));
    fprintf('  repair >= observed - 0.05 in %d of %d events (%.1f%%; >= 90%%: %s)\n', x.nB, x.n, 100 * x.shareB, yn(x.shareB >= 0.90));
    fprintf('C-start: %s\n', pf(x.pass));
end


function x = report_targets(T)
    T = T(~[T.excluded]);
    mB = arrayfun(@(t) mean(t.dB), T); mS = arrayfun(@(t) mean(t.dS), T); m0 = arrayfun(@(t) mean(t.dBstored), T);
    cB = vertcat(T.cB); cS = vertcat(T.cS); sB = vertcat(T.srcB); sS = vertcat(T.srcS);
    x = struct('rB', [median(mB ./ m0), min(mB ./ m0), max(mB ./ m0)], 'rS', [median(mS ./ m0), min(mS ./ m0), max(mS ./ m0)], ...
        'cB', [min(cB), median(cB)], 'cS', [min(cS), median(cS)], ...
        'srcB', arrayfun(@(k) mean(sB == k), 1:3), 'srcS', arrayfun(@(k) mean(sS == k), 1:3));
    fprintf('Also reported, %d included events (main and the five run-time events):\n', numel(T));
    fprintf('  per-event mean distance / stored repair mean distance: repair median %.3f [%.3f, %.3f], sham median %.3f [%.3f, %.3f]\n', x.rB, x.rS);
    fprintf('  repair distance equals stored repair distance in %d of %d events\n', sum(abs(mB - m0) < 1e-12), numel(T));
    fprintf('  individuals replaced per event: median %g [%d, %d]\n', median([T.nReplaced]), min([T.nReplaced]), max([T.nReplaced]));
    fprintf('  c of chosen targets (per replaced individual): repair min %.4f median %.4f; sham min %.4f median %.4f\n', x.cB, x.cS);
    fprintf('  target source shares (restart point / earlier snapshot / control), per replaced individual: repair %.3f / %.3f / %.3f; sham %.3f / %.3f / %.3f\n', ...
        x.srcB, x.srcS);
    fprintf('  pool size per event: median %d genotypes [%d, %d]; c computed per event: median %d\n', ...
        median([T.nPool]), min([T.nPool]), max([T.nPool]), median([T.nCcomputed]));
end


function B = load_build(cfg)
% The call 'main' reads the build and refuses to run unless the stop rule did not fire and the
% distance and starting-cooperation checks passed.
    S = load(build_file(cfg), 'out'); B = S.out;
    assert(B.cfg.tau == cfg.tau && strcmp(B.cfg.pool, cfg.pool) && B.cfg.matchSource == cfg.matchSource && B.cfg.rngSeed == cfg.rngSeed, ...
        'build made with another tau, pool or seed');
    assert(~B.exclusion.stop && B.cdist.pass && B.cstart.pass, 'build gates did not pass: do not run');
end


function job = make_job(e, M, B, arm, cfg, k)
% Job construction of run_repair_test_observed (make_jobs) for one event and one population: the
% added copies at slots idx from rng(rngSeed + num), left out in the arms without them (B0p,
% S0p); k replicates, replicate r under the run seed seedBase + seedOffset + r.
% run_timed draws each job's starting cooperation under rng(rngSeed + 600000 + 10 num + arm code,
% 'twister').
    N = e.config.populationSize;
    rng(cfg.rngSeed + e.num);
    idx = randperm(N, cfg.m);
    seedBase = 1000000 + 1000 * e.num;
    invG = e.invG; M0 = M; B0 = B;
    if any(strcmp(arm, {'B0p', 'S0p'})), invG = []; end
    if ~isempty(invG)
        M(:,:,idx) = repmat(reshape(logical(invG(1:16)), 4, 4), 1, 1, cfg.m);
        B(:,:,idx) = repmat(reshape(single(invG(17:32)), 4, 4), 1, 1, cfg.m);
    end
    c = e.config;
    c.numGenerations = cfg.nGen; c.numReplicates = k;
    c.historySampleInterval = 10 * cfg.nGen;
    c.alwaysSaveGenerations = [1 cfg.nGen];
    c.mutationRate = cfg.lambda;
    c.initPopulationM = M; c.initPopulationB = B; c.rngSeedBase = seedBase + cfg.seedOffset;
    job = struct('event', e.num, 'arm', arm, 'idx', idx, 'invG', invG, 'config', c, ...
        'M0', M0, 'B0', B0, 'coopSeed', cfg.rngSeed + 600000 + 10 * e.num + arm_code(arm));
end


function c = arm_code(arm)
    % Codes 5 to 8 are the arms run here; the first four keep the stored runs' seeds
    c = find(strcmp(arm, {'A', 'B', 'S', 'A0', 'Bp', 'Sp', 'B0p', 'S0p'}));
end


function r = run_timed(job, cfg)
% run_job plus worker time and the starting cooperation, under the runs' encounter rule, of the
% job's population without the added copies.
    t = tic;
    r = run_job(job, cfg);
    rng(job.coopSeed, 'twister');
    r.coop1 = coop_now(job.M0, job.B0, job.config, cfg.nCheck);
    r.sec = toc(t);
end


function r = run_job(job, cfg)
% As run_job of run_repair_test_observed: per run the cooperation trace, collapse (cooperation below
% collapseFloor at any generation) and the share descended from the slots of the added copies.
    R = evolve_population(job.config);
    k = numel(R); nG = cfg.nGen;
    r = struct('hc', zeros(k, nG, 'single'), 'collapse', false(k, 1), 'firstBelow', nan(k, 1), ...
        'coopEnd', nan(k, 1), 'minCoop', nan(k, 1), 'lineage', nan(k, 1), 'invFreq', nan(k, 1), ...
        'seedUsed', false(k, 1));
    M0 = job.config.initPopulationM; B0 = job.config.initPopulationB;
    for q = 1:k
        x = R{q}; hc = x.historyCooperation(:)';
        r.hc(q, :) = hc;
        fb = find(hc < cfg.collapseFloor, 1, 'first');
        r.collapse(q) = ~isempty(fb);
        if ~isempty(fb), r.firstBelow(q) = fb; end
        r.coopEnd(q) = hc(end); r.minCoop(q) = min(hc);
        flag = false(size(x.historyParents, 1), 1); flag(job.idx) = true;
        for g = 1:nG - 1
            flag = flag(double(x.historyParents(:, g)));
        end
        r.lineage(q) = mean(flag);
        r.seedUsed(q) = isequal(x.historyM(:,:,:,1), logical(M0)) && isequal(x.historyB(:,:,:,1), single(B0));
        if ~isempty(job.invG)
            Mg = x.historyM(:,:,:,end); Bg = x.historyB(:,:,:,end); N = size(Mg, 3);
            Gg = [reshape(double(Mg), 16, N); reshape(double(Bg), 16, N)]';
            r.invFreq(q) = mean(all(Gg == double(single(job.invG)), 2));
        end
    end
end


function out = run_arms(ev, Bd, arms, cfg)
% Arms Bp Sp B0p S0p, k seeds, for the included main events.
% Job labels: 'Bp' the repair and 'B0p' the same without the added defector; 'Sp' and 'S0p' the
% sham, another cooperator that does not retaliate, from the source of the repair target, at the
% distance closest to the repair's (ties to the most abundant). Arms A and A0 are not rerun.
    T = Bd.targets;
    T = T([T.inRun] & ~[T.excluded]);
    jobs = {};
    for t = T
        e = ev([ev.num] == t.num);
        assert(strcmp(e.id, t.id), 'event bookkeeping');
        for a = 1:numel(arms)
            if any(strcmp(arms{a}, {'Bp', 'B0p'})), M = t.Mb; B = t.Bb; else, M = t.Ms; B = t.Bs; end
            jobs{end+1} = make_job(e, M, B, arms{a}, cfg, cfg.k); %#ok<AGROW>
        end
    end
    fprintf('main: %d events, %d jobs, k = %d\n', numel(T), numel(jobs), cfg.k);
    ck = fullfile(cfg.outDir, 'repair_matched_main_ckpt.mat');
    if exist(ck, 'file') && ~cfg.resume
        error('run_repair_test:checkpoint', '%s exists: pass resume true to continue from it, or delete it for a fresh run', ck);
    end
    t1 = tic;
    res = run_pool(jobs, cfg, ck);
    out = struct('runMin', toc(t1) / 60, 'nJobs', numel(jobs), 'events', {{T.id}}, 'arms', {arms});
    out.runs = struct('event', {}, 'arm', {}, 'r', {});
    for j = 1:numel(jobs)
        out.runs(end+1) = struct('event', jobs{j}.event, 'arm', jobs{j}.arm, 'r', res{j});
    end
    sec = cellfun(@(r) r.sec, res);
    out.secPerRun = sum(sec) / (numel(jobs) * cfg.k);
    fprintf('main done in %.1f min; %.2f worker-s per run\n', out.runMin, out.secPerRun);
end


function res = run_pool(jobs, cfg, ck)
% parfeval over jobs; every finished job is saved to the checkpoint file, and a call with resume
% true continues from it (jobs are keyed by event and arm).
    key = cellfun(@(j) sprintf('%d-%s', j.event, j.arm), jobs, 'uni', 0);
    res = cell(numel(jobs), 1);
    if cfg.resume && exist(ck, 'file')
        C = load(ck, 'ckKey', 'ckRes');
        [tf, loc] = ismember(key, C.ckKey);
        res(tf) = C.ckRes(loc(tf));
        fprintf('resumed %d of %d jobs from %s\n', sum(tf), numel(jobs), ck);
    end
    todo = find(cellfun(@isempty, res));
    p = gcp;
    for n = numel(todo):-1:1
        F(n) = parfeval(p, @run_timed, 1, jobs{todo(n)}, cfg); %#ok<AGROW>
    end
    t0 = tic;
    for n = 1:numel(todo)
        [f, r] = fetchNext(F);
        res{todo(f)} = r;
        ckKey = key(~cellfun(@isempty, res)); ckRes = res(~cellfun(@isempty, res)); %#ok<NASGU>
        save(ck, 'ckKey', 'ckRes');
        fprintf('job %d/%d done (%s, %.0f s) at %.1f min\n', n, numel(todo), key{todo(f)}, r.sec, toc(t0) / 60);
    end
end


function s = yn(t), if t, s = 'YES'; else, s = 'NO'; end, end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
