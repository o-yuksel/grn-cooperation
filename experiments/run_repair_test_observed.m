function out = run_repair_test_observed(cfg)
% RUN_REPAIR_TEST_OBSERVED  How often do failing populations collapse when restarted as they were?
% The observed arm of the repair test: each failing population of the census before failure
% restarted from its restart point (the last stored generation at or above 0.90 before the fall),
% with and without the defector that later took over.
%
% Produces : the observed arm and the populations the repair test uses, for Table S10, Fig. S1,
%            SI S4.3 and main text "How cooperating populations lose cooperation" (read by
%            run_repair_restart_points, run_repair_test, repair_test_summary,
%            repair_test_observed_arm, repair_test_without_defector and figS1_repair_test)
% Inputs   : under the data root's publication folder: collapse_census/collapse_census.mat
%              (census_before_failure); <architecture>_100reps.mat for the four architectures
%              (run_signaling_architectures)
% Outputs  : repair_test_observed/repair_test_observed.mat   (the events and the runs of arms A and A0)
% Run      : run_repair_test_observed()   (not timed)
% Seed     : replicate r of event e under rng(1000000 + 1000 e + r, 'twister') (in evolve_population),
%            the same in both arms; default generator: slots of the added copies rng(20260928 + e),
%            the random choice of the five events kept out of the runs rng(20260928), starting cooperation
%            rng(20260928 + 200000 + e)
%            on a pool worker
%
% Arm A restarts each population as observed, with ten added copies of the defector, and A0 the
% same without them (make_jobs).

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    assert(any(strcmp(which('evolve_population'), ...
        fullfile(root, 'model', 'evolve_population.m'))), 'evolve_population not from this tree');

    t0 = tic;
    ev = build_events(cfg);
    fprintf('events built in %.1f min\n', toc(t0) / 60);
    [ev, pilotIds] = select_events(ev, cfg);
    jobs = make_jobs(ev, cfg);
    fprintf('main: %d events, %d jobs, k = %d, nGen = %d\n', sum([ev.inRun]), numel(jobs), cfg.k, cfg.nGen);
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
    res = cell(numel(jobs), 1);
    t1 = tic;
    parfor j = 1:numel(jobs)
        res{j} = run_job(jobs{j}, cfg);
    end
    runMin = toc(t1) / 60;
    fprintf('jobs done in %.1f min\n', runMin);
    out = struct('cfg', cfg, 'events', rmfield(ev, {'Mpre', 'Bpre', 'Mrep', 'Brep'}), ...
        'pilotIds', {pilotIds}, 'runMin', runMin, 'nJobs', numel(jobs));
    out.runs = assemble(res, jobs);
    save(fullfile(cfg.outDir, 'repair_test_observed.mat'), 'out', '-v7.3');
end


function cfg = fill_defaults(cfg)
    def = struct('dataDir', grn_data_file('publication'), 'outDir', '', 'censusFile', '', 'rngSeed', 20260928, ...
        'nPilot', 5, 'k', 10, 'm', 10, 'nGen', 100, 'lambda', 0.005, ...
        'arms', {{'A', 'A0'}}, 'nWorkers', 3, ...
        'collapseFloor', 0.5, 'maxEvents', Inf, 'nCheck', 5);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'repair_test_observed'); end
    if isempty(cfg.censusFile)
        cfg.censusFile = fullfile(cfg.dataDir, 'collapse_census', 'collapse_census.mat');
    end
end


function ev = build_events(cfg)
% One record per failing population of the census with a defector that later took over: the
% population at its restart point, its classes by the membership test and its repair.
% The defector is the commonest wiring that defects after the fall; a population with no member at
% its restart point has nothing to repair toward and is skipped (select_events then leaves out
% five drawn at random, used to measure run time). The starting cooperation stored here is what
% run_repair_test builds on.
    C = load(cfg.censusFile, 'out'); C = C.out;
    CODE = {'GS', 'GC', 'SS', 'SC'};
    ev = struct([]);
    for a = 1:numel(C)
        d = load(fullfile(cfg.dataDir, C(a).file), 'results', 'config');
        for q = 1:numel(C(a).coll)
            inv = C(a).coll(q).inv;
            if ~inv.found, continue; end
            e = C(a).events(C(a).collEvent(q));
            R = d.results{e.rep};
            x = struct('id', sprintf('%s-r%03d', CODE{a}, e.rep), 'arch', a, 'archName', C(a).name, ...
                'rep', e.rep, 'kPre', e.kPre, 'gPre', C(a).sg(e.kPre), 'invG', inv.G, ...
                'onGridCls', inv.onGridCls, 'invPresentPre', inv.presentPre, 'usedPost2', C(a).coll(q).usedPost2, ...
                'ctrlRep', e.ctrlRep, 'config', d.config, 'sp', C(a).sp, ...
                'Mpre', R.historyM(:,:,:,e.kPre), 'Bpre', R.historyB(:,:,:,e.kPre), ...
                'Mrep', [], 'Brep', [], 'shares', nan(1, 3), ...
                'nReplaced', NaN, 'usable', false, 'skipReason', '', 'num', 0, 'inRun', false, ...
                'isPilot', false, 'meanRepairDist', NaN);
            ev = [ev, x]; %#ok<AGROW>
        end
        clear d
    end
    for i = 1:numel(ev), ev(i).num = i; end
    rep = cell(numel(ev), 1);
    parfor_pool(cfg);
    parfor i = 1:numel(ev)
        rep{i} = repair(ev(i).Mpre, ev(i).Bpre, ev(i).sp);
    end
    for i = 1:numel(ev)
        r = rep{i};
        ev(i).shares = r.shares; ev(i).nReplaced = r.nReplaced; ev(i).meanRepairDist = r.meanDist;
        if r.nMember == 0
            ev(i).skipReason = 'no member at the restart point';
        else
            ev(i).usable = true; ev(i).Mrep = r.M; ev(i).Brep = r.B;
        end
    end
    % Starting cooperation: cooperation under the runs' encounter rule of the population at its
    % restart point and of the repaired one, without the added defector, mean of nCheck
    % one-generation evaluations; no evolution
    [c1A, c1B] = deal(nan(numel(ev), 1));
    parfor i = 1:numel(ev)
        if ~ev(i).usable, continue; end
        rng(cfg.rngSeed + 200000 + i);
        c1A(i) = coop_now(ev(i).Mpre, ev(i).Bpre, ev(i).config, cfg.nCheck);
        c1B(i) = coop_now(ev(i).Mrep, ev(i).Brep, ev(i).config, cfg.nCheck);
    end
    for i = 1:numel(ev), ev(i).coop1A = c1A(i); ev(i).coop1B = c1B(i); end
end


function parfor_pool(cfg)
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
end


function [ev, pilotIds] = select_events(ev, cfg)
    usable = find([ev.usable]);
    rng(cfg.rngSeed);
    pick = usable(randperm(numel(usable), cfg.nPilot));
    for i = pick, ev(i).isPilot = true; end
    pilotIds = {ev(pick).id};
    run = setdiff(usable, pick);
    run = run(1:min(numel(run), cfg.maxEvents));
    for i = run, ev(i).inRun = true; end
end


function jobs = make_jobs(ev, cfg)
% One job per event and arm. The added copies take the same slots idx in both arms of an event,
% and replicate r of event e runs under the run seed 1000000 + 1000 e + r in both arms.
% Arms: A as observed, with ten copies of the defector, and A0 the same without them. Each
% restarts k = 10 times for 100 generations. A run has collapsed once population cooperation falls
% below 0.5 at any of its generations, a rule that differs from the one that selected the
% populations.
    jobs = {};
    N = ev(1).config.populationSize;
    for i = find([ev.inRun])
        e = ev(i);
        rng(cfg.rngSeed + e.num);
        idx = randperm(N, cfg.m);
        seedBase = 1000000 + 1000 * e.num;
        for a = 1:numel(cfg.arms)
            arm = cfg.arms{a};
            invG = e.invG; M = e.Mpre; B = e.Bpre;
            switch arm
                case 'A'
                case 'A0', invG = [];
                otherwise, error('unknown arm %s', arm);
            end
            if ~isempty(invG)
                M(:,:,idx) = repmat(reshape(logical(invG(1:16)), 4, 4), 1, 1, cfg.m);
                B(:,:,idx) = repmat(reshape(single(invG(17:32)), 4, 4), 1, 1, cfg.m);
            end
            c = e.config;
            c.numGenerations = cfg.nGen; c.numReplicates = cfg.k;
            c.historySampleInterval = 10 * cfg.nGen;
            c.alwaysSaveGenerations = [1 cfg.nGen];
            c.mutationRate = cfg.lambda;
            c.initPopulationM = M; c.initPopulationB = B; c.rngSeedBase = seedBase;
            jobs{end+1} = struct('event', e.num, 'arm', arm, 'idx', idx, 'invG', invG, ...
                'config', c); %#ok<AGROW>
        end
    end
end


function r = run_job(job, cfg)
% Per run: the cooperation trace, collapse (cooperation below collapseFloor at any generation),
% and the share of the last generation descended from the slots of the added copies.
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


function runs = assemble(res, jobs)
    runs = struct('event', {}, 'arm', {}, 'r', {});
    for j = 1:numel(jobs)
        runs(end+1) = struct('event', jobs{j}.event, 'arm', jobs{j}.arm, 'r', res{j}); %#ok<AGROW>
    end
end
