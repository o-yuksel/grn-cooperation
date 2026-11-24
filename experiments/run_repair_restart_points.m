function out = run_repair_restart_points(cfg)
% RUN_REPAIR_RESTART_POINTS  Rebuilds the populations of the repair test at their restart points,
% with a repair and a sham chosen by the membership test alone.
%
% Produces : a repair and sham not filtered on clone cooperation, for Table S10 and SI S4.3 (read
%            by run_repair_test, which builds on and checks against them, and by
%            repair_vs_sham_checks)
% Inputs   : under the data root's publication folder: collapse_census/collapse_census.mat
%              (census_before_failure); <architecture>_100reps.mat for the four architectures
%              (run_signaling_architectures); repair_test_observed/repair_test_observed.mat
%              (run_repair_test_observed)
% Outputs  : repair_restart_points/repair_restart_points_events.mat ('events': the events with their repair and sham)
% Run      : run_repair_restart_points()   (about 1 min on 2 workers)
% Seed     : the fallback sham's permutation from rng(20260928 + 400000 + e, 'twister'); the sham's
%            starting cooperation from rng(20260928 + 300000 + e) on a pool worker
%
% The sham moves each cooperator that does not retaliate to another such wiring at about the
% repair's distance, keeping the class shares (details at sham()).

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    assert(any(strcmp(which('evolve_population'), ...
        fullfile(root, 'model', 'evolve_population.m'))), 'evolve_population not from this tree');
    % cap the pool before any parfor opens a default-size one
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
    evFile = fullfile(cfg.outDir, 'repair_restart_points_events.mat');
    switch cfg.mode
        case 'events'
            out = make_events(cfg); save(evFile, 'out', '-v7.3');
        otherwise
            error('unknown mode %s', cfg.mode);
    end
end


function cfg = fill_defaults(cfg)
    def = struct('dataDir', grn_data_file('publication'), 'outDir', '', 'rngSeed', 20260928, 'mode', 'events', ...
        'nWorkers', 2);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'repair_restart_points'); end
    cfg.injDir = fullfile(cfg.dataDir, 'repair_test_observed');
    cfg.censusFile = fullfile(cfg.dataDir, 'collapse_census', 'collapse_census.mat');
end


function out = make_events(cfg)
% Events: build_events, repair and the job construction are those of run_repair_test_observed.
% Every rebuilt event must equal its stored event in repair_test_observed.mat (class shares,
% number replaced, mean repair distance, the defector that later took over), and the stored event
% gives which events are main or left out and the starting cooperation.
    t0 = tic;
    ev = build_events(cfg);
    S = load(fullfile(cfg.injDir, 'repair_test_observed.mat'), 'out'); st = S.out.events;
    assert(numel(st) == numel(ev), 'event count differs from the stored run');
    for i = 1:numel(ev)
        assert(strcmp(ev(i).id, st(i).id) && isequaln(ev(i).shares, st(i).shares) && ...
            isequaln(ev(i).nReplaced, st(i).nReplaced) && isequaln(ev(i).meanRepairDist, st(i).meanRepairDist) && ...
            isequaln(ev(i).invG, st(i).invG), 'rebuilt event %s differs from the stored one', ev(i).id);
        ev(i).inRun = st(i).inRun; ev(i).isPilot = st(i).isPilot;
        ev(i).coop1A = st(i).coop1A; ev(i).coop1B = st(i).coop1B;
    end
    fprintf('all %d rebuilt events match the stored events exactly\n', numel(ev));
    shm = cell(numel(ev), 1);
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
    parfor i = 1:numel(ev)
        if ~ev(i).usable, continue; end
        shm{i} = sham(ev(i).Mpre, ev(i).Bpre, ev(i).sp, ev(i).nReplaced, cfg.rngSeed + 400000 + i);
    end
    c1S = nan(numel(ev), 1);
    parfor i = 1:numel(ev)
        if ~ev(i).usable, continue; end
        rng(cfg.rngSeed + 300000 + i);
        c1S(i) = coop_now(shm{i}.M, shm{i}.B, ev(i).config, 5);
    end
    for i = 1:numel(ev)
        ev(i).sham = shm{i}; ev(i).coop1S = c1S(i);
        if ~isempty(shm{i}), ev(i).Msham = shm{i}.M; ev(i).Bsham = shm{i}.B; end
        if ~isempty(ev(i).sham), ev(i).sham = rmfield(ev(i).sham, {'M', 'B'}); end
    end
    out = struct('cfg', cfg, 'events', ev, 'buildMin', toc(t0) / 60);
end


function ev = build_events(cfg)
% build_events of run_repair_test_observed without the starting cooperation and the defector of
% another event; make_events takes the starting cooperation from the stored run.
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
                'config', d.config, 'sp', C(a).sp, ...
                'Mpre', R.historyM(:,:,:,e.kPre), 'Bpre', R.historyB(:,:,:,e.kPre), ...
                'Mrep', [], 'Brep', [], 'shares', nan(1, 3), 'nReplaced', NaN, 'usable', false, ...
                'num', 0, 'inRun', false, 'isPilot', false, 'meanRepairDist', NaN, ...
                'coop1A', NaN, 'coop1B', NaN, 'coop1S', NaN, 'sham', [], 'Msham', [], 'Bsham', []);
            ev = [ev, x]; %#ok<AGROW>
        end
        clear d
    end
    for i = 1:numel(ev), ev(i).num = i; end
    rep = cell(numel(ev), 1);
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
    parfor i = 1:numel(ev)
        rep{i} = repair(ev(i).Mpre, ev(i).Bpre, ev(i).sp);
    end
    for i = 1:numel(ev)
        r = rep{i};
        ev(i).shares = r.shares; ev(i).nReplaced = r.nReplaced; ev(i).meanRepairDist = r.meanDist;
        if r.nMember > 0
            ev(i).usable = true; ev(i).Mrep = r.M; ev(i).Brep = r.B;
        end
    end
end


function s = sham(M, B, sp, nRep, seed)
% The sham: the distinct genotypes at the restart point and their classes by the membership test
% are those of repair(). For each genotype q that cooperates but does not retaliate, the repair
% moved its individuals a Frobenius wiring distance d_q, to the nearest member; the sham moves every
% individual of q to the genotype of the same class, with a different wiring, whose distance from
% q's wiring W_q is closest to d_q (ties to the most abundant), so the class shares stay unchanged.
% When that class holds a single distinct wiring, the fallback sham moves the same number of
% individuals (at most as many as exist), drawn by a permutation under rng(seed, 'twister') from
% the members and the other class, each to its nearest member of a different wiring. Unlike the
% sham of run_repair_test, these targets are not filtered on clone cooperation.
    N = size(M, 3);
    G = [reshape(double(M), 16, N); reshape(double(B), 16, N)]';
    [Gu, ~, ig] = unique(G, 'rows'); cnt = accumarray(ig, 1);
    Wg = Gu(:, 1:16) .* Gu(:, 17:32);
    [Wu, ~, iw] = unique(Wg, 'rows');
    cw = zeros(size(Wu, 1), 1);
    for u = 1:size(Wu, 1)
        rr = membership_test(reshape(Wu(u, :), 4, 4), sp);
        cw(u) = 3 - rr.is_weak - rr.is_css;
    end
    cls = cw(iw); mem = find(cls == 1); cooperate_only = find(cls == 2);
    s = struct('type', 'S', 'M', M, 'B', B, 'nReplaced', 0, 'ratio', [], 'dRepair', [], ...
        'dSham', [], 'nC1oWirings', numel(unique(iw(cooperate_only))), 'shares', [], 'sharesOut', []);
    s.shares = [sum(cnt(cls == 1)), sum(cnt(cls == 2)), sum(cnt(cls == 3))] / N;
    target = zeros(size(Gu, 1), 1); dR = nan(size(target)); dS = nan(size(target));
    for q = cooperate_only'
        dd = sqrt(sum((Wg(mem, :) - Wg(q, :)).^2, 2)); dR(q) = min(dd);
        dc = sqrt(sum((Wg(cooperate_only, :) - Wg(q, :)).^2, 2));
        cand = find(dc > 0);
        if isempty(cand), continue; end
        gap = abs(dc(cand) - dR(q));
        best = cand(gap == min(gap));
        [~, b] = max(cnt(cooperate_only(best)));
        target(q) = cooperate_only(best(b)); dS(q) = dc(best(b));
    end
    if s.nC1oWirings >= 2
        ind = find(cls(ig) == 2);
    else                                        % second sham, S2
        s.type = 'S2';
        pool = find(cls(ig) ~= 2);
        rng(seed, 'twister');
        pool = pool(randperm(numel(pool)));
        ind = pool(1:min(nRep, numel(pool)));
        for q = unique(ig(ind))'
            dd = sqrt(sum((Wg(mem, :) - Wg(q, :)).^2, 2));
            ok = find(dd > 0);
            if isempty(ok), continue; end
            best = ok(dd(ok) == min(dd(ok)));
            [~, b] = max(cnt(mem(best)));
            target(q) = mem(best(b)); dS(q) = min(dd(ok)); dR(q) = NaN;
        end
    end
    ind = ind(target(ig(ind)) > 0);
    for i = ind'
        t = target(ig(i));
        s.M(:,:,i) = reshape(logical(Gu(t, 1:16)), 4, 4);
        s.B(:,:,i) = reshape(single(Gu(t, 17:32)), 4, 4);
    end
    s.nReplaced = numel(ind);
    s.dRepair = dR(ig(ind)); s.dSham = dS(ig(ind)); s.ratio = s.dSham ./ s.dRepair;
    G2 = [reshape(double(s.M), 16, N); reshape(double(s.B), 16, N)]';
    [~, ia] = ismember(G2, Gu, 'rows');
    s.sharesOut = [sum(cls(ia) == 1), sum(cls(ia) == 2), sum(cls(ia) == 3)] / N;
end
