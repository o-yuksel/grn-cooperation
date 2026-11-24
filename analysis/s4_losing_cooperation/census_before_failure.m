function out = census_before_failure(cfg)
% CENSUS_BEFORE_FAILURE  Just before cooperation falls, do failing populations hold fewer members
% and more cooperators that do not retaliate than populations that kept cooperation, and how do
% their cooperators answer the defector that later took over?
%
% Produces : main text, "How cooperating populations lose cooperation"; Methods, "Failure route";
%            SI S4.2 and Table S9, bottom rows (failed after 0.90, events, matched, the shares
%            before the fall); prints them itself
% Inputs   : cfg (optional struct; defaults in fill_defaults)
%            publication/<architecture>_100reps.mat   (the four architectures: cooperation
%              history, parents, stored snapshots; written by run_signaling_architectures)
% Outputs  : publication/collapse_census/collapse_census.mat   (per-population results, read by
%              defect_back_rates and the repair-test runs) and collapse_census.txt (the report);
%              a full run overwrites both
% Run      : census_before_failure()   (6 workers; not timed);
%            census_before_failure(struct('reportOnly', true)) reprints the report (about 1 min)
% Seed     : rng(20260927) (cfg.rngSeed, the session's generator type) for the control draws;
%            job j reseeds with 20260927 + j, the shuffled pairing of failing population i of
%            architecture a with 20260927 + 100000 a + i, the report's bootstraps with
%            20260927 + 7; membership_test seeds each wiring itself
%
% A population fails after 0.90 if its cooperation exceeded 0.90 and ends below 0.95. It is read
% at the restart point, the last stored generation at or above 0.90 before the fall, against
% matched controls at that generation. The defector that later took over is the most common
% wiring at the first stored generation after the fall (else the next) that cooperates in under
% half of its encounters with that population under the runs' encounter rule.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
    matFile = fullfile(cfg.outDir, 'collapse_census.mat');
    if cfg.reportOnly
        S = load(matFile, 'out', 'cfg'); out = S.out; cfgRun = S.cfg;
    else
        out = run_all(cfg); cfgRun = cfg;
        save(matFile, 'out', 'cfg');
    end
    txt = fullfile(cfg.outDir, 'collapse_census.txt');
    if exist(txt, 'file'), delete(txt); end
    diary(txt);
    print_report(out, cfgRun);
    diary off;
end


% =========================================================================
function cfg = fill_defaults(cfg)
    pub = fullfile(grn_data_file('publication'));
    def = struct('dataDir', pub, 'outDir', '', 'rngSeed', 20260927, 'engineDraws', 10, ...
        'exploitDraws', 50, 'nCtrlClass', 20, 'nShuffle', 5, 'nBoot', 2000, ...
        'nWorkers', 6, 'arrive', 0.90, 'successThreshold', 0.95, 'collapseFloor', 0.50, ...
        'maxReps', Inf, 'reportOnly', false);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'collapse_census'); end
end


% =========================================================================
function out = run_all(cfg)
    ARCH = {'graded separate',     'graded_separate_100reps.mat'
            'graded contrast',     'graded_contrast_100reps.mat'
            'saturating separate', 'saturating_separate_100reps.mat'
            'saturating contrast', 'saturating_contrast_100reps.mat'};
    rng(cfg.rngSeed);
    [jobs, meta] = build_jobs(ARCH, cfg);
    fprintf('%d jobs\n', numel(jobs));
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
    res = cell(numel(jobs), 1);
    t0 = tic;
    parfor j = 1:numel(jobs)
        rng(cfg.rngSeed + j);
        res{j} = run_job(jobs{j}, cfg);
    end
    fprintf('jobs done in %.1f min\n', toc(t0) / 60);
    out = assemble(res, jobs, meta);
    out = shuffle_baseline(out, cfg);
end


% =========================================================================
function [jobs, meta] = build_jobs(ARCH, cfg)
    jobs = {}; meta = struct([]);
    for a = 1:size(ARCH, 1)
        f = fullfile(cfg.dataDir, ARCH{a, 2});
        d = load(f, 'results', 'config'); info = dir(f); c = d.config;
        nR = min(numel(d.results), cfg.maxReps); R = d.results(1:nR);
        hc = cell2mat(cellfun(@(r) double(r.historyCooperation(:))', R, 'uni', 0));
        sg = double(R{1}.historyGenerations(:))';
        [sp, spE] = sim_params(c);
        m = test_membership_reps(hc, sg, cfg);
        m.name = ARCH{a, 1}; m.file = ARCH{a, 2}; m.fileDate = info.date;
        m.sigma = c.mutationStdWeight; m.engineSteps = c.developmentSteps;
        m.nOpp = c.numOpponents; m.payoffMode = c.payoffMode; m.sg = sg; m.hc = hc;
        m.sp = sp; m.spE = spE;
        m = pick_controls(m, cfg);
        jobs = [jobs, make_jobs(a, m, R, sp, spE)]; %#ok<AGROW>
        meta = [meta, m]; %#ok<AGROW>
    end
end


% =========================================================================
function [sp, spE] = sim_params(c)
    sp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, ...
        'discreteOutput', c.discreteOutput, 'michaelisConstant', c.michaelisConstant, ...
        'degradationRate', c.degradationRate, 'timeStep', c.timeStep, ...
        'developmentSteps', 600, 'payoffDC', 0, ...
        'sensorIntegrationRate', c.sensorIntegrationRate, 'backend', 'simulate_pair', ...
        'sigma_mut', c.mutationStdWeight);
    spE = sp; spE.developmentSteps = c.developmentSteps;   % the runs' step count (300)
end


% =========================================================================
function m = test_membership_reps(hc, sg, cfg)
% Per replicate: success (final cooperation >= 0.95); failing after 0.90 (not successful, and
% cooperation above 0.90 at some generation); never collapsing (successful, and after first
% reaching 0.90 never below 0.50). Per failing replicate: the fall begins at gOn, one after the
% last generation with cooperation >= 0.90 (none if final cooperation is still >= 0.90,
% status 'no_final_fall'); kPre is the last stored generation before gOn with cooperation
% >= 0.90 (status 'no_pre_snapshot' if none), kPost the first stored generation at or after
% gOn, kPost2 the next one. Events with status 'ok' enter the census.
    nR = size(hc, 1); fin = hc(:, end); pk = max(hc, [], 2);
    m.success = fin >= cfg.successThreshold;
    m.collapse = ~m.success & pk > cfg.arrive;
    m.neverCollapse = false(nR, 1);
    for r = 1:nR
        ga = find(hc(r, :) >= cfg.arrive, 1, 'first');
        m.neverCollapse(r) = m.success(r) && ~isempty(ga) && all(hc(r, ga:end) >= cfg.collapseFloor);
    end
    ev = struct('rep', {}, 'status', {}, 'gOn', {}, 'kPre', {}, 'kPost', {}, 'kPost2', {}, ...
        'pool', {}, 'ctrlRep', {});
    for r = find(m.collapse)'
        e = struct('rep', r, 'status', 'ok', 'gOn', NaN, 'kPre', NaN, 'kPost', NaN, 'kPost2', NaN, ...
            'pool', [], 'ctrlRep', NaN);
        if hc(r, end) >= cfg.arrive
            e.status = 'no_final_fall';
        else
            e.gOn = find(hc(r, :) >= cfg.arrive, 1, 'last') + 1;
            kp = find(sg < e.gOn & hc(r, sg) >= cfg.arrive, 1, 'last');
            if isempty(kp), e.status = 'no_pre_snapshot';
            else
                e.kPre = kp; e.kPost = find(sg >= e.gOn, 1, 'first');
                if e.kPost < numel(sg), e.kPost2 = e.kPost + 1; end
            end
        end
        ev(end+1) = e; %#ok<AGROW>
    end
    m.events = ev;
end


% =========================================================================
function m = pick_controls(m, cfg)
% Per 'ok' event, its pool: the never-collapsing replicates cooperating >= 0.90 at both its kPre
% and kPost generations. One pool replicate, drawn at random, is the event's control for the
% defector pipeline; for the class shares, up to nCtrlClass replicates are drawn from the union
% of the pools of the events sharing a kPre generation and tested for membership at that generation.
    ok = strcmp({m.events.status}, 'ok');
    for i = find(ok)
        e = m.events(i);
        pool = find(m.neverCollapse & m.hc(:, m.sg(e.kPre)) >= cfg.arrive ...
                    & m.hc(:, m.sg(e.kPost)) >= cfg.arrive)';
        m.events(i).pool = pool;
        if isempty(pool), m.events(i).ctrlRep = NaN;
        else, m.events(i).ctrlRep = pool(randi(numel(pool)));
        end
    end
    krestart = unique([m.events(ok).kPre]);
    m.classGen = struct('k', {}, 'reps', {});
    for k = krestart
        pool = unique([m.events(ok & [m.events.kPre] == k).pool]);
        pool = pool(randperm(numel(pool)));
        m.classGen(end+1) = struct('k', k, 'reps', pool(1:min(cfg.nCtrlClass, numel(pool))));
    end
end


% =========================================================================
function jobs = make_jobs(a, m, R, sp, spE)
    jobs = {};
    for i = find(strcmp({m.events.status}, 'ok'))
        e = m.events(i);
        jobs{end+1} = defector_job(a, i, true, R{e.rep}, e, m.sg, sp, spE); %#ok<AGROW>
        if isfinite(e.ctrlRep)
            jobs{end+1} = defector_job(a, i, false, R{e.ctrlRep}, e, m.sg, sp, spE); %#ok<AGROW>
        end
    end
    for q = 1:numel(m.classGen)
        k = m.classGen(q).k;
        for r = m.classGen(q).reps
            jobs{end+1} = struct('type', 'class', 'arch', a, 'rep', r, 'k', k, ...
                'M', R{r}.historyM(:,:,:,k), 'B', R{r}.historyB(:,:,:,k), 'sp', sp); %#ok<AGROW>
        end
    end
end


function j = defector_job(a, i, isColl, rep, e, sg, sp, spE)
    j = struct('type', 'inv', 'arch', a, 'event', i, 'isCollapse', isColl, ...
        'kPre', e.kPre, 'kPost', e.kPost, 'kPost2', e.kPost2, 'sg', sg, 'sp', sp, 'spE', spE, ...
        'M', rep.historyM, 'B', rep.historyB, 'P', rep.historyParents, ...
        'hc', double(rep.historyCooperation(:))');
end


% =========================================================================
function r = run_job(job, cfg)
    if strcmp(job.type, 'class')
        s = test_membership_snapshot(job.M, job.B, job.sp);
        r = struct('share', s.share, 'c1AmongC1', s.c1AmongC1);
        return;
    end
    restart = test_membership_snapshot(job.M(:,:,:,job.kPre), job.B(:,:,:,job.kPre), job.sp);
    r.share = restart.share; r.c1AmongC1 = restart.c1AmongC1;
    r.post = later_defector_at(job, job.kPost, restart, cfg);
    r.usedPost2 = false;
    if isfinite(job.kPost2)
        r.post2 = later_defector_at(job, job.kPost2, restart, cfg);
    else
        r.post2 = later_defector_at(job, NaN, restart, cfg);
    end
    r.inv = r.post;
    if ~r.post.found && r.post2.found, r.inv = r.post2; r.usedPost2 = true; end
    res = restart.cls <= 2;                         % the cooperators: genotypes meeting (C1)
    r.resG = restart.G(res, :); r.resCnt = restart.cnt(res);
end


% =========================================================================
function s = test_membership_snapshot(M, B, sp)
% Distinct (M, B) genotypes (one row [M(:)' B(:)'] per individual, unique rows), the class of
% each one's wiring by the membership test, and the class shares over all N individuals.
    N = size(M, 3);
    G = [reshape(double(M), 16, N); reshape(double(B), 16, N)]';
    [s.G, ~, s.ig] = unique(G, 'rows');
    s.cnt = accumarray(s.ig, 1);
    Wg = s.G(:, 1:16) .* s.G(:, 17:32);
    [Wu, ~, iw] = unique(Wg, 'rows');
    cw = zeros(size(Wu, 1), 1);
    for u = 1:size(Wu, 1)
        rr = membership_test(reshape(Wu(u, :), 4, 4), sp);
        cw(u) = 3 - rr.is_weak - rr.is_css;     % 1 member, 2 cooperates, does not retaliate, 3 other
    end
    s.cls = cw(iw);
    s.share = [sum(s.cnt(s.cls == 1)), sum(s.cnt(s.cls == 2)), sum(s.cnt(s.cls == 3))] / N;
    s.c1AmongC1 = s.share(2) / max(s.share(1) + s.share(2), eps);
    if s.share(1) + s.share(2) == 0, s.c1AmongC1 = NaN; end
end


% =========================================================================
function v = later_defector_at(job, k, restart, cfg)
% The defector that later took over, at stored generation k: under the runs' encounter rule,
% the most common wiring (cW, from the genotype counts cg) whose cooperation against a partner
% drawn uniformly from the other N - 1 is below 0.5; its genotype is the most common one
% carrying that wiring. Then its distances, lineage and encounters with the kPre cooperators.
    v = empty_defector();
    if ~isfinite(k), return; end
    g = job.sg(k); M = job.M(:,:,:,k); B = job.B(:,:,:,k); N = size(M, 3);
    G = [reshape(double(M), 16, N); reshape(double(B), 16, N)]';
    [Gu, ~, ig] = unique(G, 'rows'); cg = accumarray(ig, 1);
    [Wu, ~, iw] = unique(Gu(:, 1:16) .* Gu(:, 17:32), 'rows');
    cW = accumarray(iw, cg);
    A = encounter_rule_table(Wu, job.spE, cfg.engineDraws);
    Pm = (repmat(cW', numel(cW), 1) - eye(numel(cW))) / (N - 1);
    PC = sum(Pm .* A, 2);
    v.gen = g; v.coopRec = job.hc(g); v.coopPred = sum(cW .* PC) / N;
    def = PC < 0.5;
    v.defShare = sum(cW(def)) / N;
    if ~any(def), return; end
    [~, u] = max(cW .* def);
    v.found = true; v.share = cW(u) / N; v.shareOfDef = cW(u) / sum(cW(def));
    v.PCinv = PC(u);
    gi = find(iw == u); [~, b] = max(cg(gi)); gI = gi(b);
    v.G = Gu(gI, :); v.W = Wu(u, :);
    v = distances(v, restart, job.sp);
    carriers = find(ig == gI);
    v = lineage(v, job, k, carriers, restart);
    v = exploitation(v, restart, job, cfg);
    rr = membership_test(reshape(v.W, 4, 4), job.sp);
    v.invClass = 3 - rr.is_weak - rr.is_css;
end


function v = empty_defector()
    v = struct('found', false, 'gen', NaN, 'coopRec', NaN, 'coopPred', NaN, 'defShare', NaN, ...
        'share', NaN, 'shareOfDef', NaN, 'PCinv', NaN, 'G', nan(1, 32), 'W', nan(1, 16), ...
        'eMin', NaN, 'eMinAll', NaN, 'typeMin', '', 'wSlotsMin', NaN, 'nearCls', NaN, ...
        'presentPre', NaN, 'onGrid', NaN, 'onGridCls', NaN, 'onGridMemberShare', NaN, ...
        'eAnc', NaN, 'ancCls', NaN, 'ancIsDominant', NaN, 'originGap', NaN, ...
        'originSnap', NaN, 'eParent', NaN, 'parentCls', NaN, 'parentFreq', NaN, ...
        'E', NaN, 'pID', NaN, 'pDD', NaN, 'pCC', NaN, 'Edet', NaN, 'pIDdet', NaN, ...
        'pDDdet', NaN, 'invClass', NaN);
end


% =========================================================================
function A = encounter_rule_table(Wu, spE, nDraws)
% A(u,v) = share of encounters in which u cooperates when it meets v, under the runs' encounter
% rule (each agent started from U(0, 0.01)); nDraws encounters per unordered pair, both agents'
% actions used.
    nu = size(Wu, 1); C = zeros(nu); n = zeros(nu); E = ones(4);
    for u = 1:nu
        Wa = reshape(Wu(u, :), 4, 4);
        for w = u:nu
            Wb = reshape(Wu(w, :), 4, 4);
            for k = 1:nDraws
                spE.initialState1 = rand(4, 1) * 0.01;
                spE.initialState2 = rand(4, 1) * 0.01;
                oc = simulate_pair(E, Wa, [], E, Wb, [], spE);
                C(u, w) = C(u, w) + (oc(1) == 'C'); n(u, w) = n(u, w) + 1;
                C(w, u) = C(w, u) + (oc(2) == 'C'); n(w, u) = n(w, u) + 1;
            end
        end
    end
    A = C ./ n;
end


% =========================================================================
function v = distances(v, restart, sp)
% Edit distance of the defector's genotype to the kPre cooperators, nnz(M ~= M') + nnz(B ~= B')
% (a toggle changes M only, a weight change B only; a lower bound on mutations), and whether it
% lies on a cooperator's test grid: its wiring differs in one entry by at most 3 sigma_mut.
    res = restart.cls <= 2;
    e = edit_dist(v.G, restart.G);
    v.eMinAll = min(e);
    v.presentPre = sum(restart.cnt(e == 0));
    if ~any(res), return; end
    er = e; er(~res) = Inf;
    [v.eMin, j] = min(er);
    dM = nnz(v.G(1:16) ~= restart.G(j, 1:16)); dB = nnz(v.G(17:32) ~= restart.G(j, 17:32));
    v.typeMin = sprintf('%dT%dW', dM, dB);
    Wr = restart.G(:, 1:16) .* restart.G(:, 17:32);
    dW = abs(Wr - v.W);
    v.wSlotsMin = min(sum(dW(res, :) > 0, 2));
    v.nearCls = restart.cls(j);
    onG = res & sum(dW > 0, 2) == 1 & max(dW, [], 2) <= 3 * sp.sigma_mut + 1e-9;
    v.onGrid = any(onG);
    if v.onGrid
        v.onGridMemberShare = sum(restart.cnt(onG & restart.cls == 1)) / sum(restart.cnt(onG));
        [~, jj] = max(restart.cnt .* onG); v.onGridCls = restart.cls(jj);
    end
end


function e = edit_dist(g, Gs)
    e = sum(Gs(:, 1:16) ~= g(1:16), 2) + sum(Gs(:, 17:32) ~= g(17:32), 2);
end


% =========================================================================
function v = lineage(v, job, k, carriers, restart)
% Ancestors of the defector's carriers at kPre through historyParents (slot i of generation
% g + 1 descends from slot P(i, g) of generation g), and the stored generation where its exact
% genotype first appears on that lineage. Genotypes are seen only at stored generations, so a
% defector that arose and spread between two of them appears only as a genotype absent at kPre.
    P = job.P; gPost = job.sg(k); gPre = job.sg(job.kPre);
    idx = carriers(:);
    for g = gPost - 1:-1:gPre, idx = double(P(idx, g)); end
    ga = restart.ig(idx);
    anc = mode(ga);
    v.eAnc = edit_dist(v.G, restart.G(anc, :));
    v.ancCls = restart.cls(anc);
    [~, dom] = max(restart.cnt .* (restart.cls <= 2));
    v.ancIsDominant = double(anc == dom);
    i = carriers(find(ga == anc, 1));           % representative carrier
    slot = nan(1, k); slot(k) = i;
    for g = gPost - 1:-1:1
        i = double(P(i, g));
        kk = find(job.sg == g, 1);
        if ~isempty(kk) && kk < k, slot(kk) = i; end
    end
    gen = @(kk) [reshape(double(job.M(:,:,slot(kk),kk)), 1, 16), reshape(double(job.B(:,:,slot(kk),kk)), 1, 16)];
    kk = k;
    while kk > 1 && isequal(gen(kk - 1), v.G), kk = kk - 1; end
    v.originSnap = kk; v.originGap = gPost - job.sg(kk);
    if kk == 1, return; end
    pg = gen(kk - 1);
    v.eParent = edit_dist(v.G, pg);
    rr = membership_test(reshape(pg(1:16) .* pg(17:32), 4, 4), job.sp);
    v.parentCls = 3 - rr.is_weak - rr.is_css;
    Mk = job.M(:,:,:,kk-1); Bk = job.B(:,:,:,kk-1); N = size(Mk, 3);
    Gk = [reshape(double(Mk), 16, N); reshape(double(Bk), 16, N)]';
    v.parentFreq = mean(all(Gk == pg, 2));
end


% =========================================================================
function v = exploitation(v, restart, job, cfg)
% The kPre cooperators (wirings meeting (C1), weighted by individuals) against the defector,
% cooperator first: exploitDraws encounters under the runs' encounter rule, and one as the
% membership test develops a pair (600 steps, start 0.005). E = P(cooperator C, defector D);
% pID = P(defector D); pDD = P(both defect), so defecting back = pDD / pID.
    [Wr, cr] = member_wirings(restart);
    E = ones(4); WI = reshape(v.W, 4, 4); spE = job.spE;
    x = zeros(numel(cr), 4); xd = zeros(numel(cr), 4);
    for u = 1:numel(cr)
        WR = reshape(Wr(u, :), 4, 4);
        for k = 1:cfg.exploitDraws
            spE.initialState1 = rand(4, 1) * 0.01; spE.initialState2 = rand(4, 1) * 0.01;
            oc = simulate_pair(E, WR, [], E, WI, [], spE);
            x(u, :) = x(u, :) + outcome_row(oc) / cfg.exploitDraws;
        end
        xd(u, :) = outcome_row(simulate_pair(E, WR, [], E, WI, [], job.sp));
    end
    w = cr / sum(cr);
    m = w' * x; md = w' * xd;              % columns: CD, DD, CC, DC (cooperator first)
    v.E = m(1); v.pID = m(1) + m(2); v.pDD = m(2); v.pCC = m(3);
    v.Edet = md(1); v.pIDdet = md(1) + md(2); v.pDDdet = md(2);
end


function [Wr, cr] = member_wirings(restart)
    res = restart.cls <= 2;
    W = restart.G(res, 1:16) .* restart.G(res, 17:32);
    [Wr, ~, iw] = unique(W, 'rows');
    cr = accumarray(iw, restart.cnt(res));
end


function x = outcome_row(oc)
    x = [strcmp(oc, 'CD'), strcmp(oc, 'DD'), strcmp(oc, 'CC'), strcmp(oc, 'DC')];
end


% =========================================================================
function out = assemble(res, jobs, meta)
    out = meta;
    for a = 1:numel(meta)
        J = find(cellfun(@(j) j.arch == a, jobs));
        isDefector = cellfun(@(j) strcmp(j.type, 'inv'), jobs(J));
        JI = J(isDefector); JC = J(~isDefector);
        isC = cellfun(@(j) j.isCollapse, jobs(JI));
        out(a).coll = [res{JI(isC)}];
        out(a).collEvent = cellfun(@(j) j.event, jobs(JI(isC)));
        out(a).ctrl = [res{JI(~isC)}];
        out(a).ctrlEvent = cellfun(@(j) j.event, jobs(JI(~isC)));
        out(a).cls = [res{JC}];
        out(a).clsK = cellfun(@(j) j.k, jobs(JC));
        out(a).clsRep = cellfun(@(j) j.rep, jobs(JC));
    end
end


% =========================================================================
function out = shuffle_baseline(out, cfg)
% A shuffled pairing: each event's kPre cooperators against up to nShuffle defectors of other
% events of the same architecture (distance and E).
    for a = 1:numel(out)
        C = out(a).coll; n = numel(C);
        f = arrayfun(@(c) c.inv.found, C);
        sh = nan(n, 2);
        parfor i = 1:n
            rng(cfg.rngSeed + 100000 * a + i);
            if ~f(i), continue; end
            others = find(f); others(others == i) = [];
            if isempty(others), continue; end
            others = others(randperm(numel(others), min(cfg.nShuffle, numel(others))));
            [eS, dS] = deal(nan(numel(others), 1));
            for q = 1:numel(others)
                [eS(q), dS(q)] = cross_pair(C(i), C(others(q)).inv, out(a), cfg);
            end
            sh(i, :) = [mean(eS), mean(dS <= 1)];
        end
        out(a).shuffleE = sh(:, 1); out(a).shuffleDist = sh(:, 2);
    end
end


function [E, dmin] = cross_pair(ci, inv, oa, cfg)
    dmin = min(edit_dist(inv.G, ci.resG));
    W = ci.resG(:, 1:16) .* ci.resG(:, 17:32);
    [Wr, ~, iw] = unique(W, 'rows'); cr = accumarray(iw, ci.resCnt);
    spE = oa.spE;
    Eo = ones(4); WI = reshape(inv.W, 4, 4); x = zeros(numel(cr), 1);
    nd = max(10, round(cfg.exploitDraws / 5));
    for u = 1:numel(cr)
        WR = reshape(Wr(u, :), 4, 4);
        for k = 1:nd
            spE.initialState1 = rand(4, 1) * 0.01; spE.initialState2 = rand(4, 1) * 0.01;
            x(u) = x(u) + strcmp(simulate_pair(Eo, WR, [], Eo, WI, [], spE), 'CD') / nd;
        end
    end
    E = (cr' * x) / sum(cr);
end


% =========================================================================
% REPORT
% =========================================================================
function print_report(out, cfg)
    rng(cfg.rngSeed + 7);                       % bootstrap stream, reproducible on reportOnly
    fprintf('CENSUS BEFORE FAILURE  (seed %d; engineDraws %d per pair per role; exploitDraws %d; nBoot %d)\n', ...
        cfg.rngSeed, cfg.engineDraws, cfg.exploitDraws, cfg.nBoot);
    fprintf('The restart point is the last stored generation at or above 0.90 before the final fall; POST is the first stored generation after it.\n');
    print_inputs(out);
    print_events(out);
    print_item2(out, cfg);
    print_item3(out, cfg);
    print_item4(out, cfg);
    print_item5(out);
    print_later_defector_extras(out, cfg);
end


function print_inputs(out)
    fprintf('\n[0] INPUTS (file date | nOpp | payoffMode | development steps | success at >= 0.95)\n');
    for a = 1:numel(out)
        fprintf('  %-20s %s | %d | %s | %d | %d/%d\n', out(a).name, out(a).fileDate, out(a).nOpp, ...
            out(a).payoffMode, out(a).engineSteps, sum(out(a).success), numel(out(a).success));
    end
end


function print_events(out)
    fprintf('\n[1] COLLAPSES (failed: below 0.95 at generation 500 after exceeding 0.90)\n');
    fprintf('  %-20s %5s %8s %8s %4s %6s | %13s %15s %15s | %11s %11s\n', 'architecture', 'coll', ...
        'noFall', 'noPre', 'ok', 'ctrlP', 'onset med[rng]', 'restart->POST med[rng]', 'onset->POST', ...
        'coop restart point', 'coop POST');
    for a = 1:numel(out)
        ev = out(a).events; st = {ev.status}; ok = strcmp(st, 'ok'); e = ev(ok);
        gOn = [e.gOn]; gPre = out(a).sg([e.kPre]); gPost = out(a).sg([e.kPost]);
        cRestart = arrayfun(@(x) out(a).hc(x.rep, out(a).sg(x.kPre)), e);
        cPost = arrayfun(@(x) out(a).hc(x.rep, out(a).sg(x.kPost)), e);
        fprintf('  %-20s %5d %8d %8d %4d %6d | %4d [%3d,%3d] %5.0f [%3d,%3d] %5.0f [%3d,%3d] | %11.3f %11.3f\n', ...
            out(a).name, numel(ev), sum(strcmp(st, 'no_final_fall')), sum(strcmp(st, 'no_pre_snapshot')), ...
            sum(ok), sum(out(a).neverCollapse), round(median(gOn)), min(gOn), max(gOn), ...
            median(gPost - gPre), min(gPost - gPre), max(gPost - gPre), median(gPost - gOn), ...
            min(gPost - gOn), max(gPost - gOn), median(cRestart), median(cPost));
    end
    fprintf('  (ctrlP = never-collapse successful replicates: after first reaching 0.90 never below 0.50)\n');
    fprintf('\n  ENCOUNTER-RULE CHECK at POST: predicted focal cooperation (the runs'' encounter rule, pair table) vs recorded\n');
    for a = 1:numel(out)
        [pc, rc] = pred_rec(defector_of(out(a).coll, 'post')); [pk, rk] = pred_rec(defector_of(out(a).ctrl, 'post'));
        fprintf('  %-20s collapse: median |pred-rec| %.4f (median rec %.3f, pred %.3f) | control: %.4f (rec %.3f, pred %.3f)\n', ...
            out(a).name, median(abs(pc - rc)), median(rc), median(pc), median(abs(pk - rk)), median(rk), median(pk));
    end
end


function [p, r] = pred_rec(P)
    p = [P.coopPred]; r = [P.coopRec];
end


% -------------------------------------------------------------------------
function print_item2(out, cfg)
    fprintf('\n[2] CLASSES AT THE RESTART POINT (membership test, clone-weighted), collapse vs matched never-collapse control\n');
    fprintf('    D = mean over events of (collapse share - control mean at the same generation), 95%% bootstrap CI\n');
    names = {'member', 'cooperate_only', 'other'};
    fprintf('  %-20s %4s | %-17s %8s %8s %24s\n', 'architecture', 'n', 'class', 'collapse', 'control', 'D [95% CI]');
    for a = [1:numel(out), 0]
        sel = arch_sel(out, a);
        for c = 1:3
            [D, lo, hi, xm, cm, n] = item2(out, sel, @(s) s.share(c), cfg);
            fprintf('  %-20s %4d | %-17s %8.3f %8.3f %+8.3f [%+.3f, %+.3f]\n', arch_name(out, a), n, ...
                names{c}, xm, cm, D, lo, hi);
        end
        [D, lo, hi, xm, cm, n] = item2(out, sel, @(s) s.c1AmongC1, cfg);
        fprintf('  %-20s %4d | %-17s %8.3f %8.3f %+8.3f [%+.3f, %+.3f]  (cooperate_only among wirings that cooperate with their clone)\n', ...
            arch_name(out, a), n, 'cooperate_only|C1', xm, cm, D, lo, hi);
    end
end


function [D, lo, hi, xm, cm, nE] = item2(out, sel, fx, cfg)
% D = mean over events of (the event's share at kPre - the mean of its matched controls at that
% generation); events without a matched control are skipped. 95% bootstrap interval resampling
% events and, within each matched generation, control replicates.
    X = []; pid = []; pools = {};
    for a = sel
        C = out(a).coll; ks = [out(a).events(out(a).collEvent).kPre];
        for i = 1:numel(C)
            y = arrayfun(fx, out(a).cls(out(a).clsK == ks(i)));
            y = y(isfinite(y)); x = fx(C(i));
            if isempty(y) || ~isfinite(x), continue; end
            X(end+1) = x; pools{end+1} = y; pid(end+1) = numel(pools); %#ok<AGROW>
        end
    end
    nE = numel(X);
    if nE == 0, [D, lo, hi, xm, cm] = deal(NaN); return; end
    ym = cellfun(@mean, pools);
    D = mean(X - ym(pid)); xm = mean(X); cm = mean(ym(pid));
    Db = zeros(cfg.nBoot, 1);
    for b = 1:cfg.nBoot
        ie = randi(nE, nE, 1);
        yb = cellfun(@(y) mean(y(randi(numel(y), numel(y), 1))), pools);
        Db(b) = mean(X(ie) - yb(pid(ie)));
    end
    q = quantile(Db, [0.025 0.975]); lo = q(1); hi = q(2);
end


% -------------------------------------------------------------------------
function print_item3(out, cfg)
    fprintf('\n[3] DEFECTOR THAT TOOK OVER = dominant wiring at POST with focal cooperation < 0.5 (the runs'' encounter rule); POST2 if POST has none\n');
    fprintf('    distance = mutation edit count to the nearest cooperator at the restart point (a genotype that cooperates with its clone)\n');
    fprintf('  %-20s %-9s %4s %5s | %5s %6s %6s %5s %5s | %22s\n', 'architecture', 'arm', 'n', 'found', ...
        'd=0', '1 tog', '1 wt', 'd=2', 'd>=3', 'P(d<=1) [95% CP]');
    for a = [1:numel(out), 0]
        sel = arch_sel(out, a);
        [C, K] = pooled(out, sel);
        for arm = 1:2
            if arm == 1, S = C; lab = 'collapse'; else, S = K; lab = 'control'; end
            I = defector_of(S, 'inv'); f = [I.found]; I = I(f);
            dist_row(arch_name(out, a), lab, numel(S), I);
        end
        sh = [];
        for aa = sel, sh = [sh; out(aa).shuffleDist]; end %#ok<AGROW>
        sh = sh(isfinite(sh));
        fprintf('  %-20s %-9s %4d %5s | shuffled pairs with d<=1: %.3f (mean over events)\n', arch_name(out, a), ...
            'shuffled', numel(sh), '', mean(sh));
    end
    print_lineage(out, cfg);
end


function dist_row(name, lab, n, I)
    e = [I.eMin]; t = {I.typeMin};
    tog = sum(e == 1 & strcmp(t, '1T0W')); wt = sum(e == 1 & strcmp(t, '0T1W'));
    r = frac_ci(e <= 1);
    fprintf('  %-20s %-9s %4d %5d | %5d %6d %6d %5d %5d | %.3f [%.3f, %.3f]\n', name, lab, n, numel(I), ...
        sum(e == 0), tog, wt, sum(e == 2), sum(e >= 3), r(1), r(2), r(3));
end


function print_lineage(out, ~)
    fprintf('\n  LINEAGE of the defectors that took over (historyParents back to the restart point and to the genotype''s origin)\n');
    fprintf('  %-20s %4s | %11s %9s %9s | %11s %26s | %14s %29s | %10s | %11s\n', 'architecture', 'n', 'atRestart>0', ...
        'eAnc=0', 'eAnc<=1', 'anc=domCoop', 'anc mem/cooperate_only/oth', 'origin gap', 'parent mem/cooperate_only/oth', ...
        'eParent<=1', 'parent freq');
    for a = [1:numel(out), 0]
        C = pooled(out, arch_sel(out, a)); I = defector_of(C, 'inv'); I = I([I.found]);
        if isempty(I), continue; end
        ac = [I.ancCls]; pc = [I.parentCls]; pcv = pc(isfinite(pc));
        fprintf('  %-20s %4d | %11.3f %9.3f %9.3f | %11.3f %7d/%14d/%3d | %4.0f [%3.0f,%3.0f] %10d/%14d/%3d | %10.3f | %11.3f\n', ...
            arch_name(out, a), numel(I), mean([I.presentPre] > 0), mean([I.eAnc] == 0), mean([I.eAnc] <= 1), ...
            mean([I.ancIsDominant]), sum(ac == 1), sum(ac == 2), sum(ac == 3), median([I.originGap]), ...
            min([I.originGap]), max([I.originGap]), sum(pcv == 1), sum(pcv == 2), sum(pcv == 3), ...
            mean([I.eParent] <= 1 & isfinite([I.eParent])), median([I.parentFreq], 'omitnan'));
    end
    fprintf('  (atRestart>0: defector genotype already present at the restart point; eAnc: distance to its own ancestor at the restart point;\n');
    fprintf('   origin gap: generations from the snapshot where the exact genotype first appears on its lineage to POST;\n');
    fprintf('   parent: the lineage genotype at the snapshot before that, classified with the membership test)\n');
end


% -------------------------------------------------------------------------
function print_item4(out, cfg)
    fprintf('\n[4] EXPLOITATION: cooperators at the restart point (cooperate with their clone, clone-weighted) vs the defector that took over, cooperator first\n');
    fprintf('    E = P(cooperator C, defector D); ret = P(cooperator D | defector D); event = unit, 95%% bootstrap CI\n');
    fprintf('  %-20s %-9s %4s | %22s %8s %22s | %22s %8s %8s\n', 'architecture', 'arm', 'n', ...
        'E runs'' rule', 'P(def=D)', 'ret runs'' rule', 'E memb. test(600,sym)', 'P(def=D)', 'ret');
    for a = [1:numel(out), 0]
        sel = arch_sel(out, a); [C, K] = pooled(out, sel);
        for arm = 1:2
            if arm == 1, S = C; lab = 'collapse'; else, S = K; lab = 'control'; end
            I = defector_of(S, 'inv'); I = I([I.found]);
            r = exploit_row(I, cfg);
            fprintf('  %-20s %-9s %4d | %.3f [%.3f, %.3f] %8.3f %.3f [%.3f, %.3f] | %.3f [%.3f, %.3f] %8.3f %8.3f\n', ...
                arch_name(out, a), lab, numel(I), r.E, r.Eci, r.pID, r.ret, r.retci, r.Ed, r.Edci, r.pIDd, r.retd);
        end
        sh = [];
        for aa = sel, sh = [sh; out(aa).shuffleE]; end %#ok<AGROW>
        sh = sh(isfinite(sh)); b = boot_mean(sh, cfg.nBoot);
        fprintf('  %-20s %-9s %4d | %.3f [%.3f, %.3f]  (cooperators vs defectors of OTHER collapses)\n', ...
            arch_name(out, a), 'shuffled', numel(sh), b(1), b(2), b(3));
        I = defector_of(C, 'inv'); I = I([I.found]);
        fprintf('  %-20s disagreement: runs rule E<0.5 but membership test CD in %d of %d; runs rule E>=0.5 but membership test not CD in %d\n', ...
            '', sum([I.E] < 0.5 & [I.Edet] >= 0.5), numel(I), sum([I.E] >= 0.5 & [I.Edet] < 0.5));
    end
end


function r = exploit_row(I, cfg)
    E = [I.E]; pid = [I.pID]; ret = [I.pDD] ./ pid; ret(pid == 0) = NaN;
    Ed = [I.Edet]; pidd = [I.pIDdet]; retd = [I.pDDdet] ./ pidd; retd(pidd == 0) = NaN;
    b = boot_mean(E, cfg.nBoot); r.E = b(1); r.Eci = b(2:3);
    b = boot_mean(ret(isfinite(ret)), cfg.nBoot); r.ret = b(1); r.retci = b(2:3);
    b = boot_mean(Ed, cfg.nBoot); r.Ed = b(1); r.Edci = b(2:3);
    r.pID = mean(pid); r.pIDd = mean(pidd); r.retd = mean(retd, 'omitnan');
end


% -------------------------------------------------------------------------
function print_item5(out)
    fprintf('\n[5] GRID AGREEMENT: the defector on the C2 grid of some cooperator at the restart point (one slot, |dW| <= 3 sigma = 0.6)\n');
    fprintf('  %-20s %4s | %20s | %34s | %37s\n', 'architecture', 'n', 'on grid [95% CP]', ...
        'grid cooperator mem/cooperate_only', 'nearest cooperator mem/cooperate_only');
    for a = [1:numel(out), 0]
        C = pooled(out, arch_sel(out, a)); I = defector_of(C, 'inv'); I = I([I.found]);
        if isempty(I), continue; end
        g = [I.onGrid] == 1; r = frac_ci(g); gc = [I.onGridCls]; nc = [I.nearCls];
        fprintf('  %-20s %4d | %.3f [%.3f, %.3f] | %19d/%14d | %22d/%14d\n', arch_name(out, a), ...
            numel(I), r(1), r(2), r(3), sum(gc(g) == 1), sum(gc(g) == 2), sum(nc == 1), sum(nc == 2));
    end
    fprintf('  (a MEMBER with the defector on its grid is a case the membership test certifies as uninvadable)\n');
end


% -------------------------------------------------------------------------
function print_later_defector_extras(out, cfg)
    fprintf('\n[6] THE DEFECTOR THAT TOOK OVER, and the POST2 sensitivity\n');
    fprintf('  %-20s %4s %5s | %8s %8s %8s | %26s | %9s | %s\n', 'architecture', 'n', 'POST2', ...
        'share', 'ofDef', 'defPOST', 'def mem/cooperate_only/oth', 'PC(def)', 'POST2 defector: n, P(d<=1), mean E');
    for a = [1:numel(out), 0]
        C = pooled(out, arch_sel(out, a)); I = defector_of(C, 'inv'); f = [I.found]; I = I(f);
        if isempty(I), continue; end
        ic = [I.invClass]; P2 = defector_of(C, 'post2'); P2 = P2([P2.found]);
        b = boot_mean([P2.E], cfg.nBoot);
        fprintf('  %-20s %4d %5d | %8.3f %8.3f %8.3f | %7d/%14d/%3d | %9.3f | %d, %.3f, %.3f\n', arch_name(out, a), ...
            numel(I), sum([C.usedPost2]), median([I.share]), median([I.shareOfDef]), median([I.defShare]), ...
            sum(ic == 1), sum(ic == 2), sum(ic == 3), median([I.PCinv]), numel(P2), mean([P2.eMin] <= 1), b(1));
    end
    fprintf('  (share = the defector''s share of the POST population; ofDef = its share of defecting individuals;\n');
    fprintf('   defPOST = share of POST individuals whose focal cooperation < 0.5; medians over events)\n');
end


% -------------------------------------------------------------------------
function sel = arch_sel(out, a)
    if a == 0, sel = 1:numel(out); else, sel = a; end
end


function s = arch_name(out, a)
    if a == 0, s = 'POOLED'; else, s = out(a).name; end
end


function [C, K] = pooled(out, sel)
    C = [out(sel).coll]; K = [out(sel).ctrl];
end


function I = defector_of(S, fld)
% Records of the defector that later took over, for an arm; a 1x0 array when the arm is empty.
    if isempty(S), I = repmat(empty_defector(), 1, 0); else, I = [S.(fld)]; end
end


function r = frac_ci(x)
    x = logical(x(:)); n = numel(x); k = sum(x);
    if n == 0, r = [NaN NaN NaN]; return; end
    lo = 0; hi = 1;
    if k > 0, lo = betainv(0.025, k, n - k + 1); end
    if k < n, hi = betainv(0.975, k + 1, n - k); end
    r = [k / n, lo, hi];
end


function b = boot_mean(x, nBoot)
    x = x(isfinite(x)); n = numel(x);
    if n == 0, b = [NaN NaN NaN]; return; end
    m = zeros(nBoot, 1);
    for k = 1:nBoot, m(k) = mean(x(randi(n, n, 1))); end
    q = quantile(m, [0.025 0.975]); b = [mean(x), q(1), q(2)];
end
