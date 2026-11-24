function out = run_member_source(phase, varargin)
% RUN_MEMBER_SOURCE  Do members come from the starting wirings or from mutation during the run?
% Reruns saturating-contrast populations 1 to 40 for 500 generations with mutation off, and with no
% starting response (the four weights from S_C and S_D to C and D zero in every starting
% individual), and compares them with the stored standard runs: Table S16 and its paragraph in
% SI S2.3.
%
% Produces : Table S16 (reached cooperation, mostly members, success, with Wilson intervals and
%            Fisher tests) and the numbers of the SI S2.3 paragraph "Either the starting wirings or
%            mutation can supply members", all in summary_member_source.txt
% Inputs   : publication/saturating_contrast_100reps.mat (its config and populations 1-40; written
%            by run_signaling_architectures)
% Outputs  : publication/member_source/ under the data root: check.mat; full_B_r<from>-<to>.mat and
%            full_C_r<from>-<to>.mat (each population's cooperation per generation and share of
%            members per snapshot); stored_members.mat; trace_C.mat; summary_member_source.txt; each
%            phase's printed lines in log_<phase>.txt (log_full_B.txt and log_full_C.txt for 'full')
% Run      : run_member_source('check') checks the starts and runs the response test (the 349 of
%              1,000, 48% and 24% of the SI S2.3 paragraph); ('full', 'B', 1:40) mutation off (about
%              5 min on 4 workers); ('full', 'C', 1:40) no starting response (about 6 min);
%              ('stored') the share of members of the stored runs; ('trace') the ancestry in no
%              starting response (needs 'full', 'C'); ('summary') writes the summary. Saved chunks
%              of 'full' are skipped on a rerun.
% Seed     : population k runs under rng(12345 + k, 'twister') (in evolve_population), the stored
%            runs' seeds; the 100 individuals sampled per snapshot under RandStream('mt19937ar',
%            'Seed', 777 + k)
%
% Arms (code letters; the stored standard runs are the comparison): 'B' mutation off, mutationRate = 0 from
% generation 1; 'C' no starting
% response: the standard start (the draws evolve_population makes for population k) with
% B(C,S_C), B(C,S_D), B(D,S_C), B(D,S_D) set to 0 in every individual and the connections kept,
% given as starting wirings (initPopulationM/B), so the random stream that follows is the standard
% one and mutation changes these weights as any other. Members: the membership test
% (membership_test) on 100 sampled individuals per snapshot, each distinct wiring tested once.
% Reached cooperation means mean cooperation above 0.90 at some generation, not success (at least
% 0.95 at generation 500).

    pub = grn_data_file('publication');
    odir = fullfile(pub, 'member_source');
    if ~exist(odir, 'dir'), mkdir(odir); end
    tag = phase; if strcmp(phase, 'full'), tag = [phase '_' varargin{1}]; end
    diary(fullfile(odir, ['log_' tag '.txt'])); c = onCleanup(@() diary('off'));
    fprintf('\n===== MEMBER SOURCE %s %s UTC =====\n', tag, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'config');
    base = d.config;
    switch phase
        case 'check',   out = check(base, pub, odir);
        case 'full',    out = full_run(base, odir, varargin{:});
        case 'stored',  out = stored_members(pub, odir);
        case 'trace',   out = trace_all(base, odir);
        case 'summary', out = summary(odir);
        otherwise, error('run_member_source:phase', 'unknown phase %s', phase);
    end
end


% =========================================================================
function cfg = arm_cfg(base, arm, k, ngen, snaps)
% Population k of the stored run: numReplicates = 1 runs replicate 1 under rngSeedBase + 1.
    cfg = base; cfg.numReplicates = 1; cfg.numGenerations = ngen;
    cfg.rngSeedBase = 12345 + k - 1;
    if nargin < 5, snaps = [1 10 25 50 100 150 200 300 500]; end
    cfg.alwaysSaveGenerations = snaps(snaps <= ngen); cfg.historySampleInterval = 100000;
    switch arm
        case 'S'                                   % standard (checks only)
        case 'B', cfg.mutationRate = 0;
        case 'C'
            [M, B] = standard_start(base, k);
            B(base.geneCooperator, [base.geneSensor1 base.geneSensor2], :) = 0;
            B(base.geneDefector,   [base.geneSensor1 base.geneSensor2], :) = 0;
            cfg.initPopulationM = M; cfg.initPopulationB = B;
        otherwise, error('run_member_source:arm', 'unknown arm %s', arm);
    end
end


function [M, B] = standard_start(base, k)
% The starting wirings evolve_population draws for population k, reproduced.
    rng(12345 + k, 'twister'); N = base.populationSize; G = base.numGenes;
    conn = 1.0; if isfield(base, 'initConnectivity'), conn = base.initConnectivity; end
    M = false(G, G, N); B = zeros(G, G, N, 'single');
    for i = 1:N
        M(:,:,i) = rand(G) < conn;
        B(:,:,i) = randn(G) * base.initSigma;
    end
end


function p = css_params()
% The membership test's settings: contrast sensing, saturating output, 600 development steps.
    p = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, ...
        'developmentSteps', 600, 'backend', 'simulate_pair');
end


function r = one_run(cfg, k)
% One population: its cooperation per generation and, at each snapshot, the share of members
% among 100 sampled individuals.
    t = tic; X = evolve_population(cfg); x = X{1}; r.secRun = toc(t);
    r.k = k; r.hc = double(x.historyCooperation(:))'; r.hg = double(x.historyGenerations(:))';
    nS = numel(r.hg); r.member = nan(1, nS); r.nDistinct = nan(1, nS); r.secClass = 0;
    r.W1 = double(x.historyM(:,:,:,1)) .* double(x.historyB(:,:,:,1));   % generation-1 population
    rs = RandStream('mt19937ar', 'Seed', 777 + k); p = css_params(); t = tic;
    for q = 1:nS
        sel = randperm(rs, size(x.historyM, 3), 100);
        W = reshape(double(x.historyM(:,:,sel,q)) .* double(x.historyB(:,:,sel,q)), 16, [])';
        [U, ~, ic] = unique(W, 'rows'); isM = false(size(U, 1), 1);
        for u = 1:size(U, 1), res = membership_test(reshape(U(u,:), 4, 4), p); isM(u) = res.is_css; end
        r.member(q) = mean(isM(ic)); r.nDistinct(q) = size(U, 1);
    end
    r.secClass = toc(t);
end


% =========================================================================
function O = check(base, sdir, odir)
% (1) the reproduced start equals the stored one; (2) the standard arm reproduces the stored run
% and the two arms start as described; (3) the response test of the Table S16 caption.
    S = load(fullfile(sdir, 'saturating_contrast_100reps.mat'), 'results'); R1 = S.results{1};
    st = R1.historyCooperation(:)';

    % (1) the reproduced start equals the stored generation-1 population of populations 1-4
    for k = 1:4
        Rk = S.results{k}; [M, B] = standard_start(base, k);
        fprintf('[START] rep %d: reproduced M equal %d, B equal %d (stored generation-1 snapshot)\n', k, ...
            isequal(M, Rk.historyM(:,:,:,1)), isequal(B, Rk.historyB(:,:,:,1)));
    end

    % (2) the runner's standard arm and the two arms' starts, population 1, 30 generations
    open_pool(4); arms = {'S', 'B', 'C'}; R = cell(1, 3);
    cfgs = cellfun(@(a) arm_cfg(base, a, 1, 30), arms, 'UniformOutput', false);
    parfor j = 1:3
        X = evolve_population(cfgs{j}); x = X{1};
        R{j} = struct('hc', double(x.historyCooperation(:))', 'M1', x.historyM(:,:,:,1), 'B1', x.historyB(:,:,:,1));
    end
    fprintf('[REGRESS] standard arm, rep 1, generations 1-30 identical to the stored run: %d (max |diff| %.3g)\n', ...
        isequal(single(R{1}.hc), st(1:30)), max(abs(R{1}.hc - double(st(1:30)))));
    fprintf('[REGRESS] standard arm generation-1 snapshot equal to stored: M %d, B %d\n', ...
        isequal(R{1}.M1, R1.historyM(:,:,:,1)), isequal(R{1}.B1, R1.historyB(:,:,:,1)));
    fprintf('[ARM B] gen-1 population equal to stored: M %d, B %d; coop gens 1-30 %s\n', ...
        isequal(R{2}.M1, R1.historyM(:,:,:,1)), isequal(R{2}.B1, R1.historyB(:,:,:,1)), mat2str(R{2}.hc, 3));
    ent = false(4); ent([base.geneCooperator base.geneDefector], [base.geneSensor1 base.geneSensor2]) = true;
    B0 = R1.historyB(:,:,:,1); Bc = R{3}.B1; z = reshape(Bc(repmat(ent, 1, 1, size(Bc, 3))), 4, []);
    Bo = B0; Bo(repmat(ent, 1, 1, size(B0, 3))) = 0;
    fprintf('[ARM C] gen-1: four entries zero in all individuals %d; other entries equal to stored %d; M equal to stored %d; coop gen 1 %.3f\n', ...
        all(z(:) == 0), isequal(Bc, Bo), isequal(R{3}.M1, R1.historyM(:,:,:,1)), R{3}.hc(1));

    % (3) generation-1 response test: every individual of population 1, arm C against the standard start
    pops = {'standard', double(R1.historyM(:,:,:,1)) .* double(R1.historyB(:,:,:,1)); ...
            'arm C', double(R{3}.M1) .* double(R{3}.B1)};
    O.resp = struct();
    for q = 1:2
        Wp = pops{q, 2}; n = size(Wp, 3); rc = nan(n, 1); rd = nan(n, 1); sp = sp_params(base);
        parfor i = 1:n
            [rc(i), rd(i)] = respond(Wp(:,:,i), sp);
        end
        tol = 1e-6; cl = @(r) [mean(r > tol), mean(r < -tol), mean(abs(r) <= tol)];
        flips = mean(sign_t(rc, tol) ~= sign_t(rd, tol));
        fprintf(['[RESPONSE] %-8s n=%d | partner signal +1 (cooperator): readout >0 %.3f, <0 %.3f, |.|<=1e-6 %.3f | ' ...
            'signal -1 (defector): >0 %.3f, <0 %.3f, ~0 %.3f | action differs between the two %.3f (%d individuals); ' ...
            'max |readout change| %.3g\n'], pops{q, 1}, n, cl(rc), cl(rd), flips, round(flips * n), max(abs(rc - rd)));
        % The action is readout > 0 (strict, as in the runs); the 1e-6 classes above move with the
        % stopping time, which the sensor input changes, so the action count is the test.
        fprintf('[RESPONSE] %-8s action (readout > 0) differs between cooperating and defecting partner: %d of %d\n', ...
            pops{q, 1}, sum((rc > 0) ~= (rd > 0)), n);
        O.resp.(strrep(pops{q, 1}, ' ', '_')) = struct('rc', rc, 'rd', rd);
    end
    save(fullfile(odir, 'check.mat'), 'O', 'R');
end


function [rc, rd] = respond(W, sp)
% Readout x_C - x_D of W against an all-zero partner whose C and D stay equal (partner signal 0),
% with a constant sensor input of a cooperating (+1) or defecting (-1) partner in contrast form
% (S_C gets +s, S_D gets -s, as contrast sensing in the runs).
    sp.externalSignal1 = [0; 0; 1; -1]; [~, a] = simulate_pair(W ~= 0, W, [], zeros(4), zeros(4), [], sp); rc = a(1) - a(2);
    sp.externalSignal1 = [0; 0; -1; 1]; [~, b] = simulate_pair(W ~= 0, W, [], zeros(4), zeros(4), [], sp); rd = b(1) - b(2);
end


function s = sign_t(r, tol), s = (r > tol) - (r < -tol); end


function sp = sp_params(base)
    sp = struct('hillCoefficient', 1, 'discreteOutput', false, 'sensingMode', 'contrast', ...
        'michaelisConstant', 1.0, 'degradationRate', base.degradationRate, 'timeStep', base.timeStep, ...
        'developmentSteps', base.developmentSteps, 'payoffDC', base.payoffDC);
end


% =========================================================================
% FULL RUN (arms B and C, 40 populations, 500 generations)
% =========================================================================
function O = full_run(base, odir, arm, reps)
    open_pool(4); snaps = [1 10:10:500]; tAll = tic;
    for c0 = 1:4:numel(reps)
        kk = reps(c0:min(c0 + 3, numel(reps)));
        fn = fullfile(odir, sprintf('full_%s_r%03d-%03d.mat', arm, kk(1), kk(end)));
        if exist(fn, 'file'), fprintf('skip %s (exists)\n', fn); continue; end
        jobs = arrayfun(@(k) arm_cfg(base, arm, k, 500, snaps), kk, 'UniformOutput', false);
        R = cell(numel(kk), 1); t = tic;
        parfor j = 1:numel(kk), R{j} = one_run(jobs{j}, kk(j)); end
        R = cellfun(@(r) rmfield(r, 'W1'), R, 'UniformOutput', false);
        for j = 1:numel(R)
            r = R{j}; a = find(r.hc > 0.90, 1); if isempty(a), a = NaN; end
            fprintf('[%s rep %2d] first > 0.90 at %g | final %.3f | max member share %.2f | run %.0f s, membership tests %.0f s\n', ...
                arm, r.k, a, r.hc(end), max(r.member), r.secRun, r.secClass);
        end
        save(fn, 'R');
        fprintf('chunk %s: %.1f min (total %.1f min)\n', fn, toc(t) / 60, toc(tAll) / 60);
    end
    O = [];
end


function R = load_full(odir, arm)
% The saved chunks of one arm, in population order.
    f = dir(fullfile(odir, sprintf('full_%s_r*.mat', arm))); R = {};
    for i = 1:numel(f), S = load(fullfile(odir, f(i).name)); R = [R; S.R(:)]; end %#ok<AGROW>
    k = cellfun(@(r) r.k, R); [~, o] = sort(k); R = R(o);
end


function O = stored_members(sdir, odir)
% The share of members of the stored standard populations 1-40 at their stored snapshots (same
% sampling as one_run).
    S = load(fullfile(sdir, 'saturating_contrast_100reps.mat'), 'results'); n = 40;
    HM = cell(n, 1); HB = cell(n, 1); hc = zeros(n, 500);
    for k = 1:n
        HM{k} = S.results{k}.historyM; HB{k} = S.results{k}.historyB;
        hc(k, :) = double(S.results{k}.historyCooperation(:))';
    end
    hg = double(S.results{1}.historyGenerations(:))'; clear S;
    open_pool(4); member = nan(n, numel(hg)); p = css_params();
    parfor k = 1:n
        rs = RandStream('mt19937ar', 'Seed', 777 + k); row = nan(1, numel(hg)); hm = HM{k}; hb = HB{k};
        for q = 1:numel(hg)
            sel = randperm(rs, size(hm, 3), 100);
            W = reshape(double(hm(:,:,sel,q)) .* double(hb(:,:,sel,q)), 16, [])';
            [U, ~, ic] = unique(W, 'rows'); isM = false(size(U, 1), 1);
            for u = 1:size(U, 1), res = membership_test(reshape(U(u,:), 4, 4), p); isM(u) = res.is_css; end
            row(q) = mean(isM(ic));
        end
        member(k, :) = row;
    end
    save(fullfile(odir, 'stored_members.mat'), 'hc', 'member', 'hg');
    fprintf('[STORED] members at stored generations %s computed for replicates 1-%d\n', mat2str(hg), n);
    O = [];
end


% =========================================================================
% TRACE (Arm C): ancestry of the first member wiring to reach 10 copies
% =========================================================================
function O = trace_all(base, odir)
    R = load_full(odir, 'C'); has = cellfun(@(r) any(r.member > 0), R); kk = cellfun(@(r) r.k, R(has));
    fprintf('[TRACE] Arm C populations with a sampled member at some snapshot: %d of %d (%s)\n', ...
        numel(kk), numel(R), mat2str(kk'));
    H = cellfun(@(r) r.hc, R(has), 'UniformOutput', false);
    open_pool(4); T = cell(numel(kk), 1);
    parfor j = 1:numel(kk)
        T{j} = trace_one(base, kk(j), H{j});
    end
    for j = 1:numel(T), print_trace(T{j}); end
    save(fullfile(odir, 'trace_C.mat'), 'T');
    O = T;
end


function t = trace_one(base, k, hcMain)
% Reruns population k of arm C storing every generation (it must equal the full run), finds the
% first member wiring with at least 10 copies and follows one carrier's ancestors to generation 1.
    cfg = arm_cfg(base, 'C', k, 500, []); cfg.historySampleInterval = 1; cfg.alwaysSaveGenerations = [];
    X = evolve_population(cfg); x = X{1}; t.k = k;
    t.same = isequal(double(x.historyCooperation(:))', hcMain);
    HM = x.historyM; HB = x.historyB; P = x.historyParents; nG = size(HM, 4); N = size(HM, 3);
    p = css_params(); sp = sp_params(base); cache = containers.Map();
    t.found = false; t.gStar = NaN; t.lastNRbeforeM = NaN; t.nStepsPath = 0; t.nMutPath = 0;
    t.nRespBetween = 0; t.via = 'none';
    for g = 1:nG
        W = reshape(double(HM(:,:,:,g)) .* double(HB(:,:,:,g)), 16, N)';
        [U, ia, ic] = unique(W, 'rows'); cnt = accumarray(ic, 1); cand = find(cnt >= 10);
        best = 0; bestU = 0;
        for u = cand(:)'
            c = member_class(reshape(U(u,:), 4, 4), cache, p, sp);
            if c == 2 && cnt(u) > best, best = cnt(u); bestU = u; end
        end
        if bestU > 0
            t.found = true; t.gStar = g; t.copies = best; t.carrier = ia(bestU); t.Wstar = reshape(U(bestU,:), 4, 4);
            break;
        end
    end
    if ~t.found, return; end
    % lineage of the carrier back to generation 1: individual j at g+1 has parent P(j, g) at g
    g0 = t.gStar; anc = zeros(1, g0); anc(g0) = t.carrier;
    for g = g0 - 1:-1:1, anc(g) = P(anc(g + 1), g); end
    LM = false(4, 4, g0); LB = zeros(4, 4, g0, 'single'); cls = zeros(1, g0);
    for g = 1:g0
        LM(:,:,g) = HM(:,:,anc(g),g); LB(:,:,g) = HB(:,:,anc(g),g);
        cls(g) = member_class(double(LM(:,:,g)) .* double(LB(:,:,g)), cache, p, sp);
    end
    t.cls = cls;   % 0 non-responder, 1 responder not member, 2 member
    t.lastNR = find(cls == 0, 1, 'last'); t.firstR = find(cls >= 1, 1); t.firstM = find(cls == 2, 1);
    % steps: generations g (2..g0) where the ancestor's (M, B) differs from generation g-1
    steps = struct('g', {}, 'ent', {}, 'kind', {}, 'so', {}, 'from', {}, 'to', {});
    for g = 2:g0
        dM = LM(:,:,g) ~= LM(:,:,g - 1); dB = LB(:,:,g) ~= LB(:,:,g - 1);
        if any(dM(:)) || any(dB(:))
            [ii, jj] = find(dM | dB); ent = [ii jj]; kind = cell(numel(ii), 1); so = false(numel(ii), 1);
            for e = 1:numel(ii)
                if dM(ii(e), jj(e)), kind{e} = 'toggle'; else, kind{e} = 'weight'; end
                so(e) = ii(e) <= 2 && jj(e) >= 3;
            end
            steps(end + 1) = struct('g', g, 'ent', ent, 'kind', {kind}, 'so', so, 'from', cls(g - 1), 'to', cls(g)); %#ok<AGROW>
        end
    end
    t.steps = steps;
    if isempty(t.firstM), t.firstM = NaN; end
    if isempty(t.lastNR), t.lastNR = NaN; end
    if isempty(t.firstR), t.firstR = NaN; end
    % path from the last non-responder ancestor before the first member ancestor to that member
    if ~isnan(t.firstM)
        lnr = find(cls(1:t.firstM) == 0, 1, 'last'); t.lastNRbeforeM = lnr;
        if isempty(lnr), t.lastNRbeforeM = NaN; t.path = steps([]); t.via = 'no non-responder ancestor';
        else
            inPath = arrayfun(@(s) s.g > lnr && s.g <= t.firstM, steps); t.path = steps(inPath);
            t.nMutPath = sum(arrayfun(@(s) size(s.ent, 1), t.path)); t.nStepsPath = numel(t.path);
            t.nRespBetween = sum(cls(lnr + 1:t.firstM - 1) == 1);
            if t.nRespBetween == 0, t.via = 'direct'; else, t.via = 'via responder'; end
        end
    end
end


function c = member_class(W, cache, p, sp)
% 2 member (is_css), 1 responder not member (action differs between partners), 0 non-responder;
% each wiring is tested once (cache).
    key = sprintf('%08x', typecast(single(W(:))', 'uint32'));
    if isKey(cache, key), c = cache(key); return; end
    res = membership_test(W, p);
    if res.is_css, c = 2;
    else, [rc, rd] = respond(W, sp); c = double((rc > 0) ~= (rd > 0)); end
    cache(key) = c; %#ok<NASGU>
end


function s = steps_str(st)
    names = {'C', 'D', 'Sc', 'Sd'}; cl = {'NR', 'R', 'M'}; s = '';
    for q = 1:numel(st)
        e = st(q).ent; parts = cell(1, size(e, 1));
        for i = 1:size(e, 1)
            parts{i} = sprintf('%s<-%s %s%s', names{e(i,1)}, names{e(i,2)}, st(q).kind{i}, repmat('*', 1, st(q).so(i)));
        end
        s = [s sprintf(' g%d[%s] %s->%s;', st(q).g, strjoin(parts, ','), cl{st(q).from + 1}, cl{st(q).to + 1})]; %#ok<AGROW>
    end
end


function print_trace(t)
    if ~t.found
        fprintf('[TRACE rep %2d] rerun identical %d | no member wiring reached 10 copies\n', t.k, t.same); return;
    end
    fprintf(['[TRACE rep %2d] rerun identical %d | first member >= 10 copies at gen %d (%d copies) | ancestor gens: ' ...
        'last NR %g, first R %g, first M %g | last NR before first M %g | path: %d steps, %d entry changes, ' ...
        'responder-non-member ancestors between %d -> %s |%s\n'], t.k, t.same, t.gStar, t.copies, t.lastNR, ...
        t.firstR, t.firstM, t.lastNRbeforeM, t.nStepsPath, t.nMutPath, t.nRespBetween, t.via, steps_str(t.path));
end


% =========================================================================
% SUMMARY
% =========================================================================
function O = summary(odir)
    S = load(fullfile(odir, 'stored_members.mat')); RB = load_full(odir, 'B'); RC = load_full(odir, 'C');
    n = 40; assert(numel(RB) == n && numel(RC) == n, 'full runs incomplete');
    A = struct();
    A.std = outcomes(S.hc, S.member, S.hg);
    A.B = outcomes(cell2mat(cellfun(@(r) r.hc, RB, 'UniformOutput', false)), cell2mat(cellfun(@(r) r.member, RB, 'UniformOutput', false)), RB{1}.hg);
    A.C = outcomes(cell2mat(cellfun(@(r) r.hc, RC, 'UniformOutput', false)), cell2mat(cellfun(@(r) r.member, RC, 'UniformOutput', false)), RC{1}.hg);
    fn = fullfile(odir, 'summary_member_source.txt'); fid = fopen(fn, 'w'); c = onCleanup(@() fclose(fid));
    w = @(varargin) fprintf(fid, varargin{:});
    w('MEMBER SOURCE, Table S16 (saturating contrast; replicates 1-40; seeds 12345 + k; 500 generations)\n');
    w('Arms: standard = stored saturating_contrast_100reps.mat replicates 1-40 (mutation on, standard start);\n');
    w('      B = standard start, mutation off; C = no starting response (sensor-to-output weights zero at the start), mutation on.\n');
    w('Members: membership_test is_css on 100 sampled individuals; arms B and C at generation 1 and every 10;\n');
    w('         standard at its stored snapshots %s.\n\n', mat2str(S.hg));
    labs = {'std', 'standard'; 'B', 'arm B'; 'C', 'arm C'};
    w('(a) ever mean cooperation > 0.90 | (b) ever sampled member share >= 0.5 | (c) cooperation >= 0.95 at generation 500\n');
    for q = 1:3
        o = A.(labs{q, 1});
        w('%-9s (a) %2d/40 [%.3f, %.3f] first gen median %s range %s | (b) %2d/40 [%.3f, %.3f] first snapshot median %s range %s | (c) %2d/40 [%.3f, %.3f]\n', ...
            labs{q, 2}, o.a, wilson(o.a, n), mstr(o.ga), rstr(o.ga), o.b, wilson(o.b, n), mstr(o.gb), rstr(o.gb), o.c, wilson(o.c, n));
    end
    for q = 2:3
        o = A.(labs{q, 1});
        w('Fisher two-sided, %s vs standard: (a) %d vs %d p = %.3g | (c) %d vs %d p = %.3g\n', labs{q, 2}, ...
            o.a, A.std.a, fp(o.a, A.std.a, n), o.c, A.std.c, fp(o.c, A.std.c, n));
    end
    w('\n(d) member share at generations 100, 200, 300, 500: median [range]\n');
    for q = 1:3
        o = A.(labs{q, 1}); w('%-9s', labs{q, 2});
        for g = 1:4, v = o.d(:, g); w(' g%d %.2f [%.2f, %.2f] |', o.dg(g), median(v), min(v), max(v)); end
        w('\n');
    end
    w('\nMean cooperation over populations at generations 1, 10, 25, 50, 100, 150, 200, 300, 500\n');
    gg = [1 10 25 50 100 150 200 300 500];
    for q = 1:3, o = A.(labs{q, 1}); w('%-9s %s\n', labs{q, 2}, mat2str(mean(o.hc(:, gg), 1), 3)); end
    w('\nPer population: k | first gen > 0.90 | first member >= 0.5 snapshot | coop g500 | members g100 g200 g300 g500\n');
    for k = 1:n
        w('%2d', k);
        for q = 1:3
            o = A.(labs{q, 1});
            w(' || %-8s %4s %4s %.3f %s', labs{q, 2}, nstr(o.ga(k)), nstr(o.gb(k)), o.hc(k, end), mat2str(o.d(k, :), 2));
        end
        w('\n');
    end
    tf = fullfile(odir, 'trace_C.mat');
    if exist(tf, 'file')
        L = load(tf); T = L.T;
        w('\nANCESTRY (arm C): the first member wiring to reach 10 copies (one carrier, traced to generation 1)\n');
        w('Classes: NR non-responder (action never differs between cooperating and defecting partner, direct sensor test),\n');
        w('         R responder not member, M member (is_css). Entries as target<-source; * marks a sensor-to-output entry.\n');
        w('Steps are generations at which the ancestor''s (M, B) changed; one step can hold more than one entry change.\n');
        nd = 0; nv = 0; nf = 0;
        for j = 1:numel(T)
            t = T{j};
            if ~t.found, w('rep %2d: rerun identical %d; no member wiring reached 10 copies\n', t.k, t.same); continue; end
            nf = nf + 1; nd = nd + strcmp(t.via, 'direct') * (t.nMutPath == 1); nv = nv + strcmp(t.via, 'via responder');
            w(['rep %2d: rerun identical %d; first member >= 10 copies at gen %d (%d copies); ancestors: last NR gen %g, ' ...
               'first R gen %g, first M gen %g; last NR before first M gen %g; path %d steps, %d entry changes, %d R-ancestor generations between; %s\n   path:%s\n' ...
               '   all steps from generation 1 to the first M ancestor:%s\n'], ...
                t.k, t.same, t.gStar, t.copies, t.lastNR, t.firstR, t.firstM, t.lastNRbeforeM, t.nStepsPath, ...
                t.nMutPath, t.nRespBetween, t.via, steps_str(t.path), steps_str(t.steps(arrayfun(@(s) s.g <= t.firstM, t.steps))));
        end
        nd2 = sum(cellfun(@(t) t.found && strcmp(t.via, 'direct'), T));
        w(['Traced populations with a member reaching 10 copies: %d; first member ancestor arose in ONE entry change from a ' ...
           'non-responder: %d (directly from a non-responder in any number of entry changes in one step: %d); through one or ' ...
           'more responder-non-member ancestors: %d\n'], nf, nd, nd2, nv);
    end
    fprintf('[SUMMARY] written %s\n', fn); type(fn);
    O = A;
end


function o = outcomes(hc, member, hg)
    n = size(hc, 1); o.hc = hc;
    o.ga = arrayfun(@(i) firstnan(find(hc(i, :) > 0.90, 1)), 1:n)'; o.a = sum(~isnan(o.ga));
    o.gb = arrayfun(@(i) firstnan(hg(find(member(i, :) >= 0.5, 1))), 1:n)'; o.b = sum(~isnan(o.gb));
    o.c = sum(hc(:, 500) >= 0.95);
    o.dg = [100 200 300 500]; o.d = zeros(n, 4);
    for g = 1:4, o.d(:, g) = member(:, find(hg == o.dg(g), 1)); end
end


function v = firstnan(x), if isempty(x), v = NaN; else, v = x; end, end
function s = nstr(x), if isnan(x), s = '-'; else, s = sprintf('%d', x); end, end
function s = mstr(x), x = x(~isnan(x)); if isempty(x), s = '-'; else, s = sprintf('%g', median(x)); end, end
function s = rstr(x), x = x(~isnan(x)); if isempty(x), s = '-'; else, s = sprintf('%g-%g', min(x), max(x)); end, end


function ci = wilson(x, n)
    z = 1.96; p = x / n; d = 1 + z^2 / n;
    c = (p + z^2 / (2 * n)) / d; h = z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; ci = [c - h, c + h];
end


function p = fp(a, b, n)
    [~, p] = fishertest([a, n - a; b, n - b]);
end


function open_pool(nw)
% A pool of nw workers, one thread each, opened only when none is open.
    maxNumCompThreads(1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', nw); end
    wait(parfevalOnAll(pool, @maxNumCompThreads, 0, 1));
end
