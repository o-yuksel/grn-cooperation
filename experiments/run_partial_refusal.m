function out = run_partial_refusal(phase, varargin)
% RUN_PARTIAL_REFUSAL  How much does refusing mutations that end cooperation protect it once it is
% built? Restarts the 100 saturating-contrast populations from the generation their cooperation first
% exceeded 0.90 and runs them to generation 500 with a share q of those mutations refused: the curve
% of Fig. 5.
%
% Produces : collapses per share refused, the share of these mutations left per arm and the pooled
%            Fisher tests, printed for Table S12 and drawn in Fig. 5 (fig5_refusal_curve reads the
%            runs and reports); also SI S4.6, Abstract, Significance, main text "How cooperating
%            populations lose cooperation", Discussion and Methods "Failure route"
% Inputs   : publication/saturating_contrast_100reps.mat (run_signaling_architectures);
%            publication/mutation_off/branch.mat, the start populations (run_refusal_start_populations);
%            mutation_off/runs_*.mat and report.mat, mutation on and off (run_mutation_off);
%            refused_by_class/rruns_coopend_*.mat, all refused, and rruns_response_*.mat ('report'
%              only) (run_all_refused)
% Outputs  : publication/partial_refusal/rruns_<arm>_*.mat (first batch, populations 1-100) and
%              partial_refusal_report.txt/.mat ('report': collapses and Fisher tests on these runs);
%            publication/partial_refusal_batch2/rruns_<arm>_*.mat (second batch, populations 1-50,
%              all eight arms) and partial_refusal_batch2_report.txt/.mat ('report_batch2': first 100, new 50, pooled 150);
%            each phase's printed lines in checks_log.txt, partial_refusal_log_main_<arm>.txt, main_log.txt
% Run      : first batch: run_partial_refusal('checks') runs q = 0.5 and 0.82 on populations 1-4 and
%              checks the refusal rule; ('main', arm) for each arm 'refuse050', 'refuse082', 'refuse090',
%              'refuse095', 'refuse098' (q = 0.5 to 0.98); ('report')   (about 95 min on 4 workers).
%            second batch: 'main_batch2' runs all eight arms, 'report_batch2'   (about 80 min on 4 workers)
% Seed     : rngSeedBase = 20261023 + k for population k (first batch), 30261023 + k (second),
%            Mersenne twister (in evolve_population); the refusing coin has its own stream (run_cfg)
%
% A mutation ends cooperation with the clone when the parent cooperates with its clone in at least
% 50 of 100 encounters and the mutant in fewer (mutate_with_refusal 'REFUSE_PARTIAL'); each such draw is
% refused with probability q and redrawn, up to 20 times.

    %% Parameters
    pub = grn_data_file('publication');
    sdir = fullfile(pub, 'partial_refusal');          % first batch
    sdir2 = fullfile(pub, 'partial_refusal_batch2');   % second batch
    switch phase   % the file of the printed lines: appended for runs, rewritten for reports
        case 'checks',   logFile = fullfile(sdir, 'checks_log.txt'); fmode = 'a';
        case 'main',    arm_p(varargin{1}); logFile = fullfile(sdir, ['partial_refusal_log_main_' varargin{1} '.txt']); fmode = 'a';
        case 'report',  logFile = fullfile(sdir, 'partial_refusal_report.txt'); fmode = 'w';
        case 'main_batch2',   logFile = fullfile(sdir2, 'main_log.txt'); fmode = 'a';
        case 'report_batch2', logFile = fullfile(sdir2, 'partial_refusal_batch2_report.txt'); fmode = 'w';
        otherwise, error('run_partial_refusal:phase', 'unknown phase %s', phase);
    end
    batch2 = any(strcmp(phase, {'main_batch2', 'report_batch2'}));
    stamp = char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss');

    %% Inputs
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'config');
    assert(~isfield(d.config, 'geneInputBias') || isempty(d.config.geneInputBias), 'config carries geneInputBias');
    L = load(fullfile(pub, 'mutation_off', 'branch.mat')); br = L.br;

    %% Compute
    if ~exist(fileparts(logFile), 'dir'), mkdir(fileparts(logFile)); end
    switch phase
        case 'checks',   [out, rep] = checks(d.config, br, pub, sdir);
        case 'main',    [out, rep] = main_runs(d.config, br, sdir, varargin{:});
        case 'report',  [out, rep] = report(br, pub, sdir);
        case 'main_batch2',   [out, rep] = main_batch2(d.config, br, sdir2, varargin{:});
        case 'report_batch2', [out, rep] = report_batch2(br, pub, sdir2);
    end

    %% Report
    if batch2
        head = {sprintf('\n===== FIG 5 SECOND BATCH %s %s UTC =====\n', phase, stamp), ...
            sprintf('data folder: %s\noutput .mat and printed lines: %s\n', pub, sdir2), ...
            sprintf('operator: %s\nevolve_population: %s\n', which('mutate_with_refusal'), which('evolve_population'))};
    else
        head = {sprintf('\n===== PARTIAL REFUSAL %s %s UTC =====\n', strjoin([{phase}, varargin(cellfun(@ischar, varargin))], ' '), stamp), ...
            sprintf('operator: %s\n', which('mutate_with_refusal'))};
    end
    rep = [head, rep];
    fprintf('%s', rep{:});

    %% Save
    fid = fopen(logFile, fmode); fprintf(fid, '%s', rep{:}); fclose(fid);
end


function p = arm_p(arm)
% The share q refused in each arm of the first batch.
    switch arm
        case 'refuse050', p = 0.5;
        case 'refuse082', p = 0.82;
        case 'refuse090', p = 0.90;
        case 'refuse095', p = 0.95;
        case 'refuse098', p = 0.98;
        otherwise, error('run_partial_refusal:arm', 'arm must be refuse050, refuse082, refuse090, refuse095 or refuse098');
    end
end


function cfg = run_cfg(base, b, p, ngen, s0)
% One continued run of population b.k with mutate_with_refusal with 'REFUSE_PARTIAL', refusing with
% probability p; s0 is the continuation seed base (20261023 first batch, the default; 30261023
% second batch). Populations are saved at continued generations 1, 41, 151 and the last.
% The coin that refuses (q is p here) has its own stream, seed refusalSeed =
% 1e9 + rngSeedBase + k, the same in every arm.
    if nargin < 5, s0 = 20261023; end
    cfg = base; cfg.numReplicates = 1; cfg.numGenerations = ngen;
    cfg.rngSeedBase = s0 + b.k;   % the continuation seed, same in every arm
    cfg.historySampleInterval = 100000;
    cfg.initPopulationM = b.M; cfg.initPopulationB = b.B;
    cfg.mutationOperator = @mutate_with_refusal; cfg.refusalMode = 'REFUSE_PARTIAL';
    cfg.refusalShare = p; cfg.refusalSeed = 1e9 + cfg.rngSeedBase + b.k;
    cfg.alwaysSaveGenerations = unique(min([1 41 151 ngen], ngen));
end


function r = r3_one_run(cfg)
% One refusing run: the cooperation history, the share defecting against its own clone at each
% saved generation, and the refusal's log (the same function as in run_all_refused).
    t = tic; mutate_with_refusal('reset');
    X = evolve_population(cfg); x = X{1}; L = mutate_with_refusal('log');
    r.hc = double(x.historyCooperation(:))'; r.sec = toc(t);
    r.hg = double(x.historyGenerations(:))'; r.defShare = nan(1, numel(r.hg));
    for q = 1:numel(r.hg), r.defShare(q) = def_share(x.historyM(:, :, :, q), x.historyB(:, :, :, q), cfg); end
    r.sampleP = L.sampleP; r.sampleM = L.sampleM; r.sampleS = L.sampleS;
    r.log = rmfield(L, {'sampleP', 'sampleM', 'sampleS'});
end


function sh = def_share(M, B, cfg)
% Share of individuals, among the commonest genotypes covering 95% of the population, whose
% cooperation against their clone is below 0.5.
    N = size(M, 3); [Wu, ~, ic] = unique(reshape(double(M) .* double(B), 16, N)', 'rows'); cnt = accumarray(ic, 1);
    [cs, ord] = sort(cnt, 'descend'); k = find(cumsum(cs) >= 0.95 * N, 1); sh = 0;
    for u = ord(1:k)'
        W = reshape(Wu(u, :), 4, 4);
        if mutate_with_refusal('coop', W, W, cfg) < 0.5, sh = sh + cnt(u); end
    end
    sh = sh / N;
end


function open_pool(nw)
% A pool of nw workers, one thread each, opened only when none is open.
    maxNumCompThreads(1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', nw); end
    wait(parfevalOnAll(pool, @maxNumCompThreads, 0, 1));
end


function [R, lab] = load_dir(dirn, pat)
% R and lab stacked from every file in dirn that matches pat.
    f = dir(fullfile(dirn, pat)); R = {}; lab = {};
    for i = 1:numel(f), S = load(fullfile(dirn, f(i).name)); R = [R; S.R]; lab = [lab; S.lab]; end %#ok<AGROW>
end


function s = pf(ok), if ok, s = 'PASS'; else, s = 'FAIL'; end, end


function t = first_or_inf(v)
    t = find(v, 1); if isempty(t), t = inf; end
end


function [O, rep] = checks(base, br, pub, sdir)
% One pool, 16 runs to generation 500: q = 1 and q = 0 on populations 1-2, checked against the
% stored all-refused and mutation-on runs (check1, check2); q = 0.5 and 0.82 on populations 1-4,
% saved as the first chunk of those arms, with the refused share checked against q (share_check);
% and q = 0.5 with the coin forced never to refuse, on which the same share check must fail.
% The printed lines give the checks and the refused share, not the saved runs' cooperation outcomes.
    ng = @(j) 501 - br(j).t0;
    jobs = {run_cfg(base, br(1), 1, ng(1)), run_cfg(base, br(2), 1, ng(2)), run_cfg(base, br(1), 0, ng(1)), run_cfg(base, br(2), 0, ng(2))};
    for p = [0.5 0.82], for j = 1:4, jobs{end + 1} = run_cfg(base, br(j), p, ng(j)); end, end %#ok<AGROW>
    for j = 1:4, c = run_cfg(base, br(j), 0.5, ng(j)); c.refusalNeverRefuse = true; jobs{end + 1} = c; end %#ok<AGROW>
    open_pool(4); R = cell(1, numel(jobs)); t = tic;
    parfor j = 1:numel(jobs), R{j} = r3_one_run(jobs{j}); end
    rep = {sprintf('[CHECKS] %d runs on 4 workers: wall %.1f min; seconds per run %s (generations %s)\n', numel(jobs), toc(t) / 60, ...
        mat2str(round(cellfun(@(r) r.sec, R))), mat2str(cellfun(@(r) numel(r.hc), R)))};
    [O.g1, r1] = check1(R(1:2), br, pub); [O.g2, r2] = check2(R(3:4), br, pub); rep = [rep, r1, r2];
    arms = {'refuse050', 'refuse082'};
    for a = 1:2
        Ra = R(4 + 4 * (a - 1) + (1:4)); [ok, r3] = share_check(Ra, arm_p(arms{a}), ['G3 ' arms{a}]);
        O.g3(a) = ok; rep = [rep, r3]; %#ok<AGROW>
        R0 = Ra(:); lab = [num2cell([br(1:4).k]'), repmat(arms(a), 4, 1), num2cell([br(1:4).t0]')];
        S.R = R0; S.lab = lab; save(fullfile(sdir, sprintf('rruns_%s_001-004.mat', arms{a})), '-struct', 'S');
    end
    [g4, r4] = share_check(R(13:16), 0.5, 'G4 refuse050 coin forced never to refuse'); rep = [rep, r4];
    O.g4 = ~g4;
    rep{end + 1} = sprintf('[G4] the share check returns FAIL on the forced coin -> G4 %s\n', pf(O.g4));
    rep{end + 1} = sprintf('[CHECKS] G1 %s, G2 %s, G3 refuse050 %s, G3 refuse082 %s, G4 %s\n', pf(O.g1), pf(O.g2), pf(O.g3(1)), pf(O.g3(2)), pf(O.g4));
end


function [ok, rep] = check1(R, br, pub)
% q = 1 against the stored all-refused runs ('coopend') of populations 1-2: every field both hold.
    S = load(fullfile(pub, 'refused_by_class', 'rruns_coopend_001-004.mat')); ok = true; rep = {};
    lf = {'firstState', 'calls', 'mut', 'refFirst', 'checks', 'refChecks', 'turned', 'redraws', 'fail', 'accN', 'accSq', 'firstSq', 'endFirst', 'nSeen'};
    for j = 1:2
        q = find(cell2mat(S.lab(:, 1)) == br(j).k); assert(numel(q) == 1); s0 = S.R{q}; r = R{j};
        same = [isequal(r.hc, s0.hc), isequal(r.hg, s0.hg), isequal(r.defShare, s0.defShare), ...
            isequal(r.sampleP, s0.sampleP), isequal(r.sampleM, s0.sampleM), cellfun(@(f) isequal(r.log.(f), s0.log.(f)), lf)];
        own = r.log.thinCoins == r.log.refChecks && r.log.endFirstRef == r.log.endFirst;
        ok = ok && all(same) && own;
        rep{end + 1} = sprintf('[G1] rep %d, %d generations: identical hc %d, hg %d, defShare %d, sampleP %d, sampleM %d; log fields %s identical %s; every coin refused %d\n', ...
            br(j).k, numel(r.hc), same(1:5), strjoin(lf, ','), mat2str(same(6:end)), own); %#ok<AGROW>
    end
    rep{end + 1} = sprintf('[G1] p = 1 reproduces stored all-refused reps 1-2 -> %s\n', pf(ok));
end


function [ok, rep] = check2(R, br, pub)
% q = 0 against the stored mutation-on runs ('control') of populations 1-2 (mutation_off/runs_*.mat,
% which hold hc and sec only).
    [R0, lab] = load_dir(fullfile(pub, 'mutation_off'), 'runs_*.mat'); ok = true; rep = {};
    for j = 1:2
        q = find(cell2mat(lab(:, 1)) == br(j).k & strcmp(lab(:, 2), 'control')); assert(numel(q) == 1);
        r = R{j}; same = isequal(r.hc, R0{q}.hc); none = r.log.refChecks == 0 && r.log.thinCoins > 0;
        ok = ok && same && none;
        rep{end + 1} = sprintf('[G2] rep %d, %d generations: hc identical to stored control %d; refusals %d, coins drawn %d\n', ...
            br(j).k, numel(r.hc), same, r.log.refChecks, r.log.thinCoins); %#ok<AGROW>
    end
    rep{end + 1} = sprintf('[G2] p = 0 reproduces stored mutation-on reps 1-2 -> %s\n', pf(ok));
end


function [ok, rep] = share_check(R, p, name)
% Share of the first draws that end cooperation which are refused, pooled over the runs, within
% 0.05 of q.
    lg = cellfun(@(r) r.log, R); nE = sum([lg.endFirst]); nR = sum([lg.endFirstRef]);
    sh = nR / nE; ok = abs(sh - p) <= 0.05;
    rep = {sprintf('[%s] p = %.2f; generations %s; ending first draws %d (per rep %s), refused %d: share %.4f (binomial SD at p %.4f); all coins %d, refused %d (%.4f); redraw-limit failures %d -> %s\n', ...
        name, p, mat2str(cellfun(@(r) numel(r.hc), R)), nE, mat2str([lg.endFirst]), nR, sh, sqrt(p * (1 - p) / nE), ...
        sum([lg.thinCoins]), sum([lg.refChecks]), sum([lg.refChecks]) / sum([lg.thinCoins]), sum([lg.fail]), pf(ok))};
end


function [O, rep] = main_runs(base, br, sdir, arm, chunk)
% One arm in chunks (default 16): populations 5-100 when 'checks' saved populations 1-4 (q = 0.5 and
% 0.82), else 1-100; a chunk already saved is skipped.
    if nargin < 5, chunk = 16; end
    p = arm_p(arm); open_pool(4); tAll = tic; rep = {};
    first = 5; if ~exist(fullfile(sdir, sprintf('rruns_%s_001-004.mat', arm)), 'file'), first = 1; end
    for c0 = first:chunk:numel(br)
        idx = c0:min(c0 + chunk - 1, numel(br));
        fn = fullfile(sdir, sprintf('rruns_%s_%03d-%03d.mat', arm, br(idx(1)).k, br(idx(end)).k));
        if exist(fn, 'file'), rep{end + 1} = sprintf('skip %s (exists)\n', fn); continue; end %#ok<AGROW>
        jobs = arrayfun(@(j) run_cfg(base, br(j), p, 501 - br(j).t0), idx, 'UniformOutput', false);
        lab = [num2cell([br(idx).k]'), repmat({arm}, numel(idx), 1), num2cell([br(idx).t0]')];
        R = cell(numel(jobs), 1); t = tic;
        parfor j = 1:numel(jobs), R{j} = r3_one_run(jobs{j}); end
        save(fn, 'R', 'lab');
        rep{end + 1} = sprintf('chunk %s: %.1f min (total %.1f min)\n', fn, toc(t) / 60, toc(tAll) / 60); %#ok<AGROW>
    end
    O = [];
end


function [O, rep] = report(br, pub, sdir)
% Collapses (first hc below 0.5) and lasting losses (collapsed and hc(end) below 0.95) per arm, with
% two-sided Fisher tests against the stored mutation-on, mutation-off and all-refused arms; the code
% asserts the stored arms' counts. The stored field 'location' holds run_all_refused's arm 'response' (refusing
% only the mutations that end cooperation through the response to the partner).
    kk = [br.k]'; n = numel(kk); rep = {};
    S1 = load(fullfile(pub, 'mutation_off', 'report.mat'));
    assert(isequal(S1.T.control.k(:), kk) && isequal(S1.T.off.k(:), kk), 'population order');
    A.on = mk(S1.T.control.cg, S1.T.control.fin); A.off = mk(S1.T.off.cg, S1.T.off.fin);
    A.all_refused = arm_counts(fullfile(pub, 'refused_by_class'), 'coopend', kk);
    A.location = arm_counts(fullfile(pub, 'refused_by_class'), 'response', kk);
    expect = struct('on', 36, 'off', 7, 'all_refused', 7, 'location', 26);
    for f = fieldnames(expect)'
        assert(sum(A.(f{1}).a) == expect.(f{1}), 'stored %s arm: %d collapses, expected %d', f{1}, sum(A.(f{1}).a), expect.(f{1}));
    end
    new = {'refuse050', 'refuse082', 'refuse090', 'refuse095', 'refuse098'}; got = {};
    for a = new
        x = arm_counts(sdir, a{1}, kk);
        if isempty(x), rep{end + 1} = sprintf('[REPORT] %s incomplete, not reported\n', a{1}); continue; end %#ok<AGROW>
        A.(a{1}) = x; got{end + 1} = a{1}; %#ok<AGROW>
    end
    rep{end + 1} = sprintf('\n[COUNTS] %d populations per arm; collapse = first drop below 0.5; lasting = collapsed and below 0.95 at the end\n', n);
    for a = [fieldnames(expect)', got], rep{end + 1} = sprintf('  %-10s collapses %3d  lasting %3d\n', a{1}, sum(A.(a{1}).a), sum(A.(a{1}).b)); end %#ok<AGROW>
    for a = got, [O.(a{1}), r] = read_arm(A, a{1}, n); rep = [rep, r]; end %#ok<AGROW>
    O.A = rmfield_runs(A); save(fullfile(sdir, 'partial_refusal_report.mat'), 'O');
end


function x = mk(cg, fin)
    cg = cg(:); fin = fin(:); x = struct('cg', cg, 'fin', fin, 'a', isfinite(cg), 'b', isfinite(cg) & fin < 0.95);
end


function x = arm_counts(dirn, arm, kk)
% Counts of a complete arm (empty if fewer than all populations); runs kept in x.R for the refusal's
% numbers.
    [R, lab] = load_dir(dirn, ['rruns_' arm '_*.mat']); x = [];
    if numel(R) < numel(kk), return; end
    [k, o] = sort(cell2mat(lab(:, 1))); R = R(o);
    assert(isequal(k(:), kk) && all(strcmp(stored_value('arm', lab(:, 2)), arm)));   % labels as stored in the files
    H = cellfun(@(r) r.hc, R, 'UniformOutput', false);
    x = mk(cellfun(@(h) first_or_inf(h < 0.5), H), cellfun(@(h) h(end), H)); x.R = R;
end


function A = rmfield_runs(A)
    for f = fieldnames(A)', if isfield(A.(f{1}), 'R'), A.(f{1}) = rmfield(A.(f{1}), 'R'); end, end
end


function [o, rep] = read_arm(A, arm, n)
% Fisher tests; the refusal's numbers.
    x = A.(arm); c = sum(x.a); p = arm_p(arm);
    fp = @(a, b) fisher_p(a, b, n); line = sprintf('[TESTS %s] collapses %d:', arm, c);
    for b = {'on', 'off', 'all_refused'}, line = [line sprintf(' vs %s (%d) p = %.3g;', b{1}, sum(A.(b{1}).a), fp(c, sum(A.(b{1}).a)))]; end %#ok<AGROW>
    rep = {sprintf('%s lasting %d vs on %d p = %.3g\n', line, sum(x.b), sum(A.on.b), fp(sum(x.b), sum(A.on.b)))};
    rep{end + 1} = sprintf('[ALSO REPORTED %s] vs location arm (%d): Fisher p = %.3g\n', arm, sum(A.location.a), fp(c, sum(A.location.a)));
    lg = cellfun(@(r) r.log, x.R); nE = sum([lg.endFirst]);
    o.share = sum([lg.endFirstRef]) / nE;
    rep{end + 1} = sprintf('  ending first draws %d, refused %d (share %.4f); ending draws accepted %d (%.4f of ending first draws, the count mutation on would accept); all refusals %d; redraw-limit failures %d\n', ...
        nE, sum([lg.endFirstRef]), o.share, sum([lg.endAcc]), sum([lg.endAcc]) / nE, sum([lg.refChecks]), sum([lg.fail]));
    rep{end + 1} = sprintf('  RMS step: drawn %.4f, accepted %.4f\n', sqrt(sum([lg.firstSq]) / sum([lg.mut])), sqrt(sum([lg.accSq]) / sum([lg.accN])));
    o.c = c;
end


function p = fisher_p(a, b, n)
    [~, p] = fishertest([a, n - a; b, n - b]);
end


% Second batch: populations 1-50 again with new continuation seeds, in all eight arms of Fig. 5.
% Arms: 'on' is q = 0 and 'all_refused' q = 1 (all refused), both through the same refusal step;
% 'off' is mutation off.
function a = arms2(), a = {'on', 'refuse050', 'refuse082', 'refuse090', 'refuse095', 'refuse098', 'all_refused', 'off'}; end
% second batch continuation seed base
function s = seed2base(), s = 30261023; end


function p = arm_p2(arm)
    switch arm
        case 'on', p = 0;
        case 'all_refused', p = 1;
        case 'off', p = NaN;
        otherwise, p = arm_p(arm);
    end
end


function cfg = cfg2(base, b, arm, ngen, s0)
% 'off' as run_mutation_off's arm off (mutation rate 0); every other arm REFUSE_PARTIAL at
% q = arm_p2(arm).
    if ~strcmp(arm, 'off'), cfg = run_cfg(base, b, arm_p2(arm), ngen, s0); return; end
    cfg = base; cfg.numReplicates = 1; cfg.numGenerations = ngen; cfg.mutationRate = 0;
    cfg.rngSeedBase = s0 + b.k;
    cfg.alwaysSaveGenerations = ngen; cfg.historySampleInterval = 100000;
    cfg.initPopulationM = b.M; cfg.initPopulationB = b.B;
end


function r = one_run2(cfg)
% Refusing arms as r3_one_run; the off arm as run_mutation_off's one_run (hc and sec).
    if isfield(cfg, 'refusalMode'), r = r3_one_run(cfg); return; end
    t = tic; X = evolve_population(cfg); x = X{1};
    r.hc = double(x.historyCooperation(:))'; r.sec = toc(t);
end


function [O, rep] = main_batch2(base, br, sdir, arms, chunk)
% Populations 1-50 per arm with the second-batch seeds, in chunks; a saved chunk is skipped.
    if nargin < 4 || isempty(arms), arms = arms2(); end
    if ischar(arms), arms = {arms}; end
    if nargin < 5, chunk = 16; end
    assert(all(ismember(arms, arms2())), 'unknown arm'); open_pool(4); tAll = tic; rep = {};
    for a = arms
        for c0 = 1:chunk:50
            idx = c0:min(c0 + chunk - 1, 50);
            fn = fullfile(sdir, sprintf('rruns_%s_%03d-%03d.mat', a{1}, br(idx(1)).k, br(idx(end)).k));
            if exist(fn, 'file'), rep{end + 1} = sprintf('skip %s (exists)\n', fn); continue; end %#ok<AGROW>
            jobs = arrayfun(@(j) cfg2(base, br(j), a{1}, 501 - br(j).t0, seed2base()), idx, 'UniformOutput', false);
            lab = [num2cell([br(idx).k]'), repmat(a, numel(idx), 1), num2cell([br(idx).t0]')];
            seeds = cellfun(@(c) c.rngSeedBase, jobs);
            R = cell(numel(jobs), 1); t = tic;
            parfor j = 1:numel(jobs), R{j} = one_run2(jobs{j}); end
            save(fn, 'R', 'lab', 'seeds');
            rep{end + 1} = sprintf('chunk %s: %d runs, seeds %d-%d, %.1f min (total %.1f min)\n', fn, numel(idx), seeds(1), seeds(end), toc(t) / 60, toc(tAll) / 60); %#ok<AGROW>
        end
    end
    O = [];
end


function [O, rep] = report_batch2(br, pub, sdir)
% Collapses and lasting losses per arm for the new 50, the first 100 and the pooled 150, with
% Clopper-Pearson 95% intervals (binofit); the share left (draws that end cooperation and are
% accepted over first draws that end cooperation) per arm, one minus the share removed of Table S12
% and Fig. 5; and Table S12's two-sided Fisher exact tests of each arm's pooled collapses against
% mutation on and against mutation off. The code asserts the first batch's counts.
    kk = [br.k]'; A = arms2(); expect = [36 33 27 17 21 10 7 7];
    S1 = load(fullfile(pub, 'mutation_off', 'report.mat'));
    assert(isequal(S1.T.control.k(:), kk) && isequal(S1.T.off.k(:), kk), 'population order');
    F.on = mk(S1.T.control.cg, S1.T.control.fin); F.off = mk(S1.T.off.cg, S1.T.off.fin);
    F.all_refused = arm_counts(fullfile(pub, 'refused_by_class'), 'coopend', kk);
    for a = A(2:6), F.(a{1}) = arm_counts(fullfile(pub, 'partial_refusal'), a{1}, kk); end
    for a = 1:8, N.(A{a}) = arm_counts(sdir, A{a}, kk(1:50)); assert(~isempty(N.(A{a})), 'new arm %s incomplete', A{a}); end
    for a = 1:8, assert(sum(F.(A{a}).a) == expect(a), 'first-batch %s: %d collapses, expected %d', A{a}, sum(F.(A{a}).a), expect(a)); end
    rep = {sprintf('\n[COUNTS] collapse = first drop below 0.5; lasting = collapsed and below 0.95 at the end; [Clopper-Pearson 95%%]\n')};
    rep{end + 1} = sprintf('Populations 1-50 enter the pooled count twice, 51-100 once.\n');
    rep{end + 1} = sprintf('%-10s %-5s | %-34s | %-14s | %-36s\n', 'arm', 'p', 'new 50: collapses / lasting', 'first 100', 'pooled 150: collapses / lasting');
    O = struct();
    for a = 1:8
        x = N.(A{a}); y = F.(A{a}); cn = [sum(x.a) sum(x.b)]; cf = [sum(y.a) sum(y.b)]; cp = cn + cf;
        rep{end + 1} = sprintf('%-10s %-5.2f | %2d %s / %2d %s | %3d / %3d      | %3d %s / %3d %s\n', A{a}, arm_p2(A{a}), ...
            cn(1), cpi(cn(1), 50), cn(2), cpi(cn(2), 50), cf(1), cf(2), cp(1), cpi(cp(1), 150), cp(2), cpi(cp(2), 150)); %#ok<AGROW>
        O.(A{a}) = struct('new', cn, 'first', cf, 'pooled', cp);
    end
    rep{end + 1} = sprintf('\n[SHARE LEFT] ending draws accepted / ending first draws (summed counts)\n');
    for a = 1:8
        if strcmp(A{a}, 'off'), continue; end
        ln = cellfun(@(r) r.log, N.(A{a}).R); e = [sum([ln.endAcc]) sum([ln.endFirst])];
        line = sprintf('  %-10s new 50 %d / %d = %.4f', A{a}, e(1), e(2), e(1) / e(2));
        if ismember(A{a}, A(2:6))
            lf = cellfun(@(r) r.log, F.(A{a}).R); f = [sum([lf.endAcc]) sum([lf.endFirst])]; s = e + f;
            line = [line sprintf(' | first 100 %d / %d = %.4f | pooled %d / %d = %.4f', f(1), f(2), f(1) / f(2), s(1), s(2), s(1) / s(2))]; %#ok<AGROW>
            O.(A{a}).share = [e(1) / e(2), s(1) / s(2)];
        else
            line = [line ' (check: p = 0 accepts all, p = 1 none; first batch has no REFUSE_PARTIAL log)']; %#ok<AGROW>
        end
        rep{end + 1} = sprintf('%s; redraw-limit failures %d\n', line, sum([ln.fail])); %#ok<AGROW>
    end
    % Table S12, columns p vs on and p vs off: each arm's pooled collapses (of 150) against those of
    % mutation on and of mutation off, two-sided Fisher exact tests, printed only
    on = O.on.pooled(1); off = O.off.pooled(1);
    rep{end + 1} = sprintf('\n[POOLED TESTS] two-sided Fisher exact p, pooled collapses of 150 against mutation on (%d) and mutation off (%d)\n', on, off);
    for a = 1:8
        c = O.(A{a}).pooled(1); line = sprintf('Table S12, %-10s pooled %3d of 150:', A{a}, c);
        if ~strcmp(A{a}, 'on'), line = [line sprintf(' p vs on %.2g;', fisher_p(c, on, 150))]; end %#ok<AGROW>
        if ~strcmp(A{a}, 'off'), line = [line sprintf(' p vs off %.2g;', fisher_p(c, off, 150))]; end %#ok<AGROW>
        rep{end + 1} = sprintf('%s\n', line(1:end - 1)); %#ok<AGROW>
    end
    % SI S4.6: the 90% and 95% points are out of order within chance; their pooled collapses
    % against each other, two-sided Fisher exact test, printed only
    c90 = O.refuse090.pooled(1); c95 = O.refuse095.pooled(1);
    rep{end + 1} = sprintf('SI S4.6, refuse090 against refuse095: pooled %d against %d of 150, two-sided Fisher exact p = %.2g\n', ...
        c90, c95, fisher_p(c90, c95, 150));
    save(fullfile(sdir, 'partial_refusal_batch2_report.mat'), 'O');
end


function s = cpi(x, n)
    [~, ci] = binofit(x, n); s = sprintf('[%.3f, %.3f]', ci(1), ci(2));
end
