function out = run_safety_search(phase, varargin)
% RUN_SAFETY_SEARCH  Can selection for safety leave a member that no mutation breaks? Selects 30
% members one mutation at a time while they must stay members and so retaliate (search 'S'), and
% from the same starts while they need only cooperate with their clone (search 'C').
%
% Produces : the searches behind Fig. 4B, Table S7, SI S3.4, main text "Conditional cooperation
%            stays within a mutation's reach of the threshold", Discussion and Methods "Mutation tests
%            in the game" (read by safety_search_losses('long'), fig4b_retaliation_from_rest and
%            safety_search_cooperate_only; summarized by run_safety_search('longreport'))
% Inputs   : publication/saturating_contrast_100reps.mat   (written by run_signaling_architectures)
%            publication/safety_search/long_flags.mat   ('longreport'; safety_search_losses('long'))
% Outputs  : in publication/safety_search/: search_setup.mat (the 30 starts); search_main_starts_1_4.mat and
%              search_main_starts_5_12.mat (300-step searches, starts 1-4 and 5-12; 'all' names the second
%              batch by the starts that fit its time limit, 5-12 in the stored run) and search_report.mat;
%              long_jobs/long_<k>_<S|C>.mat (one file per 1,000-step search), long_runs.mat,
%              long_report.mat; search_log_<phase>.txt (the lines each phase prints)
% Run      : run_safety_search('setup', 30); ('verify') checks the batch code; ('all', 30, 300, 150,
%              1000, 4200) runs 300-step searches on starts 1-4, re-checks every accepted wiring, then
%              adds starts within 4,200 s (56 min on 2 workers, stored log); ('longverify', 11);
%            ('long', 36000, 2000) runs the 1,000-step searches (3.8 h on 4 workers, stored log);
%              then safety_search_losses('long') and ('longreport').
% Seed     : per start k, RandStream('mt19937ar', 'Seed', 1000 + k) draws each step's proposal and
%            its 150 decision mutations, for both searches; the 1,000 test mutations, seed 99991
%
% Starts: the distinct commonest final wirings of successful saturating-contrast replicates, in
% replicate order, kept if a member with clone readout at least 1e-6 that cooperates with its clone
% in at least 90% of 100 encounters under the runs' encounter rule. Each step proposes one mutation
% from the mutation distribution. 'S' accepts a proposal that passes those three tests, has no fewer
% defecting single-weight mutants than the start, and has a breaking share no larger than the
% current wiring's on the decision mutations; 'C' drops membership and the floor on defecting
% mutants. The data folder is GRN_PUBDIR when set, else the data root's publication folder.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    pub = getenv('GRN_PUBDIR');
    if isempty(pub), pub = grn_data_file('publication'); end
    odir = fullfile(pub, 'safety_search');
    if ~exist(odir, 'dir'), mkdir(odir); end
    diary(fullfile(odir, ['search_log_' phase '.txt']));
    c = onCleanup(@() diary('off'));
    fprintf('\n===== SAFETY SEARCH %s %s UTC =====\n', phase, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    switch phase
        case 'setup',  out = setup(pub, odir, varargin{:});
        case 'verify', out = verify(odir);
        case 'main',   out = main_runs(odir, varargin{:});
        case 'check',  out = check_batch(odir, varargin{:});
        case 'report', out = report(odir, varargin{:});
        case 'all',    out = run_all(odir, varargin{:});
        case 'longverify', out = long_verify(odir, varargin{:});
        case 'long',       out = long_runs(odir, varargin{:});
        case 'longreport', out = long_report(odir, varargin{:});
    end
end


function Wd = default_defector()
% The default defector: gene D activates itself with weight 1.0, its only connection (the default defector of
%                       mutate_with_refusal, its default_defector).
    Wd = zeros(4); Wd(2, 2) = 1.0;
end


function O = setup(pub, odir, nWant)
% Distinct commonest final wirings of successful saturating-contrast replicates, in replicate order,
% kept if a member (membership_test), clone readout >= 1e-6 and cooperation with the clone >= 0.9 under
% the runs' encounter rule (mutate_with_refusal 'coop').
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'results', 'config'); P = cssp(); cfg = d.config;
    S0M = {}; S0B = {}; reps = []; keys = {};
    for r = 1:numel(d.results)
        x = d.results{r};
        if x.historyCooperation(end) < 0.95, continue; end
        N = size(x.finalNetworkM, 3);
        G = [reshape(double(x.finalNetworkM), 16, N)', reshape(double(x.finalNetworkB), 16, N)'];
        [~, ia, ic] = unique(G, 'rows'); [~, k] = max(accumarray(ic, 1)); i = ia(k);
        M = x.finalNetworkM(:, :, i); B = x.finalNetworkB(:, :, i); W = double(M) .* double(B);
        key = sprintf('%.17g,', W);
        if any(strcmp(keys, key)), fprintf('rep %3d: duplicate genotype, skipped\n', r); continue; end
        e = membership_test(W, P); [~, r1] = member_batch(W, true); cc = mutate_with_refusal('coop', W, W, cfg);
        ok = e.is_css && r1.rClone >= 1e-6 && cc >= 0.9;
        fprintf('rep %3d: share %.3f, member %d, clone readout %.4g, coop with clone (the runs'' encounter rule) %.2f, defecting grid mutants %d -> %d\n', ...
            r, mean(ic == k), e.is_css, r1.rClone, cc, r1.nDefMut, ok);
        if ok, S0M{end + 1} = M; S0B{end + 1} = B; reps(end + 1) = r; keys{end + 1} = key; end %#ok<AGROW>
        if numel(reps) == nWant, break; end
    end
    O = struct('S0M', {S0M}, 'S0B', {S0B}, 'reps', reps, 'cfg', cfg);
    save(fullfile(odir, 'search_setup.mat'), '-struct', 'O');
    fprintf('Starts (%d): reps %s\n', numel(reps), mat2str(reps));
end


function O = verify(odir)
% On the starts and 9 mutants of each: (i) coop_clone_batch equals mutate_with_refusal('coop', W, W)
% exactly; (ii) member_batch equals membership_test.
    s = load(fullfile(odir, 'search_setup.mat')); P = cssp(); rs = RandStream('mt19937ar', 'Seed', 4343);
    Ws = [];
    for k = 1:numel(s.reps)
        [Mm, Bm] = apply_muts(s.S0M{k}, s.S0B{k}, draw_muts(rs, 9));
        Ws = cat(3, Ws, double(s.S0M{k}) .* double(s.S0B{k}), double(Mm) .* double(Bm));
    end
    n = size(Ws, 3); cr = nan(n, 1); ew = false(n, 1); ec = false(n, 1); t = tic;
    for q = 1:n
        cr(q) = mutate_with_refusal('coop', Ws(:, :, q), Ws(:, :, q), s.cfg);
        e = membership_test(Ws(:, :, q), P); ew(q) = e.is_weak; ec(q) = e.is_css;
    end
    tRef = toc(t); t = tic; cb = coop_clone_batch(Ws); [bm, info] = member_batch(Ws, true); tBat = toc(t);
    O.n = n; O.sameCoop = sum(cr == cb); O.sameWeak = sum(ew == info.weak); O.sameMem = sum(ec == bm);
    O.pass = O.sameCoop == n && O.sameWeak == n && O.sameMem == n;
    fprintf(['Verify on %d wirings (starts + 9 model mutants each): coop under the runs'' encounter rule equal %d (max |diff| %.3g; ' ...
        '%d with coop >= 0.9), C1 equal %d, membership equal %d; reference %.1f s, batch %.1f s -> %s\n'], ...
        n, O.sameCoop, max(abs(cr - cb)), sum(cr >= 0.9), O.sameWeak, O.sameMem, tRef, tBat, ternary(O.pass, 'PASS', 'FAIL'));
end


function O = run_all(odir, nStarts, nSteps, nDec, nFresh, budgetSec)
% Searches on starts 1 to 4, then the check of every accepted wiring ('check'); if it passes, as many
% further starts (in setup order, up to nStarts) as fit the time limit `budgetSec` at the measured wall per
% start; stops if the check fails.
    s = load(fullfile(odir, 'search_setup.mat')); nStarts = min(nStarts, numel(s.reps)); t0 = tic;
    X = main_runs(odir, 1:4, nSteps, nDec, nFresh, 'first4');
    G = check_batch(odir, 'first4');
    if ~G.pass, fprintf('CHECK FAILED: stopping after 4 starts (%.0f s)\n', toc(t0)); O = G; return; end
    perStart = X.sec / 4; nMore = max(0, min(nStarts - 4, floor((budgetSec - toc(t0)) / (1.1 * perStart))));
    fprintf('Starts 1-4 wall %.0f s (%.0f s per start); budget %.0f s leaves room for %d more starts (of %d)\n', ...
        X.sec, perStart, budgetSec, nMore, nStarts - 4);
    tags = {'first4'};   % the second batch is named by the starts it ran (8 more: starts_5_12)
    if nMore > 0, tags{2} = sprintf('starts_%d_%d', 5, 4 + nMore); main_runs(odir, 5:(4 + nMore), nSteps, nDec, nFresh, tags{2}); end
    O = report(odir, tags);
    fprintf('run_all total wall %.0f s\n', toc(t0));
end


function O = main_runs(odir, ks, nSteps, nDec, nFresh, tag)
    s = load(fullfile(odir, 'search_setup.mat'));
    jobs = struct('k', {}, 'mode', {});
    for k = ks, jobs(end + 1) = struct('k', k, 'mode', 'S'); jobs(end + 1) = struct('k', k, 'mode', 'C'); end %#ok<AGROW>
    R = cell(numel(jobs), 1); t0 = tic;
    p = gcp('nocreate'); if isempty(p), parpool('local', 2); end
    parfor j = 1:numel(jobs), R{j} = one_search(s, jobs(j).k, jobs(j).mode, nSteps, nDec, nFresh); end
    O = struct('R', {R}, 'jobs', jobs, 'nSteps', nSteps, 'nDec', nDec, 'nFresh', nFresh, 'sec', toc(t0), 'reps', s.reps);
    for j = 1:numel(R)
        r = R{j}; u = r.fresh.(ternary(r.mode == 'S', 'uS', 'uC'));
        fprintf('start %2d (rep %3d) %s: unsafe %.3f -> %.3f (fresh), accepted %d / %d, nDef %d -> %d, %.0f s\n', ...
            r.k, s.reps(r.k), r.mode, u(1), u(end), r.nAcc, nSteps, r.nDef0, r.trace(end, 4), r.sec);
    end
    fprintf('Batch %s: wall %.0f s\n', tag, O.sec);
    save(main_file(odir, tag), '-struct', 'O');
end


function r = one_search(s, k, mode, nSteps, nDec, nFresh)
% A search of nSteps steps from start k, checkpoints every 50 steps: one_search_long from the start,
% without the stream state it keeps for a continuation.
    r = rmfield(one_search_long(s, k, mode, nSteps, nDec, nFresh, 50, []), 'rsState');
end


function u = unsafe(M, B, dr, mode)
    [Mm, Bm] = apply_muts(M, B, dr); ok = rule_ok(double(Mm) .* double(Bm), mode); u = mean(~ok);
end


function ok = rule_ok(Wm, mode)
% Search rule on mutants: (member if 'S', else C1) AND clone readout >= 1e-6 AND cooperation with
% the clone >= 0.9 under the runs' encounter rule. A logical `mode` is the full flag (true = 'S',
% false = 'C'), as one_search_long passes it.
% A mutant that fails this rule counts as broken; the floor on defecting mutants applies only to
% proposals.
    if islogical(mode), mode = ternary(mode, 'S', 'C'); end
    full = mode == 'S';
    [mem, info] = member_batch(Wm, full); base = ternary(full, mem, info.weak) & info.rClone >= 1e-6;
    cc = nan(size(base)); if any(base), cc(base) = coop_clone_batch(Wm(:, :, base)); end
    ok = base & cc >= 0.9;
end


function F = fresh_init(n)
    F = struct('uS', nan(n, 1), 'uC', nan(n, 1), 'endC1', nan(n, 1), 'unrel', nan(n, 1), 'endRet', nan(n, 1), ...
        'member', nan(n, 1), 'rClone', nan(n, 1), 'nDefMut', nan(n, 1), 'coopClone', nan(n, 1));
end


function F = fresh_put(F, q, f)
    fn = fieldnames(F); for i = 1:numel(fn), F.(fn{i})(q) = f.(fn{i}); end
end


function f = fresh_eval(M, B, fr, full)
% Both searches' breaking shares on the test sample, the S split, and the wiring's own state.
% full = false (C search, intermediate checkpoints only): C1 path only, S fields NaN.
    [Mm, Bm] = apply_muts(M, B, fr); Wm = double(Mm) .* double(Bm);
    [mem, info] = member_batch(Wm, full); base = info.weak & info.rClone >= 1e-6;
    cc = nan(size(base)); cc(base) = coop_clone_batch(Wm(:, :, base)); rel = base & cc >= 0.9;
    f.uC = mean(~rel); f.endC1 = mean(~info.weak); f.unrel = mean(info.weak & ~rel);
    f.uS = NaN; f.endRet = NaN; if full, f.uS = mean(~(mem & rel)); f.endRet = mean(rel & ~mem); end
    W = double(M) .* double(B); [m0, i0] = member_batch(W, true);
    f.member = m0; f.rClone = i0.rClone; f.nDefMut = i0.nDefMut; f.coopClone = coop_clone_batch(W);
end


function e = endpoint(M, B, cfg)
% Start and end: reference test results and responsiveness (alpha = ||(grad R)' Phi||_2 at the clone
% fixed point, as alpha_and_beta; alpha_p its sensor-channel part; ||x*||); cooperation with the
% default defector under the runs' encounter rule (coopDefRef).
    W = double(M) .* double(B); ce = membership_test(W, cssp());
    e.cssMember = ce.is_css; e.cssWeak = ce.is_weak;
    e.coopCloneRef = mutate_with_refusal('coop', W, W, cfg);
    e.coopDefRef = mutate_with_refusal('coop', W, default_defector(), cfg);
    [x, ~] = simpair_batch(W, W); x = x(:, 1);
    oC = x(1) / (x(1) + x(2) + 1e-12); oD = x(2) / (x(1) + x(2) + 1e-12); sg = [0; 0; oC - oD; -(oC - oD)];
    th = alpha_and_beta(W, x, sg, 0.1, 1.0, [1; -1; 0; 0]);
    e.alpha = th.alpha; e.alpha_p = norm(th.v(3:4)); e.xn = norm(x); e.valid = th.valid; e.x = x; e.W = W;
end


function O = check_batch(odir, tag)
% Every distinct accepted wiring re-checked with the reference implementations (membership_test,
% mutate_with_refusal coop) and the defecting-mutant floor; the trace's nDefMut never below start.
    X = load(main_file(odir, tag)); s = load(fullfile(odir, 'search_setup.mat')); P = cssp();
    nViol = 0; nChk = 0; nDefDrop = 0;
    for j = 1:numel(X.R)
        r = X.R{j}; full = r.mode == 'S';
        Wu = unique(r.Wtraj(logical(r.trace(:, 3)), :), 'rows');
        for q = 1:size(Wu, 1)
            W = reshape(Wu(q, :), 4, 4); e = membership_test(W, P); [~, inf1] = member_batch(W, true);
            cc = mutate_with_refusal('coop', W, W, s.cfg); nChk = nChk + 1;
            bad = ~ternary(full, e.is_css, e.is_weak) || inf1.rClone < 1e-6 || cc < 0.9 || (full && inf1.nDefMut < r.nDef0);
            if bad
                nViol = nViol + 1;
                fprintf('VIOLATION start %d %s: css %d weak %d rClone %.3g coop %.2f nDef %d (start %d)\n', ...
                    r.k, r.mode, e.is_css, e.is_weak, inf1.rClone, cc, inf1.nDefMut, r.nDef0);
            end
        end
        if full, nDefDrop = nDefDrop + any(r.trace(:, 4) < r.nDef0); end
    end
    O.nChk = nChk; O.nViol = nViol; O.nDefDrop = nDefDrop; O.pass = nViol == 0 && nDefDrop == 0;
    fprintf('Check %s: %d distinct accepted wirings re-checked by reference, %d violations; S searches with nDef below start: %d -> %s\n', ...
        tag, nChk, nViol, nDefDrop, ternary(O.pass, 'PASS', 'FAIL'));
end


function O = report(odir, tags)
    R = {}; sec = 0; files = cellfun(@(t) main_file(odir, t), tags, 'uni', 0);
    for i = 1:numel(tags)
        f = files{i}; if ~exist(f, 'file'), continue; end
        X = load(f); R = [R; X.R]; sec = sec + X.sec; %#ok<AGROW>
    end
    if isempty(R), error('run_safety_search:report', 'no search file found: %s', strjoin(files, ', ')); end
    s = load(fullfile(odir, 'search_setup.mat'));
    fprintf('%d searches, search wall %.0f s; steps %d, decision %d, fresh %d\n', numel(R), sec, X.nSteps, X.nDec, X.nFresh);
    for m = 'SC'
        rr = R(cellfun(@(r) r.mode == m, R)); fld = ternary(m == 'S', 'uS', 'uC');
        U = cell2mat(cellfun(@(r) r.fresh.(fld)', rr, 'UniformOutput', false)); O.(m).U = U; O.(m).ck = rr{1}.ckSteps;
        fprintf('\n--- Search %s (unsafe = %s) ---\n', m, ternary(m == 'S', 'mutant not a reliable member', 'mutant not a reliable clone cooperator'));
        fprintf(['k rep | unsafe0 unsafe1 | S split end: C1 unrel ret | nDef0 nDef1 | alpha0 alpha1 alphaP0 alphaP1 ' ...
            '||x||0 ||x||1 | coopClone0 1 | coopDef0 1 | css member1 | acc cmp | codes 1-4\n']);
        for i = 1:numel(rr)
            r = rr{i}; a = r.ends{1}; b = r.ends{2}; f = r.fresh; cd = histcounts(r.trace(:, 5), 0.5:1:4.5);
            fprintf('%2d %3d | %.3f %.3f | %.3f %.3f %.3f | %3d %3d | %6.3f %6.3f %6.3f %6.3f %5.2f %5.2f | %.2f %.2f | %.2f %.2f | %d | %3d %3d | %s\n', ...
                r.k, s.reps(r.k), U(i, 1), U(i, end), f.endC1(end), f.unrel(end), f.endRet(end), f.nDefMut(1), f.nDefMut(end), ...
                a.alpha, b.alpha, a.alpha_p, b.alpha_p, a.xn, b.xn, a.coopCloneRef, b.coopCloneRef, a.coopDefRef, b.coopDefRef, ...
                b.cssMember, r.nAcc, sum(r.trace(:, 5) == 0), mat2str(cd));
        end
        q = @(v) sprintf('%.3f [%.3f, %.3f] (IQR %.3f-%.3f)', median(v), min(v), max(v), prctile(v, 25), prctile(v, 75));
        fprintf('unsafe start %s\nunsafe end   %s\n', q(U(:, 1)), q(U(:, end)));
        fprintf('median by checkpoint: %s\n', mat2str(median(U, 1), 3));
        g = @(fn, w) cellfun(@(r) r.ends{w}.(fn), rr);
        fprintf('alpha start %s; end %s\nalpha_p start %s; end %s\n||x*|| start %s; end %s\n', q(g('alpha', 1)), q(g('alpha', 2)), ...
            q(g('alpha_p', 1)), q(g('alpha_p', 2)), q(g('xn', 1)), q(g('xn', 2)));
        fprintf('coop vs the default defector (the runs'' encounter rule): start %s; end %s; ends cooperating with it (>= 0.5): %d of %d\n', ...
            q(g('coopDefRef', 1)), q(g('coopDefRef', 2)), sum(g('coopDefRef', 2) >= 0.5), numel(rr));
        nd0 = cellfun(@(r) r.fresh.nDefMut(1), rr); nd1 = cellfun(@(r) r.fresh.nDefMut(end), rr);
        fprintf('defecting grid mutants start %s; end %s\n', q(nd0), q(nd1));
        O.(m).alpha = [g('alpha', 1), g('alpha', 2)]; O.(m).nd = [nd0, nd1];
    end
    S1 = O.S.U(:, end); C1 = O.C.U(:, end); nPos = sum(S1 > 0);
    fprintf('\nS ends above zero %d of %d; C median end %.4f\n', nPos, numel(S1), median(C1));
    save(fullfile(odir, 'search_report.mat'), '-struct', 'O');
end


function r = one_search_long(s, k, mode, nSteps, nDec, nFresh, ckEvery, R0)
% One search, 'S' or 'C', from start k: nSteps steps, the test sample scored every ckEvery steps.
% R0 empty starts at step 0; R0 a stored search continues it, replaying the proposal and decision
% stream (its draws do not depend on the search state), so it equals a fresh search of nSteps
% steps. trace columns: uCur, uCand (decision sample; NaN where the candidate failed the rules, so
% no comparison is needed), accepted, nDefMut of the current wiring after the step, code (0
% compared, 1 not member / not C1, 2 clone readout < 1e-6, 3 fewer defecting mutants, 4
% cooperation with the clone < 0.9 under the runs' encounter rule).
    t = tic; full = mode == 'S';
    rs = RandStream('mt19937ar', 'Seed', 1000 + k);
    fr = draw_muts(RandStream('mt19937ar', 'Seed', 99991), nFresh);
    if isempty(R0)
        M = s.S0M{k}; B = s.S0B{k}; s0 = 0;
        [~, i0] = member_batch(double(M) .* double(B), true); nDef0 = i0.nDefMut; nDefCur = ternary(full, nDef0, NaN);
        r.k = k; r.mode = mode; r.nDef0 = nDef0; r.trace = nan(nSteps, 5); r.nAcc = 0;
        r.Wtraj = nan(nSteps, 16); ck = 0:ckEvery:nSteps; r.ckSteps = ck;
        r.fresh = fresh_init(numel(ck)); r.fresh = fresh_put(r.fresh, 1, fresh_eval(M, B, fr, true));
        r.ends{1} = endpoint(M, B, s.cfg);
    else
        assert(R0.k == k && R0.mode == mode, 'continuation: wrong search');
        s0 = size(R0.trace, 1); nDef0 = R0.nDef0; nDefCur = R0.trace(end, 4); M = R0.M1; B = R0.B1;
        assert(isequal(reshape(double(M) .* double(B), 1, 16), R0.Wtraj(end, :)), 'continuation: end state');
        for step = 1:s0, draw_muts(rs, 1); draw_muts(rs, nDec); end
        r = R0; r.sec300 = R0.sec; r.trace = [R0.trace; nan(nSteps - s0, 5)]; r.Wtraj = [R0.Wtraj; nan(nSteps - s0, 16)];
        ckNew = (s0 + ckEvery):ckEvery:nSteps; ck = [R0.ckSteps, ckNew]; r.ckSteps = ck;
        fn = fieldnames(r.fresh); for i = 1:numel(fn), r.fresh.(fn{i}) = [r.fresh.(fn{i}); nan(numel(ckNew), 1)]; end
    end
    for step = s0 + 1:nSteps
        pr = draw_muts(rs, 1); dec = draw_muts(rs, nDec);
        [Mc, Bc] = apply_muts(M, B, pr); Wc = double(Mc) .* double(Bc);
        [okM, inf1] = member_batch(Wc, full); code = 0;
        if ~ternary(full, okM, inf1.weak), code = 1;
        elseif inf1.rClone < 1e-6, code = 2;
        elseif full && inf1.nDefMut < nDef0, code = 3;
        elseif coop_clone_batch(Wc) < 0.9, code = 4; end
        uCur = NaN; uCand = NaN; acc = false;
        if code == 0
            uCur = unsafe(M, B, dec, full); uCand = unsafe(Mc, Bc, dec, full); acc = uCand <= uCur;
        end
        if acc, M = Mc; B = Bc; r.nAcc = r.nAcc + 1; if full, nDefCur = inf1.nDefMut; end, end
        r.trace(step, :) = [uCur, uCand, acc, nDefCur, code];
        r.Wtraj(step, :) = reshape(double(M) .* double(B), 1, 16);
        q = find(ck == step, 1); if ~isempty(q), r.fresh = fresh_put(r.fresh, q, fresh_eval(M, B, fr, full || step == nSteps)); end
    end
    r.rsState = rs.State; r.ends{2} = endpoint(M, B, s.cfg); r.M1 = M; r.B1 = B; r.sec = toc(t);
end


function st = replay_state(k, nSteps, nDec)
    rs = RandStream('mt19937ar', 'Seed', 1000 + k);
    for step = 1:nSteps, draw_muts(rs, 1); draw_muts(rs, nDec); end
    st = rs.State;
end


function R = results_300steps(odir)
    R = [getfield(load(fullfile(odir, 'search_main_starts_1_4.mat'), 'R'), 'R'); ...
         getfield(load(fullfile(odir, 'search_main_starts_5_12.mat'), 'R'), 'R')]; %#ok<GFLD>
end


function r = pick(R, k, mode)
    r = R{cellfun(@(x) x.k == k && x.mode == mode, R)};
end


function O = long_verify(odir, k)
% Fresh 300-step S and C searches from start k (one_search_long) against the stored 300-step
% searches; the replayed stream state after 300 steps against the fresh run's.
    s = load(fullfile(odir, 'search_setup.mat')); R = results_300steps(odir);
    if isempty(gcp('nocreate')), parpool('local', 4); end
    modes = 'SC'; F = parallel.FevalFuture.empty;
    for m = modes, F(end + 1) = parfeval(@one_search_long, 1, s, k, m, 300, 150, 1000, 50, []); end %#ok<AGROW>
    O.pass = true;
    for i = 1:2
        r = fetchOutputs(F(i)); v = pick(R, k, modes(i));
        e = [isequaln(r.trace, v.trace), isequaln(r.Wtraj, v.Wtraj), isequaln(r.fresh, v.fresh), ...
             isequal(r.M1, v.M1) && isequal(r.B1, v.B1), r.nAcc == v.nAcc, isequal(r.rsState, replay_state(k, 300, 150))];
        fprintf(['verify start %d (rep %d) %s: trace %d, Wtraj %d, fresh %d, end wiring %d, accepted %d (%d vs %d), ' ...
            'replayed stream state %d; %.0f s (stored %.0f s)\n'], k, s.reps(k), modes(i), e(1:5), r.nAcc, v.nAcc, e(6), r.sec, v.sec);
        O.pass = O.pass && all(e); O.sec(i) = r.sec;
    end
    fprintf('long_verify -> %s\n', ternary(O.pass, 'PASS', 'FAIL'));
end


function r = long_job(s, j, odir)
    r = one_search_long(s, j.k, j.mode, j.nSteps, 150, 1000, j.ckEvery, j.R0);
    save(fullfile(odir, 'long_jobs', sprintf('long_%02d_%s.mat', j.k, j.mode)), 'r');
end


function O = long_runs(odir, budgetSec, estS)
% Queue in priority order: the 12 S searches continued from 300 to 1,000 steps, then for each of
% starts 13 to 30 S (1,000) and C (100). Queued new starts whose S could not finish within
% the time limit `budgetSec` (at estS s) are cancelled.
    s = load(fullfile(odir, 'search_setup.mat')); R = results_300steps(odir); t0 = tic;
    if ~exist(fullfile(odir, 'long_jobs'), 'dir'), mkdir(fullfile(odir, 'long_jobs')); end
    J = struct('k', {}, 'mode', {}, 'nSteps', {}, 'ckEvery', {}, 'R0', {});
    for k = 1:12, J(end + 1) = struct('k', k, 'mode', 'S', 'nSteps', 1000, 'ckEvery', 100, 'R0', pick(R, k, 'S')); end %#ok<AGROW>
    for k = 13:numel(s.reps)
        J(end + 1) = struct('k', k, 'mode', 'S', 'nSteps', 1000, 'ckEvery', 100, 'R0', []); %#ok<AGROW>
        J(end + 1) = struct('k', k, 'mode', 'C', 'nSteps', 100, 'ckEvery', 50, 'R0', []); %#ok<AGROW>
    end
    fprintf('Queue: %d jobs; new starts reps %s\n', numel(J), mat2str(s.reps(13:end)));
    if isempty(gcp('nocreate')), parpool('local', 4); end
    F = parallel.FevalFuture.empty;
    for j = 1:numel(J), F(j) = parfeval(@long_job, 1, s, J(j), odir); end
    done = false(1, numel(J)); cut = false(1, numel(J));
    while ~all(done | cut)
        pause(30);
        for j = find(~done & ~cut)
            st = F(j).State;
            if strcmp(st, 'finished')
                done(j) = true;
                if ~isempty(F(j).Error), fprintf('JOB %d ERROR: %s\n', j, F(j).Error.message); continue; end
                r = fetchOutputs(F(j)); u = r.fresh.(ternary(r.mode == 'S', 'uS', 'uC'));
                fprintf('[%5.0f s] start %2d (rep %3d) %s: unsafe %.3f -> %.3f, accepted %d / %d, %.0f s\n', ...
                    toc(t0), r.k, s.reps(r.k), r.mode, u(1), u(end), r.nAcc, size(r.trace, 1), r.sec);
            elseif strcmp(st, 'queued') && isempty(J(j).R0) && toc(t0) + estS > budgetSec
                ks = J(j).k; for i = find([J.k] == ks & ~done), cancel(F(i)); cut(i) = true; end
                fprintf('[%5.0f s] budget: start %d (rep %d) cut\n', toc(t0), ks, s.reps(ks));
            end
        end
    end
    O.J = rmfield(J, 'R0'); O.done = done; O.cut = cut; O.wall = toc(t0);
    fprintf('long_runs wall %.0f s; jobs done %d, cut %d\n', O.wall, sum(done), sum(cut));
    save(fullfile(odir, 'long_runs.mat'), '-struct', 'O');
end


function O = long_report(odir)
% The summary of the 1,000-step searches from the job files: the breaking share by checkpoint, its
% change from step 700 to 1,000, the C ends at zero, and the S ends above zero; the share that end
% cooperation comes from safety_search_losses('long') (long_flags.mat) when it exists.
    s = load(fullfile(odir, 'search_setup.mat')); V = results_300steps(odir);
    L = {}; d = dir(fullfile(odir, 'long_jobs', 'long_*.mat'));
    for i = 1:numel(d), L{end + 1} = getfield(load(fullfile(d(i).folder, d(i).name)), 'r'); end %#ok<AGROW,GFLD>
    RS = L(cellfun(@(r) r.mode == 'S' && size(r.trace, 1) == 1000, L));
    [~, o] = sort(cellfun(@(r) r.k, RS)); RS = RS(o); kS = cellfun(@(r) r.k, RS);
    ck = 0:100:1000; U = nan(numel(RS), numel(ck));
    for i = 1:numel(RS), for c = 1:numel(ck), U(i, c) = RS{i}.fresh.uS(RS{i}.ckSteps == ck(c)); end, end
    cont = kS <= 12;
    fprintf('S searches to 1,000 steps: %d (%d continued from the 300-step searches, %d new)\n', numel(RS), sum(cont), sum(~cont));
    fprintf('S median by checkpoint %s\n  all       %s\n  continued %s\n  new       %s\n', mat2str(ck), ...
        mat2str(median(U, 1), 3), mat2str(median(U(cont, :), 1), 3), mat2str(median(U(~cont, :), 1), 3));
    fprintf('S IQR at 1000: %.3f-%.3f; range %.3f-%.3f\n', prctile(U(:, end), 25), prctile(U(:, end), 75), min(U(:, end)), max(U(:, end)));
    fprintf(['k rep | uS 0 300 700 1000 | end split C1 unrel ret | nDef0 nDef1 | alpha0 alpha1 alphaP0 alphaP1 ||x||0 ||x||1 ' ...
        '| css member1 coopClone1 | acc/1000 | nDef below start | sec\n']);
    for i = 1:numel(RS)
        r = RS{i}; a = r.ends{1}; b = r.ends{2}; f = r.fresh;
        fprintf('%2d %3d | %.3f %.3f %.3f %.3f | %.3f %.3f %.3f | %3d %3d | %6.3f %6.3f %6.3f %6.3f %5.2f %5.2f | %d %.2f | %3d | %d | %.0f\n', ...
            r.k, s.reps(r.k), U(i, 1), U(i, 4), U(i, 8), U(i, end), f.endC1(end), f.unrel(end), f.endRet(end), ...
            f.nDefMut(1), f.nDefMut(end), a.alpha, b.alpha, a.alpha_p, b.alpha_p, a.xn, b.xn, b.cssMember, b.coopCloneRef, ...
            r.nAcc, any(r.trace(:, 4) < r.nDef0), r.sec);
    end
    q = @(v) sprintf('%.3f [%.3f, %.3f]', median(v, 'omitnan'), min(v), max(v));
    g = @(R, fn, w) cellfun(@(r) r.ends{w}.(fn), R);
    fprintf('S alpha start %s; end %s; ratio end/start median %.2f\nS alpha_p start %s; end %s\nS ||x*|| start %s; end %s\n', ...
        q(g(RS, 'alpha', 1)), q(g(RS, 'alpha', 2)), median(g(RS, 'alpha', 2) ./ g(RS, 'alpha', 1), 'omitnan'), ...
        q(g(RS, 'alpha_p', 1)), q(g(RS, 'alpha_p', 2)), q(g(RS, 'xn', 1)), q(g(RS, 'xn', 2)));
    % C ends: the 300-step searches for starts 1 to 12, the 100-step searches for starts 13 to 30.
    RC = {};
    for k = kS(:)'
        iC = find(cellfun(@(r) r.k == k && r.mode == 'C', L));
        if isempty(iC), RC{end + 1} = pick(V, k, 'C'); else, RC{end + 1} = L{iC}; end %#ok<AGROW>
    end
    uC1 = cellfun(@(r) r.fresh.uC(end), RC)'; nStepC = cellfun(@(r) size(r.trace, 1), RC)';
    newC = ~cont;
    fprintf('C ends: %d of %d at exactly 0 (median %.4f, max %.3f); new starts (100 steps): %d of %d at 0, values %s\n', ...
        sum(uC1 == 0), numel(uC1), median(uC1), max(uC1), sum(uC1(newC) == 0), sum(newC), mat2str(uC1(newC)', 3));
    if any(newC)
        uCnew = cell2mat(cellfun(@(r) r.fresh.uC', RC(newC)', 'UniformOutput', false));
        fprintf('C new median at steps 0, 50, 100: %s\n', mat2str(median(uCnew, 1), 3));
    end
    fprintf('C alpha start %s; end %s\n', q(g(RC, 'alpha', 1)), q(g(RC, 'alpha', 2)));
    % Summary lines (i) to (iii).
    dMed = median(U(:, ck == 700)) - median(U(:, end));
    O.levels = ternary(abs(dMed) < 0.01, 'LEVELS OFF', 'STILL FALLING');
    fprintf('\n(ii) S median 700 -> 1000: %.4f -> %.4f (change %.4f) -> %s\n', median(U(:, ck == 700)), median(U(:, end)), -dMed, O.levels);
    fprintf('(iii) C ends at 0: %d of %d\n', sum(uC1 == 0), numel(uC1));
    posS = U(:, end) > 0; fprintf('(i) S ends with uS > 0: %d of %d', sum(posS), numel(posS));
    ff = fullfile(odir, 'long_flags.mat');
    if exist(ff, 'file')
        F = load(ff, 'lostC', 'grp', 'k'); lc = arrayfun(@(k) F.lostC(F.grp == 2 & F.k == k), kS);
        O.immune = ternary(all(posS) && all(lc > 0), 'NEVER IMMUNE holds', 'NOT SHOWN');
        fprintf('; S ends with lostC > 0: %d of %d (min %.3f) -> %s\n', sum(lc > 0), numel(lc), min(lc), O.immune);
    else
        fprintf('; lostC pending (run safety_search_losses(''long''))\n');
    end
    O.U = U; O.ck = ck; O.kS = kS; O.uC1 = uC1; O.nStepC = nStepC;
    save(fullfile(odir, 'long_report.mat'), '-struct', 'O');
end

function f = main_file(odir, tag)
% The file of one batch of searches: batch 'first4' holds starts 1 to 4, batch 'rest' starts 5 to 12;
% any other tag names its file itself ('all' tags its second batch by its starts, e.g. 'starts_5_12').
    names = struct('first4', 'starts_1_4', 'rest', 'starts_5_12');
    if isfield(names, tag), tag = names.(tag); end
    f = fullfile(odir, ['search_main_' tag '.mat']);
end
