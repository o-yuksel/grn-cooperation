function out = run_safety_search_without_rules(phase, varargin)
% RUN_SAFETY_SEARCH_WITHOUT_RULES  How safe does selecting members for safety make them without the
% two rules (reliable cooperation with the clone, a floor on defecting mutants)? Runs, from each of
% 10 starts, a search that keeps only membership.
%
% Produces : SI S3.4, "Two rules close loopholes that made a search without them look eight times
%            safer" (the lines run_safety_search_without_rules('report') prints)
% Inputs   : publication/saturating_contrast_100reps.mat   (written by run_signaling_architectures)
% Outputs  : publication/safety_search/member_setup.mat   (the starts: wirings and replicates)
%            publication/safety_search/member_main_300steps.mat   (one record per search; 'report'
%              reads this file)
%            publication/safety_search/member_log_<phase>.txt   (the lines each phase prints)
% Run      : run_safety_search_without_rules('setup') picks the starts; ('main', 300, 150, 1000, 10)
%              runs the searches (300 steps, 150 decision mutations per step, 1,000 test mutations,
%              10 starts); ('report') prints the medians
% Seed     : per start k, RandStream('mt19937ar', 'Seed', 1000 + k) draws the proposals and the
%            decision mutations; the 1,000 test mutations come from seed 99991, the same for every start
%
% Starts: the commonest final wiring of each successful saturating-contrast replicate (cooperation
% at least 0.95 at the end), in replicate order, kept if a member with clone readout at least 1e-6,
% the first 10. A wiring's breaking share is the share of its own mutations (the mutation
% distribution) whose result is not a member, by the membership test computed in batch
% (shared/member_batch, shared/simpair_batch).

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    pub = grn_data_file('publication');
    odir = fullfile(pub, 'safety_search'); if ~exist(odir, 'dir'), mkdir(odir); end
    diary(fullfile(odir, ['member_log_' phase '.txt'])); c = onCleanup(@() diary('off'));
    fprintf('\n===== MEMBER SEARCH %s %s UTC =====\n', phase, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    switch phase
        case 'setup',  out = setup(pub, odir);
        case 'main',   out = main_runs(odir, varargin{:});
        case 'report', out = report(odir);
    end
end


function O = setup(pub, odir)
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'results'); P = cssp();
    S0M = {}; S0B = {}; reps = [];
    for r = 1:numel(d.results)
        x = d.results{r};
        if x.historyCooperation(end) < 0.95, continue; end
        N = size(x.finalNetworkM, 3);
        G = [reshape(double(x.finalNetworkM), 16, N)', reshape(double(x.finalNetworkB), 16, N)'];
        [~, ia, ic] = unique(G, 'rows'); [~, k] = max(accumarray(ic, 1)); i = ia(k);
        M = x.finalNetworkM(:, :, i); B = x.finalNetworkB(:, :, i); W = double(M) .* double(B);
        e = membership_test(W, P); [m1, r1] = member_batch(W, true);
        ok = e.is_css && r1.rClone >= 1e-6;
        fprintf('rep %3d: dominant share %.3f, membership_test weak %d member %d; batch member %d, clone readout %.4f -> %d\n', ...
            r, mean(ic == k), e.is_weak, e.is_css, m1, r1.rClone, ok);
        if ok, S0M{end + 1} = M; S0B{end + 1} = B; reps(end + 1) = r; end %#ok<AGROW>
        if numel(reps) == 10, break; end
    end
    O = struct('S0M', {S0M}, 'S0B', {S0B}, 'reps', reps);
    save(fullfile(odir, 'member_setup.mat'), '-struct', 'O');
    fprintf('Starts: reps %s\n', mat2str(reps));
end


function O = main_runs(odir, nSteps, nDec, nFresh, nStarts)
    s = load(fullfile(odir, 'member_setup.mat')); nStarts = min(nStarts, numel(s.reps));
    jobs = struct('k', num2cell(1:nStarts), 'mode', 'S');
    R = cell(numel(jobs), 1); t0 = tic;
    p = gcp('nocreate'); if isempty(p), parpool('local', 4); end
    parfor j = 1:numel(jobs), R{j} = one_search(s, jobs(j).k, nSteps, nDec, nFresh); end
    O = struct('R', {R}, 'jobs', jobs, 'nSteps', nSteps, 'nDec', nDec, 'nFresh', nFresh, 'sec', toc(t0), 'reps', s.reps);
    for j = 1:numel(R)
        r = R{j};
        fprintf('start %2d %s: notMember %.3f -> %.3f, endC1 %.3f -> %.3f (fresh), accepted %d / %d (trait-fail %d), %.0f s\n', ...
            r.k, r.mode, r.f0.endMem, r.f1.endMem, r.f0.endC1, r.f1.endC1, r.nAcc, nSteps, r.nTraitFail, r.sec);
    end
    fprintf('Total wall %.0f s\n', O.sec);
    save(fullfile(odir, sprintf('member_main_%dsteps.mat', nSteps)), '-struct', 'O');
end


function r = one_search(s, k, nSteps, nDec, nFresh)
% Each step proposes one mutation; a proposal that stays a member replaces the current wiring when
% its breaking share on the step's nDec decision mutations (applied to both) is no larger. trace
% columns: current share, proposal share, accepted, ||W||_F.
    t = tic; M = s.S0M{k}; B = s.S0B{k};
    rs = RandStream('mt19937ar', 'Seed', 1000 + k);
    fr = draw_muts(RandStream('mt19937ar', 'Seed', 99991), nFresh);
    r.k = k; r.mode = 'S'; r.f0 = fresh_eval(M, B, fr); r.W0 = double(M) .* double(B);
    r.trace = nan(nSteps, 4); r.nAcc = 0; r.nTraitFail = 0;
    for step = 1:nSteps
        pr = draw_muts(rs, 1); dec = draw_muts(rs, nDec);
        [Mc, Bc] = apply_muts(M, B, pr); Wc = double(Mc) .* double(Bc);
        ok = member_batch(Wc, true);
        uCur = unsafe(M, B, dec); uCand = NaN; acc = false;
        if ~ok, r.nTraitFail = r.nTraitFail + 1;
        else, uCand = unsafe(Mc, Bc, dec); acc = uCand <= uCur; end
        if acc, M = Mc; B = Bc; r.nAcc = r.nAcc + 1; end
        r.trace(step, :) = [uCur, uCand, acc, norm(double(M) .* double(B), 'fro')];
    end
    r.f1 = fresh_eval(M, B, fr); r.M1 = M; r.B1 = B; r.W1 = double(M) .* double(B); r.sec = toc(t);
end


function u = unsafe(M, B, dr)
    [Mm, Bm] = apply_muts(M, B, dr); mem = member_batch(double(Mm) .* double(Bm), true);
    u = mean(~mem);
end


function f = fresh_eval(M, B, fr)
% Shares on the 1,000 test mutations (always the full membership test: ends membership, ends
% cooperation, ends retaliation only), plus the wiring's own readouts and counts.
    [Mm, Bm] = apply_muts(M, B, fr); [mem, info] = member_batch(double(Mm) .* double(Bm), true);
    f.endMem = mean(~mem); f.endC1 = mean(~info.weak); f.endRet = mean(info.weak & ~mem);
    W = double(M) .* double(B); [m0, i0] = member_batch(W, true);
    f.member = m0; f.weak = i0.weak; f.rClone = i0.rClone; f.worstRet = i0.worstRet; f.nDefMut = i0.nDefMut;
    f.minRM = i0.minRM; f.normW = norm(W, 'fro'); f.nEdges = nnz(W);
end


function O = report(odir)
% The median breaking share at the start and the end of the searches in member_main_300steps.mat,
% per start, with the split of the test mutations and the accepted counts.
    f = 'member_main_300steps.mat'; X = load(fullfile(odir, f)); R = X.R;
    fprintf('File %s: %d steps, decision sample %d, fresh sample %d, wall %.0f s\n', f, X.nSteps, X.nDec, X.nFresh, X.sec);
    rr = R(cellfun(@(r) r.mode == 'S', R)); v = @(fn) cellfun(fn, rr);
    u0 = v(@(r) r.f0.endMem); u1 = v(@(r) r.f1.endMem); O.S.u0 = u0; O.S.u1 = u1; O.S.acc = v(@(r) r.nAcc);
    fprintf('\n--- Search S (unsafe = mutant not a member) ---\n');
    fprintf(['rep   unsafe0 unsafe1 | endMem0 endMem1 endC1_0 endC1_1 endRet0 endRet1 | member1 rClone0 rClone1 ' ...
        'worstRet0 worstRet1 nDefMut0 nDefMut1 | ||W|| 0->1 edges 0->1 | acc tfail\n']);
    for i = 1:numel(rr)
        a = rr{i}.f0; b = rr{i}.f1;
        fprintf('%3d   %7.3f %7.3f | %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f | %7d %7.3f %7.3f %9.3f %9.3f %8d %8d | %5.2f->%5.2f %2d->%2d | %3d %5d\n', ...
            X.reps(rr{i}.k), u0(i), u1(i), a.endMem, b.endMem, a.endC1, b.endC1, a.endRet, b.endRet, b.member, a.rClone, b.rClone, ...
            a.worstRet, b.worstRet, a.nDefMut, b.nDefMut, a.normW, b.normW, a.nEdges, b.nEdges, rr{i}.nAcc, rr{i}.nTraitFail);
    end
    fprintf('median unsafe %.3f [%.3f, %.3f] -> %.3f [%.3f, %.3f]; accepted median %g [%d, %d]\n', median(u0), min(u0), max(u0), ...
        median(u1), min(u1), max(u1), median(O.S.acc), min(O.S.acc), max(O.S.acc));
    e1 = v(@(r) r.f1.endC1); e0 = v(@(r) r.f0.endC1); t0 = v(@(r) r.f0.endRet); t1 = v(@(r) r.f1.endRet);
    fprintf('census split: end C1 median %.3f -> %.3f; end retaliation only median %.3f -> %.3f\n', median(e0), median(e1), median(t0), median(t1));
end
