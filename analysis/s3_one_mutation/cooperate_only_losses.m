function cooperate_only_losses(stage)
% COOPERATE_ONLY_LOSSES  Do cooperators that do not retaliate lose cooperation to one mutation less
% often than members? Compares the share of mutations that end cooperation in the two groups,
% pooled and by architecture.
%
% Produces : SI S3.3, "Cooperators that do not retaliate lose cooperation less often than
%            members" (stage 'full', which prints the numbers)
% Inputs   : publication/rho_by_architecture.mat   (the census of members; per_mutation_census)
%            publication/cooperate_only_census/: full_cooperate_only.mat ('full'; written by
%              cooperate_only_census, run by run_cooperate_only_census) and control_member_10reps.mat
%              ('control'; the census code rerun on members of 10 replicates under a new seed, part
%              of the data archive)
% Outputs  : publication/cooperate_only_census/<stage>.txt   (the lines this run prints). 'full': the
%              share of mutations ending cooperation in each group with 95% intervals, their ratio,
%              by architecture; weighted medians of the cooperation margin r, beta, r / beta, alpha
%              and ||x*||; the shares by fifth of r / beta. 'control' checks the measuring code
%              against the census (no paper number)
% Run      : cooperate_only_losses('full')   (about 40 s)
% Seed     : bootstraps under rng(20261016 + 9, 'twister') pooled and 20261016 + 9 + m for
%            architecture m ('full' takes 20261016 from the stored run's rngSeed), around a saved
%            and restored stream
%
% Members: the census rows marked ok (cooperation margin at least 1e-6, ||x*|| in [0.5, 14]), on
% the replicates the other group's run covered; the others: the ok rows of full_cooperate_only.mat (the
% same rule, set by cooperate_only_census). Shares are weighted by the individuals carrying each
% wiring; intervals resample whole replicates within architecture, jointly for both groups.

    D = grn_data_file('publication');
    O = fullfile(D, 'cooperate_only_census');
    modeNames = {'graded_separate','graded_contrast','saturating_separate','saturating_contrast'};
    C = load_census(D);

    %% Report
    switch stage
        case 'control', report = @() do_control(C, O, modeNames); %#ok<NASGU> called inside evalc
        case 'full',    report = @() do_full(C, O, modeNames); %#ok<NASGU> called inside evalc
        otherwise, error('unknown stage');
    end
    rep = evalc('report();');
    fprintf('%s', rep);

    %% Save
    fid = fopen(fullfile(O, [stage '.txt']), 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


% =========================================================================
function T = load_census(D)
% The per-mutation census of members, all four architectures as one table, with its seed
    S = load(fullfile(D, 'rho_by_architecture.mat'), 'W_all', 'metadata');
    T = flatten(S.W_all{1}, 1);
    for m = 2:numel(S.W_all)
        F = flatten(S.W_all{m}, m);
        f = fieldnames(T);
        for i = 1:numel(f), T.(f{i}) = [T.(f{i}); F.(f{i})]; end
    end
    T.seed = S.metadata.cfg.rngSeed;
end

function T = flatten(R, m)
% One row per wiring. Each wiring stores its K draws' outcomes as a K-by-1 column; placing the
% columns side by side and transposing gives lost with wirings in rows and draws in columns.
% pLost is the wiring's share of draws ending cooperation (stored field 'pLostReadout').
    n = numel(R);
    if n == 0, T = []; return; end
    if nargin < 2 || isfield(R, 'mode') && ~isempty(R(1).mode) && R(1).mode > 0
        mv = double([R.mode])';
    else
        mv = repmat(m, n, 1);
    end
    T = struct('mode', mv, 'rep', double([R.rep])', 'gen', double([R.gen])', ...
        'count', double([R.count])', 'ok', logical([R.ok])', 'r', double([R.r])', ...
        'xn', double([R.xn])', 'beta', double([R.beta_SI])', 'alpha', double([R.alpha_SI])', ...
        'pLost', double([R.pLostReadout])', 'pExit', double([R.pExit])', ...
        'lost', double([R.lostReadout])');
    if isfield(R, 'alpha_p'), T.alpha_p = double([R.alpha_p])'; else, T.alpha_p = nan(n,1); end
    if isfield(R, 'is_weak'), T.is_weak = logical([R.is_weak])'; else, T.is_weak = true(n,1); end
end

function T = subset(T, keep)
    f = fieldnames(T);
    for i = 1:numel(f)
        if size(T.(f{i}), 1) == numel(keep), T.(f{i}) = T.(f{i})(keep, :); end
    end
end


% =========================================================================
function do_control(C, O, modeNames)
% Check: members rerun under a new seed on 10 replicates per architecture; the rows
% (architecture, replicate, generation, count, r) must equal the census's as sets, and the
% pooled weighted share of draws ending cooperation must lie inside the census's 95% interval on
% the same replicates
    S = load(fullfile(O, 'control_member_10reps.mat'));
    N = flatten(S.W);
    sub = unique(N.rep)';
    c = subset(C, ismember(C.rep, sub) & C.ok);
    n = subset(N, N.ok);
    fprintf('[CONTROL] member branch, new seed %d, replicates %s, all four architectures\n', ...
        S.cfg.rngSeed, mat2str(sub));
    fprintf('  snapshot generations: %s\n', strjoin(cellfun(@mat2str, S.snapGens, 'UniformOutput', false), ' | '));
    % deterministic quantities
    ca = subset(C, ismember(C.rep, sub));
    keyN = [N.mode N.rep N.gen N.count N.r]; keyC = [ca.mode ca.rep ca.gen ca.count ca.r];
    detSame = isequal(sortrows(keyN), sortrows(keyC));
    fprintf('  member rows: new %d, census %d; (arch, rep, gen, count, r) identical as sets: %d\n', ...
        numel(N.r), numel(ca.r), detSame);
    pN = wmean(n.pLost, n.count); pC = wmean(c.pLost, c.count);
    [lo, hi] = boot_rate(c, 2000, 20261016 + 9);
    pass = pN >= lo && pN <= hi;
    fprintf('  pooled clone-weighted P(lose cooperation): new %.4f | census %.4f [%.4f, %.4f] (ok rows %d / %d)\n', ...
        pN, pC, lo, hi, numel(n.r), numel(c.r));
    for m = 1:4
        nm = subset(n, n.mode == m); cm = subset(c, c.mode == m);
        [l2, h2] = boot_rate(cm, 2000, 20261016 + 9 + m);
        fprintf('    %-20s new %.4f | census %.4f [%.4f, %.4f]  n %d\n', modeNames{m}, ...
            wmean(nm.pLost, nm.count), wmean(cm.pLost, cm.count), l2, h2, numel(cm.r));
    end
    fprintf('  pExit pooled: new %.4f | census %.4f\n', wmean(n.pExit, n.count), wmean(c.pExit, c.count));
    fprintf('  wall %.1f min on 2 workers; job seconds sum %.0f\n', S.wallSec/60, sum(S.jobSec));
    if pass && detSame, fprintf('CONTROL: PASS\n'); else, fprintf('CONTROL: FAIL\n'); end
end


% =========================================================================
function do_full(C, O, modeNames)
    S = load(fullfile(O, 'full_cooperate_only.mat'));
    N = flatten(S.W);
    q = subset(N, N.ok);                 % cooperators that do not retaliate
    repsRun = unique([N.mode N.rep], 'rows');
    % members restricted to the replicates the other group's run covered (all 100 by default)
    inRun = ismember([C.mode C.rep], repsRun, 'rows');
    mbr = subset(C, C.ok & inRun);
    K = size(q.lost, 2);
    fprintf('[FULL] cooperate_only seed %d (member census seed %d), K = %d, replicates run per architecture %s\n', ...
        S.cfg.rngSeed, C.seed, K, mat2str(accumarray(repsRun(:,1), 1, [4 1])'));
    fprintf('  snapshot generations: %s\n', strjoin(cellfun(@mat2str, S.snapGens, 'UniformOutput', false), ' | '));
    fprintf('  cooperate_only rows %d (kept %d; r < 1e-6 %d; ||x*|| out %d; is_weak false among kept %d)\n', ...
        numel(N.r), numel(q.r), sum(N.r < 1e-6), sum(N.r >= 1e-6 & ~N.ok), sum(~q.is_weak));
    fprintf('  wall %.1f min on 2 workers\n', S.wallSec/60);

    % Ratio R of the share of draws ending cooperation in the others to that in members, pooled
    % and per architecture
    [R, Rlo, Rhi, pM, pQ, ciM, ciQ] = boot_ratio(mbr, q, 2000, S.cfg.rngSeed + 9);
    fprintf('\nR = P(lose | cooperate_only) / P(lose | member), pooled clone-weighted\n');
    fprintf('  P member %.4f [%.4f, %.4f] (%d wirings, %d draws) | P cooperate_only %.4f [%.4f, %.4f] (%d wirings, %d draws)\n', ...
        pM, ciM(1), ciM(2), numel(mbr.r), numel(mbr.r)*K, pQ, ciQ(1), ciQ(2), numel(q.r), numel(q.r)*K);
    fprintf('  R = %.3f [%.3f, %.3f]\n', R, Rlo, Rhi);
    lower = false(4,1);
    for m = 1:4
        a = subset(mbr, mbr.mode == m); b = subset(q, q.mode == m);
        [Rm, l, h, pa, pb, ca, cb] = boot_ratio(a, b, 2000, S.cfg.rngSeed + 9 + m);
        lower(m) = pb < pa;
        fprintf('  %-20s member %.4f [%.4f, %.4f] n %5d (%6d draws) | cooperate_only %.4f [%.4f, %.4f] n %5d (%6d draws) | ratio %.3f [%.3f, %.3f] | cooperate_only lower %d\n', ...
            modeNames{m}, pa, ca(1), ca(2), numel(a.r), numel(a.r)*K, pb, cb(1), cb(2), numel(b.r), numel(b.r)*K, Rm, l, h, lower(m));
    end

    % Weighted medians [IQR] by group of r (the cooperation margin), beta_SI, r / beta_SI,
    % alpha_SI and ||x*||, then the shares within fifths of r / beta_SI
    fprintf('\nAlso reported\n');
    mbr.rb = mbr.r ./ mbr.beta; q.rb = q.r ./ q.beta;
    fprintf('  clone-weighted median [IQR] (n finite), by group\n');
    fields = {'r', 'beta', 'rb', 'alpha', 'xn'};
    labels = {'r (margin)', 'beta_SI', 'r / beta_SI', 'alpha_SI', '||x*||'};
    for i = 1:numel(fields)
        fprintf('  %-12s', labels{i});
        for g = 1:2
            if g == 1, T = mbr; nm = 'member'; else, T = q; nm = 'cooperate_only'; end
            st = wstats(T.(fields{i}), T.count);
            fprintf('  %s %.4g [%.4g, %.4g] (%d)', nm, st(1), st(2), st(3), st(6));
        end
        fprintf('\n');
        for m = 1:4
            a = subset(mbr, mbr.mode == m); b = subset(q, q.mode == m);
            sa = wstats(a.(fields{i}), a.count); sb = wstats(b.(fields{i}), b.count);
            fprintf('      %-20s member %.4g [%.4g, %.4g] | cooperate_only %.4g [%.4g, %.4g]\n', ...
                modeNames{m}, sa(1), sa(2), sa(3), sb(1), sb(2), sb(3));
        end
    end
    fprintf('  cooperate_only alpha_p (membership_test): %s\n', mat2str(wstats(q.alpha_p, q.count), 4));

    % Share ending cooperation within fifths of r / beta_SI; edges are the weighted quintiles
    % of both groups pooled
    allrb = [mbr.rb; q.rb]; allw = [mbr.count; q.count];
    ed = [-Inf, arrayfun(@(p) wquantile(allrb(isfinite(allrb)), allw(isfinite(allrb)), p), 0.2:0.2:0.8), Inf];
    fprintf('  P(lose cooperation) by quintile of r/beta_SI (edges from both groups pooled, clone-weighted): %s\n', mat2str(ed(2:end-1), 4));
    for b = 1:5
        ia = mbr.rb >= ed(b) & mbr.rb < ed(b+1); ib = q.rb >= ed(b) & q.rb < ed(b+1);
        fprintf('    bin %d  member %.4f (n %5d) | cooperate_only %.4f (n %5d)\n', b, ...
            wmean(mbr.pLost(ia), mbr.count(ia)), sum(ia), wmean(q.pLost(ib), q.count(ib)), sum(ib));
    end
    fprintf('    beta_SI undefined (theorem invalid): member %d, cooperate_only %d rows, excluded from the bins\n', ...
        sum(~isfinite(mbr.rb)), sum(~isfinite(q.rb)));
    fprintf('  split by whether a one-mutation mutant of a member is itself a member: not reported (membership_test tests invasion, not the mutant''s own membership)\n');
end


% =========================================================================
function p = wmean(v, w)
    ok = isfinite(v) & w > 0;
    p = sum(v(ok) .* w(ok)) / sum(w(ok));
end

function [lo, hi] = boot_rate(T, B, seed)
    [num, den, arch] = rep_sums(T);
    rs = rng; rng(seed, 'twister');
    bs = nan(B,1);
    for b = 1:B
        idx = strat_idx(arch);
        bs(b) = sum(num(idx)) / sum(den(idx));
    end
    rng(rs);
    lo = prctile(bs, 2.5); hi = prctile(bs, 97.5);
end

function [R, lo, hi, pA, pB, ciA, ciB] = boot_ratio(A, Bt, B, seed)
% Ratio pB / pA of weighted shares, with intervals over replicates resampled within
% architecture, jointly for both groups (the same population supplies both).
    keys = unique([A.mode A.rep; Bt.mode Bt.rep], 'rows');
    [numA, denA] = sums_on(A, keys); [numB, denB] = sums_on(Bt, keys);
    pA = sum(numA) / sum(denA); pB = sum(numB) / sum(denB); R = pB / pA;
    rs = rng; rng(seed, 'twister');
    bs = nan(B,3);
    for b = 1:B
        idx = strat_idx(keys(:,1));
        a = sum(numA(idx)) / sum(denA(idx)); c = sum(numB(idx)) / sum(denB(idx));
        bs(b,:) = [c / a, a, c];
    end
    rng(rs);
    lo = prctile(bs(:,1), 2.5); hi = prctile(bs(:,1), 97.5);
    ciA = prctile(bs(:,2), [2.5 97.5]); ciB = prctile(bs(:,3), [2.5 97.5]);
end

function [num, den] = sums_on(T, keys)
    [~, loc] = ismember([T.mode T.rep], keys, 'rows');
    num = accumarray(loc, T.count .* T.pLost, [size(keys,1) 1]);
    den = accumarray(loc, T.count, [size(keys,1) 1]);
end

function [num, den, arch] = rep_sums(T)
    keys = unique([T.mode T.rep], 'rows');
    [num, den] = sums_on(T, keys); arch = keys(:,1);
end

function S = wstats(v, w)
    ok = isfinite(v) & isfinite(w) & w > 0;
    v = v(ok); w = w(ok);
    if isempty(v), S = nan(1,6); return; end
    S = [wquantile(v,w,0.5), wquantile(v,w,0.25), wquantile(v,w,0.75), min(v), max(v), numel(v)];
end
