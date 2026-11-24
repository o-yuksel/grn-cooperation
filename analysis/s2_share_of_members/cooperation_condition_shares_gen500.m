function cooperation_condition_shares_gen500(pubDir)
% COOPERATION_CONDITION_SHARES_GEN500  At generation 500, how much of the share of members is set by
% the share that cooperates with its clone (C1), and how much by the share of those that also
% retaliate (C2)? Measured per replicate in the nine signaling and spatial settings.
%
% Produces : the per-replicate counts, the (C1) share p_C1 and the intervals on the share of members
%            that share_of_members_vs_success and scope_checks_summary read for Table S2, Fig. 3,
%            SI S2.3 and Methods, "Success, share of members and statistics"
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures)
%            publication/sweep_spatial_N50_<setting>_100reps_<time>.mat   (written by
%              run_spatial_n50; the newest file of each setting)
%            publication/css_predicts_sensing_gen500.mat and css_predicts_spatial.mat
%              (gen500_share_of_members, gen500_share_of_members_spatial): the stored member
%              counts, for the check
%            pubDir   optional folder (default: the data root's publication folder)
% Outputs  : publication/gen500_c1_c2_split.mat   (per setting, for all, successful and failed
%              replicates: means and cluster-bootstrap intervals of p_C1, the share of members and
%              p_C2|C1, the share of (C1) wirings that also meet (C2); success (%); the
%              per-replicate counts; the statistics across settings, field `criterion`); printed:
%              agreement with the stored member counts, the tables, the statistics across settings
%              and the variance split
% Run      : cooperation_condition_shares_gen500   (about 8 min)
% Seed     : the default generator: rng(42) once before the four architectures (30 wirings per
%            replicate by randperm; the spatial settings test every wiring); rng(7) for a second
%            sample of distinct genotypes; rng(20260905) and rng(31415) for the bootstraps;
%            membership_test seeds each wiring itself
%
% The nine settings: the four signaling architectures and five spatial settings of Table S2; the two
% group arms of Table S2 come from group_selection_arms. Populations are nearly clonal, so every
% interval resamples replicates. Success: final cooperation at least 0.95.

    if nargin < 1 || isempty(pubDir), pubDir = grn_data_file('publication'); end

    C = setting_table();
    fprintf('pubDir: %s\n', pubDir);

    R = measure_all_settings(pubDir, C);
    ok = report_agreement(pubDir, C, R);
    S = summarise(C, R);
    S.agreement_ok = ok;
    print_table(C, S);
    D = across_settings(S);

    S.criterion = D;
    S.condition = C;
    S.per_replicate = R;
    save(fullfile(pubDir, 'gen500_c1_c2_split.mat'), '-struct', 'S');
    fprintf('\nSaved %s\n', fullfile(pubDir, 'gen500_c1_c2_split.mat'));
end


% =========================================================================
function C = setting_table()
% The nine settings, in the order the sample is drawn.
    C = struct('family', {}, 'label', {}, 'pat', {}, 'nsub', {});
    add = @(f,l,p,ns) struct('family',f,'label',l,'pat',p,'nsub',ns);
    C(1)  = add('signaling','graded-separate',    'graded_separate_100reps.mat',      30);
    C(2)  = add('signaling','graded-contrast',    'graded_contrast_100reps.mat',      30);
    C(3)  = add('signaling','saturating-separate','saturating_separate_100reps.mat',  30);
    C(4)  = add('signaling','saturating-contrast','saturating_contrast_100reps.mat',  30);
    C(5)  = add('spatial','well-mixed','sweep_spatial_N50_wellmixed_100reps_*.mat', Inf);
    C(6)  = add('spatial','ring k=1',  'sweep_spatial_N50_ring_k1_100reps_*.mat',  Inf);
    C(7)  = add('spatial','ring k=2',  'sweep_spatial_N50_ring_k2_100reps_*.mat',  Inf);
    C(8)  = add('spatial','ring k=5',  'sweep_spatial_N50_ring_k5_100reps_*.mat',  Inf);
    C(9)  = add('spatial','ring k=10', 'sweep_spatial_N50_ring_k10_100reps_*.mat', Inf);
end


% =========================================================================
function R = measure_all_settings(pubDir, C)
% Per-replicate shares of (C1) and of members for the nine settings (pVm = share of members).
    R = struct('pC1', {}, 'pVm', {}, 'coop', {}, 'succ', {}, 'nper', {}, ...
               'nVm', {}, 'nC1', {}, 'ntot', {}, 'pC1_d', {}, 'pVm_d', {}, 'nd', {});
    rng(42);   % signaling: one stream across the four architectures
    for i = 1:numel(C)
        [reps, params] = load_setting(pubDir, C(i));
        t0 = tic;
        R(i) = test_membership_replicates(reps, params, C(i).nsub);
        fprintf('%-22s share=%.4f (%d/%d)  C1=%.4f  reps=%d  %.0fs\n', C(i).label, ...
            R(i).nVm / R(i).ntot, R(i).nVm, R(i).ntot, R(i).nC1 / R(i).ntot, ...
            numel(R(i).pVm), toc(t0));
    end
end


% =========================================================================
function [reps, params] = load_setting(pubDir, c)
% Loads one setting's populations (the newest matching file) and builds the membership-test
% settings from that run's own config: 600 steps, the simulate_pair route.
    f = dir(fullfile(pubDir, c.pat));
    assert(~isempty(f), 'No file matching %s', c.pat);
    [~, mx] = max([f.datenum]);
    d = load(fullfile(pubDir, f(mx).name));
    if isfield(d, 'results'), reps = d.results; cfg = d.config;
    else,                     reps = d.res;     cfg = d.cfg;
    end
    params = struct( ...
        'sensingMode',       cfg.sensingMode, ...
        'hillCoefficient',   cfg.hillCoefficient, ...
        'discreteOutput',    cfg.discreteOutput, ...
        'michaelisConstant', cfg.michaelisConstant, ...
        'degradationRate',   cfg.degradationRate, ...
        'timeStep',          cfg.timeStep, ...
        'developmentSteps',  600, ...
        'backend',           'simulate_pair');
end


% =========================================================================
function r = test_membership_replicates(reps, params, nsub)
% Draws the sample of every replicate before testing any for membership (membership_test restores the caller's
% stream, so the order does not change the draws), then the distinct-genotype sample under rng(7)
% inside a saved and restored stream, so the main stream is left as it was.
    n = numel(reps);
    pick = cell(n, 1);
    for k = 1:n
        nAg = size(reps{k}.finalNetworkB, 3);
        if isinf(nsub), pick{k} = 1:nAg;
        else,           pick{k} = randperm(nAg, min(nsub, nAg));
        end
    end
    st = rng();  rng(7);  dpick = distinct_picks(reps, 12);  rng(st);
    r = struct('pC1', nan(n,1), 'pVm', nan(n,1), 'coop', nan(n,1), 'succ', false(n,1), ...
               'nper', zeros(n,1), 'nVm', 0, 'nC1', 0, 'ntot', 0, ...
               'pC1_d', nan(n,1), 'pVm_d', nan(n,1), 'nd', zeros(n,1));
    for k = 1:n
        [c1, is_member] = test_membership_picks(reps{k}, pick{k}, params);
        r.pC1(k)  = mean(c1);  r.pVm(k) = mean(is_member);
        r.nper(k) = numel(c1);
        r.nC1 = r.nC1 + sum(c1);  r.nVm = r.nVm + sum(is_member);  r.ntot = r.ntot + numel(c1);
        r.coop(k) = double(reps{k}.historyCooperation(end));
        r.succ(k) = r.coop(k) >= 0.95;
        [c1d, is_member_d] = test_membership_picks(reps{k}, dpick{k}, params);
        r.pC1_d(k) = mean(c1d);  r.pVm_d(k) = mean(is_member_d);  r.nd(k) = numel(c1d);
    end
end


% =========================================================================
function dpick = distinct_picks(reps, nWant)
% Up to nWant distinct genotypes per replicate (wirings equal after rounding to 6 decimals count
% once): a sample without the copy-number weighting of sampled individuals.
    dpick = cell(numel(reps), 1);
    for k = 1:numel(reps)
        M = reps{k}.finalNetworkM;  B = reps{k}.finalNetworkB;
        nAg = size(B, 3);
        flat = zeros(nAg, 16);
        for j = 1:nAg
            flat(j,:) = reshape(double(M(:,:,j)) .* double(B(:,:,j)), 1, 16);
        end
        [~, ia] = unique(round(flat, 6), 'rows', 'stable');
        dpick{k} = ia(randperm(numel(ia), min(nWant, numel(ia))))';
    end
end


% =========================================================================
function [c1, is_member] = test_membership_picks(rep, pick, params)
% One call of the membership test gives both results: is_weak is (C1), is_css is membership.
    B = rep.finalNetworkB; M = rep.finalNetworkM;
    c1 = false(numel(pick), 1); is_member = c1;
    for q = 1:numel(pick)
        a = pick(q);
        res = membership_test(double(M(:,:,a)) .* double(B(:,:,a)), params);
        c1(q) = res.is_weak;  is_member(q) = res.is_css;
    end
end


% =========================================================================
function ok = report_agreement(pubDir, C, R)
% Agreement check: members and wirings tested here against the counts stored by
% gen500_share_of_members and gen500_share_of_members_spatial. Prints AGREES when every count is
% equal and DIFFERS otherwise; the run continues either way.
    Ss = load(fullfile(pubDir, 'css_predicts_sensing_gen500.mat'), 'n_in_Vm', 'n_subsample');
    Sp = load(fullfile(pubDir, 'css_predicts_spatial.mat'), 'n_in_Vm', 'n_evolved_total');
    storedMem = [Ss.n_in_Vm(:); Sp.n_in_Vm(:)];
    storedN = [Ss.n_subsample(:); Sp.n_evolved_total(:)];
    fprintf('\n--- agreement with the stored member counts ---\n');
    fprintf('%-22s %14s %14s\n', 'setting', 'stored', 'here');
    same = false(numel(C), 1);
    for i = 1:numel(C)
        same(i) = R(i).nVm == storedMem(i) && R(i).ntot == storedN(i);
        fprintf('%-22s %6d/%-7d %6d/%-7d\n', C(i).label, storedMem(i), storedN(i), R(i).nVm, R(i).ntot);
    end
    ok = all(same);
    fprintf('member counts -> %s\n', ternary(ok, 'AGREES', 'DIFFERS'));
end


% =========================================================================
function S = summarise(C, R)
% Condition-level means, cluster-bootstrap intervals over replicates, and the spread of p_C2|C1
% across the replicates of each stratum. Strata: 1 = all, 2 = successful, 3 = failed.
    rng(20260905);
    nB = 10000;
    n = numel(C);
    S = struct();
    names = {'all', 'succ', 'fail'};
    for s = 1:3
        f = names{s};
        S.(['pC1_' f])  = nan(n,1);  S.(['pC1_' f '_ci']) = nan(n,2);
        S.(['pVm_' f])  = nan(n,1);  S.(['pVm_' f '_ci']) = nan(n,2);
        S.(['pC2_' f])  = nan(n,1);  S.(['pC2_' f '_ci']) = nan(n,2);
        S.(['nrep_' f]) = zeros(n,1);
        S.(['pC2_' f '_spread']) = nan(n,5);   % median IQRlo IQRhi min max
        S.(['pC2_' f '_nspread']) = zeros(n,1);
    end
    S.pC1_d_all = nan(n,1); S.pVm_d_all = nan(n,1); S.pC2_d_all = nan(n,1);
    S.pC1_d_succ = nan(n,1); S.pVm_d_succ = nan(n,1); S.pC2_d_succ = nan(n,1);
    for i = 1:n
        m = {true(size(R(i).succ)), R(i).succ, ~R(i).succ};
        for s = 1:3
            S = fill_stratum(S, i, names{s}, R(i), m{s}, nB);
        end
        S.pC1_d_all(i)  = mean(R(i).pC1_d);  S.pVm_d_all(i)  = mean(R(i).pVm_d);
        S.pC2_d_all(i)  = S.pVm_d_all(i) / S.pC1_d_all(i);
        S.pC1_d_succ(i) = mean(R(i).pC1_d(R(i).succ));
        S.pVm_d_succ(i) = mean(R(i).pVm_d(R(i).succ));
        S.pC2_d_succ(i) = S.pVm_d_succ(i) / S.pC1_d_succ(i);
    end
    S.success_rate = arrayfun(@(r) 100 * sum(r.succ) / numel(r.succ), R(:));
    S.labels       = {C.label}';
    S.families     = {C.family}';
end


% =========================================================================
function S = fill_stratum(S, i, f, r, mask, nB)
% One setting, one stratum: means, the ratio of means p_C2|C1, their cluster bootstraps, and the
% spread of the per-replicate p_C2|C1 (median, IQR, range, count) over replicates with p_C1 > 0.
    x = [r.pC1(mask), r.pVm(mask)];
    S.(['nrep_' f])(i) = size(x, 1);
    if isempty(x), return; end
    S.(['pC1_' f])(i) = mean(x(:,1));
    S.(['pVm_' f])(i) = mean(x(:,2));
    S.(['pC2_' f])(i) = mean(x(:,2)) / mean(x(:,1));
    S.(['pC1_' f '_ci'])(i,:) = boot_ci(x, @(z) mean(z(:,1)), nB);
    S.(['pVm_' f '_ci'])(i,:) = boot_ci(x, @(z) mean(z(:,2)), nB);
    S.(['pC2_' f '_ci'])(i,:) = boot_ci(x, @(z) mean(z(:,2)) / mean(z(:,1)), nB);
    v = x(x(:,1) > 0, 2) ./ x(x(:,1) > 0, 1);
    S.(['pC2_' f '_nspread'])(i) = numel(v);
    if ~isempty(v)
        q = quantile(v, [0.25 0.75]);
        S.(['pC2_' f '_spread'])(i,:) = [median(v), q(1), q(2), min(v), max(v)];
    end
end


% =========================================================================
function ci = boot_ci(x, statfun, nB)
% Percentile cluster bootstrap: the resampling unit is the replicate, because populations are
% nearly clonal and the wirings of one population are not independent draws.
    n = size(x, 1);
    if n < 2, ci = [NaN NaN]; return; end
    v = nan(nB, 1);
    for b = 1:nB
        v(b) = statfun(x(randi(n, n, 1), :));
    end
    ci = prctile(v(~isnan(v)), [2.5 97.5]);
end


% =========================================================================
function print_table(C, S)
% The nine-row table per stratum, with counts and the spread of p_C2|C1 across replicates.
    for f = {'all', 'succ', 'fail'}
        s = f{1};
        fprintf('\n--- stratum: %s ---\n', s);
        fprintf('%-22s %5s %5s %8s %8s %8s %26s %20s\n', 'setting', 'suc%', 'nrep', ...
            'p_C1', 'p_mem', 'p_C2|C1', 'p_C2|C1 med [IQR]', 'range');
        for i = 1:numel(C)
            sp = S.(['pC2_' s '_spread'])(i,:);
            fprintf('%-22s %5d %5d %8.4f %8.4f %8.4f   %7.4f [%6.4f,%6.4f] %8.4f-%.4f (n=%d)\n', ...
                C(i).label, S.success_rate(i), S.(['nrep_' s])(i), S.(['pC1_' s])(i), ...
                S.(['pVm_' s])(i), S.(['pC2_' s])(i), sp(1), sp(2), sp(3), sp(4), sp(5), ...
                S.(['pC2_' s '_nspread'])(i));
        end
    end
    fprintf('\n--- 95%% cluster-bootstrap CIs, successful stratum ---\n');
    for i = 1:numel(C)
        fprintf('%-22s p_C1 %.4f [%.4f,%.4f]  p_mem %.4f [%.4f,%.4f]  p_C2|C1 %.4f [%.4f,%.4f]\n', ...
            C(i).label, S.pC1_succ(i), S.pC1_succ_ci(i,1), S.pC1_succ_ci(i,2), ...
            S.pVm_succ(i), S.pVm_succ_ci(i,1), S.pVm_succ_ci(i,2), ...
            S.pC2_succ(i), S.pC2_succ_ci(i,1), S.pC2_succ_ci(i,2));
    end
end


% =========================================================================
function D = across_settings(S)
% Statistics across the settings: the smallest p_C1 and the fold ranges of p_C1 and p_C2|C1 among
% successful replicates, the correlations with success, and the decomposition of the log share of
% members; the same on the distinct-genotype sample; and whether the share of members tracks
% success more closely than p_C1 does.
    suc = S.success_rate;
    fold = @(v) max(v) / min(v);
    D.pC1_succ_min      = min(S.pC1_succ);
    D.pC1_succ_fold     = fold(S.pC1_succ);
    D.pC2_succ_fold     = fold(S.pC2_succ);
    D.r_pC2_succ        = corr(S.pC2_succ, suc);
    D.r_pC1_all         = corr(S.pC1_all,  suc);
    D.r_pC1_succ        = corr(S.pC1_succ, suc);
    D.r_pVm_all         = corr(S.pVm_all,  suc);
    D.pC2_all_fold      = fold(S.pC2_all);
    D.pC1_all_fold      = fold(S.pC1_all);
    D.decomp_all  = log_decomposition(S.pC1_all,  S.pC2_all,  S.pVm_all);
    D.decomp_succ = log_decomposition(S.pC1_succ, S.pC2_succ, S.pVm_succ);
    % The same quantities on the distinct-genotype sample
    D.dist.pC1_succ_min  = min(S.pC1_d_succ);
    D.dist.pC1_succ_fold = fold(S.pC1_d_succ);
    D.dist.pC2_succ_fold = fold(S.pC2_d_succ);
    D.dist.r_pC2_succ    = corr(S.pC2_d_succ, suc);
    D.dist.r_pC1_all     = corr(S.pC1_d_all,  suc);
    D.dist.r_pVm_all     = corr(S.pVm_d_all,  suc);
    D.dist.pC2_all_fold  = fold(S.pC2_d_all);
    D.dist.decomp_all    = log_decomposition(S.pC1_d_all, S.pC2_d_all, S.pVm_d_all);
    % Does the share of members track success more closely than p_C1 does: the difference of
    % the two correlations over settings
    D.posthoc = paired_r_gap(S.pVm_all, S.pC1_all, suc);
    print_across(D);
end


% =========================================================================
function p = paired_r_gap(share, pc1, suc)
% r(share of members, success) - r(p_C1, success), with a percentile bootstrap over the
% settings (the unit of the correlation); resamples with fewer than 3 distinct settings are
% skipped.
    rng(31415);
    p.r_occ = corr(share, suc);
    p.r_pC1 = corr(pc1, suc);
    p.gap   = p.r_occ - p.r_pC1;
    nB = 10000; n = numel(suc); v = nan(nB, 1);
    for b = 1:nB
        s = randi(n, n, 1);
        if numel(unique(s)) < 3, continue; end
        v(b) = corr(share(s), suc(s)) - corr(pc1(s), suc(s));
    end
    v = v(~isnan(v));
    p.ci = prctile(v, [2.5 97.5]);
    p.frac_positive = mean(v > 0);
end


% =========================================================================
function d = log_decomposition(pC1, pC2, pVm)
% log(share of members) = log p_C1 + log p_C2|C1 exactly, so its variance across settings
% splits into Var(A) + Var(B) + 2 Cov(A, B), with shares summing to 1.
    A = log(pC1);  B = log(pC2);  T = log(pVm);
    d.max_identity_err = max(abs((A + B) - T));
    V = var(T);
    d.var_total  = V;
    d.share_C1   = var(A) / V;
    d.share_C2   = var(B) / V;
    d.share_cov  = 2 * cov_scalar(A, B) / V;
    d.sum_shares = d.share_C1 + d.share_C2 + d.share_cov;
    assert(abs(d.sum_shares - 1) < 1e-9, 'variance shares must sum to 1');
end

function c = cov_scalar(a, b)
    m = cov(a, b);  c = m(1,2);
end


% =========================================================================
function print_across(D)
    fprintf('\n--- across the settings ---\n');
    fprintf('min p_C1 among successful   = %.4f\n', D.pC1_succ_min);
    fprintf('fold-range p_C1 successful  = %.3f\n', D.pC1_succ_fold);
    fprintf('fold-range p_C2|C1 succ     = %.3f\n', D.pC2_succ_fold);
    fprintf('r(p_C2|C1 succ, success)    = %+.4f\n', D.r_pC2_succ);
    fprintf('r(p_C1 all, success)        = %+.4f\n', D.r_pC1_all);
    fprintf('fold-range p_C2|C1 all      = %.3f\n', D.pC2_all_fold);
    fprintf('r(share of members, success) = %+.4f\n', D.r_pVm_all);
    for f = {'all', 'succ'}
        d = D.(['decomp_' f{1}]);
        fprintf('decomposition (%s): Var(log p_mem) = %.4f | C1 %.3f  C2|C1 %.3f  2Cov %.3f\n', ...
            f{1}, d.var_total, d.share_C1, d.share_C2, d.share_cov);
    end
    fprintf('\n--- 12 distinct genotypes per replicate ---\n');
    fprintf('min p_C1 succ %.4f | fold p_C1 succ %.3f | fold p_C2|C1 succ %.3f | r(p_C2|C1 succ) %+.4f\n', ...
        D.dist.pC1_succ_min, D.dist.pC1_succ_fold, D.dist.pC2_succ_fold, D.dist.r_pC2_succ);
    fprintf('r(p_C1 all) %+.4f | fold p_C2|C1 all %.3f | r(p_mem) %+.4f | shares C1 %.3f  C2|C1 %.3f  2Cov %.3f\n', ...
        D.dist.r_pC1_all, D.dist.pC2_all_fold, D.dist.r_pVm_all, ...
        D.dist.decomp_all.share_C1, D.dist.decomp_all.share_C2, D.dist.decomp_all.share_cov);
    p = D.posthoc;
    fprintf('\n--- does the share of members track success more closely than p_C1? ---\n');
    fprintf('r(p_mem) %+.4f  r(p_C1) %+.4f  gap %+.4f  [%+.4f, %+.4f]  positive in %.1f%% of resamples\n', ...
        p.r_occ, p.r_pC1, p.gap, p.ci(1), p.ci(2), 100*p.frac_positive);
end
