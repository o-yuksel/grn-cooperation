function out = share_of_members_vs_success(opts)
% SHARE_OF_MEMBERS_VS_SUCCESS  How closely does the share of members track success over the 11
% settings, and more closely than the share meeting only the cooperation condition (C1)? Also counts
% replicates whose success and sampled members disagree, tests the two group arms against the line
% through the nine other settings, and gives intervals on success.
%
% Produces : the correlation r over the 11 settings (main text, "Mutual sensing lets gene networks
%            escape the dilemma"; Fig. 3); the rows of Table S2 (printed as ROW); SI S2.3 (the
%            (C1) comparison, replicates without a member); SI S2.7 (the group arms); Table S3
% Inputs   : publication/gen500_c1_c2_split.mat   (cooperation_condition_shares_gen500;
%              per-replicate counts of the nine signaling and spatial settings)
%            publication/group_two_level/group_fission_analysis_N1000.mat
%              (group_selection_arms; group selection, the random-split control and the undivided
%              arm, and the stored nine-setting line)
%            opts   optional struct: write (default true; false prints only)
% Outputs  : publication/share_of_members_checks_11settings.mat (out) and
%            share_of_members_checks_11settings_report.txt (the lines this run prints)
% Run      : share_of_members_vs_success   (about 1 min)
% Seed     : rng(31415) (MATLAB's default generator) before each bootstrap over settings of the
%            correlation difference; rng(20261007) before the bootstrap over runs of the
%            differences between arms (10,000 resamples each)
%
% The 11 settings: the nine signaling and spatial settings and the two group arms; the undivided
% arm is reported outside them. The run stops unless the line through the nine settings equals the
% line stored with the group arms (A.line).
    if nargin < 1, opts = struct(); end
    doWrite = ~isfield(opts, 'write') || opts.write;
    pub = fullfile(grn_data_file('publication'));
    S = load(fullfile(pub, 'gen500_c1_c2_split.mat'));
    G = load(fullfile(pub, 'group_two_level', 'group_fission_analysis_N1000.mat'));
    A = G.A;

    old = from_split(S.per_replicate(1:9));

    new = [old, from_arm(A.fission), from_arm(A.fission_sham)];
    lab = cellstr(S.labels);
    out.labels = [lab(1:9); {'group selection'; 'random splits'}];
    out.new = summarize(new);
    out.four = summarize(old(1:4));
    out.occ_ci = [S.pVm_all_ci(1:9, :); A.fission.occ_ci(:)'; A.fission_sham.occ_ci(:)'];
    out.line = line_check(S, A);
    out.diff = arm_differences(A, S);
    out.cp = success_cp([new, from_arm(A.control)]);
    out.pc1_undivided = mean(double(A.control.occ.nC1(:)) ./ double(A.control.occ.nSampled(:)));

    print_report(1, out);
    if doWrite
        save(fullfile(pub, 'share_of_members_checks_11settings.mat'), 'out');
        fid = fopen(fullfile(pub, 'share_of_members_checks_11settings_report.txt'), 'w');
        cleaner = onCleanup(@() fclose(fid));
        print_report(fid, out);
    end
end


function R = from_split(P)
% Per-replicate counts (members, sampled wirings, (C1)) and success of gen500_c1_c2_split's
% settings, rebuilt from the stored shares and checked against the stored totals.
    R = struct('n_members_rep', {}, 'nt', {}, 'nc1', {}, 'succ', {});
    for i = 1:numel(P)
        nt = double(P(i).nper(:));
        n_members_rep = round(double(P(i).pVm(:)) .* nt);
        nc1 = round(double(P(i).pC1(:)) .* nt);
        assert(sum(n_members_rep) == P(i).nVm && sum(nc1) == P(i).nC1, 'count reconstruction, setting %d', i);
        R(i) = struct('n_members_rep', n_members_rep, 'nt', nt, 'nc1', nc1, 'succ', logical(P(i).succ(:)));
    end
end


function r = from_arm(v)
% Per-replicate counts of one group arm (or the undivided control); success is final cooperation
% >= 0.95, checked against the arm's stored count.
    o = v.occ;
    succ = double(o.finalCoop(:)) >= 0.95;
    assert(sum(succ) == v.k, 'success count of an arm does not match its analysis file');
    r = struct('n_members_rep', double(o.nVm(:)), 'nt', double(o.nSampled(:)), ...
               'nc1', double(o.nC1(:)), 'succ', succ);
end


function s = summarize(R)
% Over settings: the share of members (pooled ratio), the mean per-replicate (C1) share, success
% (%), both correlations with success and their difference; over replicates: successful ones
% without a sampled member, failed ones with one.
    n = numel(R);
    s.occ = arrayfun(@(r) sum(r.n_members_rep) / sum(r.nt), R)';
    s.pc1 = arrayfun(@(r) mean(r.nc1 ./ r.nt), R)';
    s.suc = arrayfun(@(r) 100 * mean(r.succ), R)';
    s.r_occ = corr(s.occ, s.suc);
    s.r_c1 = corr(s.pc1, s.suc);
    s.gap = s.r_occ - s.r_c1;
    s.gap_ci = gap_bootstrap(s.occ, s.pc1, s.suc);
    succ = vertcat(R.succ);
    hasm = vertcat(R.n_members_rep) > 0;
    s.n_succ = sum(succ);  s.succ_nomember = sum(succ & ~hasm);
    s.n_fail = sum(~succ); s.fail_member = sum(~succ & hasm);
    s.n = n;
end


function ci = gap_bootstrap(share, pc1, suc)
% Percentile bootstrap over settings of the correlation difference, built as paired_r_gap in
% cooperation_condition_shares_gen500 (resamples with fewer than 3 distinct settings skipped).
    rng(31415);
    nB = 10000; n = numel(suc); v = nan(nB, 1);
    for b = 1:nB
        k = randi(n, n, 1);
        if numel(unique(k)) < 3, continue; end
        v(b) = corr(share(k), suc(k)) - corr(pc1(k), suc(k));
    end
    ci = prctile(v(~isnan(v)), [2.5 97.5]);
end


function L = line_check(S, A)
% Least-squares line of success (%) on the share of members (%, rounded to one decimal) over the
% nine signaling and spatial settings, each arm's observed success against its 95% prediction
% interval, and the slope between the two group arms; stops unless the line equals A.line.
    x = round(100 * S.pVm_all(1:9), 1); y = double(S.success_rate(1:9)); x = x(:); y = y(:);
    p = polyfit(x, y, 1); L.slope = p(1); L.intercept = p(2);
    res = y - polyval(p, x); n = numel(x); L.residSD = sqrt(sum(res.^2) / (n - 2));
    t = tinv(0.975, n - 2); sxx = sum((x - mean(x)).^2);
    arms = {A.fission, A.fission_sham, A.control};
    for i = 1:3
        x0 = 100 * arms{i}.occ.occupancy;
        L.pred(i) = polyval(p, x0);
        h = t * L.residSD * sqrt(1 + 1/n + (x0 - mean(x)).^2 / sxx);
        L.pi(i, :) = [L.pred(i) - h, L.pred(i) + h];
        L.obs(i) = 100 * arms{i}.rate;
    end
    L.slope_between_arms = (L.obs(1) - L.obs(2)) / (100 * (A.fission.occ.occupancy - A.fission_sham.occ.occupancy));
    if abs(L.slope - A.line.slope) > 1e-6 || abs(L.intercept - A.line.intercept) > 1e-6
        error('the nine-setting line differs from the line stored in group_fission_analysis_N1000.mat (A.line)');
    end
end


function D = arm_differences(A, S)
% Differences in the share of members between arms (pp), with a bootstrap over runs within each
% arm (pooled ratio, as for the share itself), and rank-sum tests on the per-run shares.
    rng(20261007); B = 10000;
    arms = {A.fission, A.fission_sham, A.control};
    for i = 1:3
        n_members_rep = double(arms{i}.occ.nVm(:)); nt = double(arms{i}.occ.nSampled(:));
        k = randi(numel(n_members_rep), numel(n_members_rep), B);
        bs{i} = sum(n_members_rep(k), 1) ./ sum(nt(k), 1);  %#ok<AGROW>
        per{i} = n_members_rep ./ nt;                        %#ok<AGROW>
        pt(i) = sum(n_members_rep) / sum(nt);                %#ok<AGROW>
    end
    D.gs_rs = 100 * (pt(1) - pt(2)); D.gs_rs_ci = 100 * prctile(bs{1} - bs{2}, [2.5 97.5]);
    D.gs_rs_p = ranksum(per{1}, per{2});
    D.gs_un = 100 * (pt(1) - pt(3)); D.gs_un_ci = 100 * prctile(bs{1} - bs{3}, [2.5 97.5]);
    D.gs_un_p = ranksum(per{1}, per{3});
    D.un_vs_published_p = ranksum(per{3}, double(S.per_replicate(4).pVm(:)));
end


function C = success_cp(R)
% Clopper-Pearson exact 95% intervals on success, in percent: [k, n, lower, upper].
    C = zeros(numel(R), 4);
    for i = 1:numel(R)
        k = sum(R(i).succ); n = numel(R(i).succ);
        lo = 0; hi = 1;
        if k > 0, lo = betainv(0.025, k, n - k + 1); end
        if k < n, hi = betainv(0.975, k + 1, n - k); end
        C(i, :) = [k, n, 100 * lo, 100 * hi];
    end
end


function print_report(fid, o)
% One line per quantity; the ROW lines are the rows of Table S2.
    s = o.new; f = o.four; L = o.line; D = o.diff;
    fprintf(fid, 'R_SHARE_11 %+.3f  R_C1_11 %+.3f  GAP_11 %+.3f  GAP_CI_11 [%+.3f, %+.3f]\n', s.r_occ, s.r_c1, s.gap, s.gap_ci);
    fprintf(fid, 'SUCC_NOMEMBER_11 %d of %d (%.1f%%)  FAIL_MEMBER_11 %d of %d (%.1f%%)\n', ...
        s.succ_nomember, s.n_succ, 100 * s.succ_nomember / s.n_succ, s.fail_member, s.n_fail, 100 * s.fail_member / s.n_fail);
    fprintf(fid, 'FOUR_FAIL_NONE %d of %d\n', f.n_fail - f.fail_member, f.n_fail);
    fprintf(fid, 'LINE success = %.3f + %.4f x share of members (%%), residSD %.2f\n', L.intercept, L.slope, L.residSD);
    nm = {'group selection', 'random splits', 'undivided'};
    for i = 1:3
        fprintf(fid, 'LINE_%s observed %.1f predicted %.1f [%.1f, %.1f] miss %+.1f\n', strrep(nm{i}, ' ', '_'), ...
            L.obs(i), L.pred(i), L.pi(i, :), L.obs(i) - L.pred(i));
    end
    fprintf(fid, 'SLOPE_BETWEEN_ARMS %.2f against %.2f\n', L.slope_between_arms, L.slope);
    fprintf(fid, 'DIFF_GS_RS %+.1f [%+.1f, %+.1f]  P_GS_RS %.2g\n', D.gs_rs, D.gs_rs_ci, D.gs_rs_p);
    fprintf(fid, 'DIFF_GS_UNDIV %+.1f [%+.1f, %+.1f]  P %.2g\n', D.gs_un, D.gs_un_ci, D.gs_un_p);
    fprintf(fid, 'P_UNDIV_VS_SETTING4 %.2f\n', D.un_vs_published_p);
    for i = 1:numel(o.labels)
        fprintf(fid, 'ROW %-22s share %5.1f [%5.1f, %5.1f]  success %d/%d [%.1f, %.1f]  pC1 %.3f\n', o.labels{i}, ...
            100 * s.occ(i), 100 * o.occ_ci(i, :), o.cp(i, 1), o.cp(i, 2), o.cp(i, 3:4), s.pc1(i));
    end
    fprintf(fid, 'ROW %-22s success %d/%d [%.1f, %.1f]  pC1 %.3f\n', 'undivided N=1000', o.cp(end, 1), o.cp(end, 2), o.cp(end, 3:4), o.pc1_undivided);
end
