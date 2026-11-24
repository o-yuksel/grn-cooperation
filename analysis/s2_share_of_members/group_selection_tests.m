function A = group_selection_tests(files, tag)
% GROUP_SELECTION_TESTS  Success, the share of members and the share of cooperators that do not
% retaliate in the three group-selection arms, with the tests between the arms.
%
% Produces : through group_selection_arms, which calls it on the N = 1,000 runs: the group-selection
%            results of SI S2.7, Table S2 and Fig. 3; Methods, "Game and evolution"; main text,
%            "Mutual sensing lets gene networks escape the dilemma"
% Inputs   : files   optional struct of the arms' run files (finalCoop, seeds, res, cfg), full paths,
%                    fields `fission`, `fission_sham`, `control` (named below); default publication/
%                    group_two_level/group_fission_<arm>_N1000_185reps_500gen.mat (N = 1,000;
%                    written by run_group_selection)
%            tag     names the output (default 'N1000')
%            publication/gen500_c1_c2_split.mat   (cooperation_condition_shares_gen500; the share
%              of members and success of the nine signaling and spatial settings, for the line)
% Outputs  : publication/group_two_level/group_fission_analysis_<tag>.mat   (struct A); prints per
%            arm success with its Wilson 95% interval, the Fisher tests, the splits per replicate,
%            the share of members, and the success a line fitted over nine settings predicts
% Run      : group_selection_arms   (the paper's call: these defaults, with the checks of the run
%            files and the fields it adds)
% Seed     : rng(20261007) (the default generator) for the bootstrap over replicates, set after the
%            share of members is measured under group_share_of_members's own rng(42)
%
% Arms, by stored name: 'fission', group selection (each generation one group, drawn in proportion
% to its mean selection weight, splits and replaces another); 'fission_sham', the random-split
% control (the splitting group drawn uniformly); 'control', the undivided population. The line:
% least-squares success (%) on the share of members (%, one decimal) over the four signaling
% architectures and the five N = 50 spatial settings, with its 95% prediction interval.

    dDir = fullfile(grn_data_file('publication'), 'group_two_level');
    if nargin < 1 || isempty(files)
        f = @(a) fullfile(dDir, sprintf('group_fission_%s_N1000_185reps_500gen.mat', a));
        files = struct('fission', f('fission'), 'fission_sham', f('random_split'), 'control', f('control'));
    end
    if nargin < 2 || isempty(tag), tag = 'N1000'; end
    arms = {'fission', 'fission_sham', 'control'};
    alpha = 0.05; nBoot = 10000;

    for a = 1:3
        d = load(files.(arms{a}), 'finalCoop', 'seeds', 'res', 'cfg');
        r.file = files.(arms{a});
        r.seeds = [d.seeds(1) d.seeds(end)];
        r.n = numel(d.finalCoop);
        r.k = sum(d.finalCoop >= 0.95);
        [r.rate, r.lo, r.hi] = wilson(r.k, r.n);
        % Generations with exactly one fission, per replicate, in the two group arms
        if a < 3
            Ng = d.cfg.groupSize;
            r.events = cellfun(@(x) count_fission_events(x.historyParents, Ng), d.res(:));
        else
            r.events = [];
        end
        clear d
        r.occ = group_share_of_members(files.(arms{a}));
        A.(arms{a}) = r;
    end

    % Two-sided Fisher exact tests on success: (i) group selection against the undivided
    % population, read at alpha = 0.05; (ii) group selection against random splitting, its raw p.
    % Random splitting against the undivided population: the difference only.
    A.p_primary = fisher2(A.fission, A.control);
    A.reject_primary = A.p_primary < alpha;
    A.p_supporting = fisher2(A.fission, A.fission_sham);
    A.diff_i = A.fission.rate - A.control.rate;
    A.diff_ii = A.fission.rate - A.fission_sham.rate;
    A.diff_sham_control = A.fission_sham.rate - A.control.rate;

    % The least-squares line of success (%) on the share of members (%, one decimal) over 9
    % settings, the four signaling architectures and the five N = 50 spatial settings, read from
    % gen500_c1_c2_split.mat, and its 95% prediction interval at each arm's share of members.
    % Bootstrap: 10,000 resamples of replicates within each arm, percentile intervals, for the
    % share of members and the share that cooperate but do not retaliate (C1 without C2).
    S9 = load(grn_data_file('publication/gen500_c1_c2_split.mat'), 'pVm_all', 'success_rate');
    share9 = round(100 * S9.pVm_all(1:9), 1); share9 = share9(:);
    suc9 = S9.success_rate(1:9); suc9 = suc9(:);
    X = [ones(9, 1) share9]; b = X \ suc9;
    s = sqrt(sum((suc9 - X * b).^2) / 7); tq = tinv(0.975, 7);
    rng(20261007);
    boot = struct();
    for a = 1:3
        o = A.(arms{a}).occ; nR = numel(o.nSampled);
        idx = randi(nR, nR, nBoot);
        S = sum(o.nSampled(idx), 1);
        A.(arms{a}).occ_ci = prctile(sum(o.nVm(idx), 1) ./ S, [2.5 97.5]);
        boot.(arms{a}) = sum(o.nC1noC2(idx), 1) ./ S;
        x0 = 100 * o.occupancy;
        se = s * sqrt(1 + [1 x0] / (X' * X) * [1; x0]);
        A.(arms{a}).line_pred = [1 x0] * b;
        A.(arms{a}).line_pi = [1 x0] * b + [-1 1] * tq * se;
    end
    % Two comparisons, stored as fields S1 and S2 (not SI sections). First: the share that cooperate
    % but do not retaliate, group selection minus undivided, with its interval (S1_met is 1 when the
    % interval lies above 0). Second: group-selection success above the upper end of the line's
    % prediction interval (S2_met is 1 then).
    dFC = boot.fission - boot.control;
    A.S1 = [A.fission.occ.shareC1noC2 - A.control.occ.shareC1noC2, prctile(dFC, [2.5 97.5])];
    A.S1_met = A.S1(2) > 0;
    A.S2_met = 100 * A.fission.rate > A.fission.line_pi(2);
    A.line = struct('intercept', b(1), 'slope', b(2), 'r', corr(share9, suc9), 'residSD', s);

    save(fullfile(dDir, sprintf('group_fission_analysis_%s.mat', tag)), 'A');
    print_report(A, arms);
end


function print_report(A, arms)
    fprintf('\n===== GROUP SELECTION BY FISSION: success per arm and tests =====\n');
    for a = 1:3
        r = A.(arms{a});
        fprintf('%-13s success %3d/%3d = %.3f [Wilson %.3f, %.3f]  seeds %d-%d\n', ...
            arms{a}, r.k, r.n, r.rate, r.lo, r.hi, r.seeds);
    end
    fprintf('Comparison 1, fission vs control: diff %+.3f, Fisher p = %.4g, below 0.05: %d\n', ...
        A.diff_i, A.p_primary, A.reject_primary);
    fprintf('Comparison 2, fission vs random splits: diff %+.3f, Fisher p = %.4g (uncorrected)\n', ...
        A.diff_ii, A.p_supporting);
    fprintf('random splits vs control (descriptive): diff %+.3f\n', A.diff_sham_control);
    for a = 1:2
        e = A.(arms{a}).events;
        fprintf('%-13s fission events per replicate: min %d, median %d, max %d\n', arms{a}, min(e), median(e), max(e));
    end
    for a = 1:3
        r = A.(arms{a}); o = r.occ;
        fprintf(['%-13s share of members %.4f (%d/%d) [boot %.4f, %.4f]; line predicts %.1f%% [PI %.1f, %.1f], ' ...
                 'observed %.1f%%; C1 not C2 %.4f of sampled (%.4f of C1)\n'], arms{a}, o.occupancy, ...
            sum(o.nVm), sum(o.nSampled), r.occ_ci, r.line_pred, r.line_pi, 100 * r.rate, ...
            o.shareC1noC2, o.shareC1noC2_givenC1);
    end
    fprintf('Cooperators that do not retaliate (C1 without C2), fission - control: %+.4f [boot %.4f, %.4f]\n', A.S1);
    fprintf('Fission success %.1f%% against the upper end of the line''s 95%% PI, %.1f%%\n', 100 * A.fission.rate, A.fission.line_pi(2));
end


function p = fisher2(x, y)
    [~, p] = fishertest([x.k, x.n - x.k; y.k, y.n - y.k]);
end


% Wilson score 95% interval of a success share k / n.
function [p, lo, hi] = wilson(k, n)
    z = 1.959963984540054;
    p = k / n;
    den = 1 + z^2 / n;
    c = (p + z^2 / (2 * n)) / den;
    h = z / den * sqrt(p * (1 - p) / n + z^2 / (4 * n^2));
    lo = max(0, c - h); hi = min(1, c + h);
end
