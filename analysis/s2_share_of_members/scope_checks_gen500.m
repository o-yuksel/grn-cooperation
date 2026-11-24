function scope_checks_gen500(which, nRepLimit)
% SCOPE_CHECKS_GEN500  Which generation-500 members of the 11 settings also resist every single
% topology toggle? Tests the same sampled wirings as the shares of members of Table S2.
%
% Produces : the input that scope_checks_summary reports (Methods, "Success, share of members and
%            statistics"; SI S2.3; Table S2)
% Inputs   : publication/<architecture>_100reps.mat   (run_signaling_architectures)
%            publication/sweep_spatial_N50_<setting>_100reps_<time>.mat   (run_spatial_n50;
%              the newest file of each)
%            publication/group_two_level/group_fission_<arm>_N1000_185reps_500gen.mat
%              (run_group_selection; group selection `fission`, the random-split control
%              `fission_sham`)
%            publication/gen500_c1_c2_split.mat (cooperation_condition_shares_gen500) and
%              group_two_level/group_fission_analysis_N1000.mat (group_selection_arms):
%              the stored member counts and sample sizes, saved beside the result for the check
%              of scope_checks_summary
%            which       settings to run, by index (default 1:11, the order of settings_table)
%            nRepLimit   test only replicates 1..nRepLimit (default Inf, the paper's run); the
%                        sample is still drawn for every replicate, so the tested wirings do not
%                        change
% Outputs  : publication/scope_checks/toggle_gen500_<tag>.mat   (one per setting: per sampled
%              wiring the replicate, the cooperation condition (C1), membership and the toggle
%              result, and the counts); one console line per setting
% Run      : scope_checks_gen500   (3 replicates of all 11 settings take about 5 min; the full run
%              is not timed)
% Seed     : rng(42) (MATLAB's default generator) once before the four signaling architectures, the
%            stream running on through them in order, and again before each group arm; the spatial
%            settings take every wiring; membership_test seeds each wiring itself
%
% Samples as the share-of-members scripts draw them: 30 wirings per replicate by randperm
% (gen500_share_of_members, group_share_of_members), every final wiring in the spatial settings
% (gen500_share_of_members_spatial). Connections and weights are kept apart, so a toggle that
% restores a removed connection gets its stored weight back. toggle_stable tests members at 600 steps.

    if nargin < 1 || isempty(which), which = 1:11; end
    if nargin < 2 || isempty(nRepLimit), nRepLimit = Inf; end
    pub = fullfile(grn_data_file('publication'));
    outDir = fullfile(pub, 'scope_checks');
    if ~isfolder(outDir), mkdir(outDir); end

    S = settings_table(pub);
    for s = which(:)'
        smp = draw_sample(pub, S(s));
        t0 = tic;
        R = test_membership_sample(smp, nRepLimit);
        R.label = S(s).label; R.tag = S(s).tag; R.family = S(s).family;
        R.pub_nvm = S(s).pub_nvm; R.pub_n = S(s).pub_n; R.pub_occ_printed = S(s).pub_occ;
        R.nRepLimit = nRepLimit; R.seconds = toc(t0);
        R.succ = smp.succ;
        save(fullfile(outDir, sprintf('toggle_gen500_%s.mat', strrep(S(s).tag, 'fission_sham', 'random_split'))), '-struct', 'R');
        fprintf(['%-20s n %5d  members %4d (paper %4d)  toggle-fail %3d  ' ...
                 'share %.4f -> %.4f  %.0fs\n'], S(s).label, R.n, R.nMem, S(s).pub_nvm, ...
                 R.nFail, R.nMem / R.n, (R.nMem - R.nFail) / R.n, R.seconds);
    end
end


% =========================================================================
function S = settings_table(pub)
% The 11 settings with the stored member count, sample size and share of members (in %, one
% decimal, as Table S2 prints it): the nine signaling and spatial settings from
% gen500_c1_c2_split.mat, the two group arms from group_fission_analysis_N1000.mat.
    P = load(fullfile(pub, 'gen500_c1_c2_split.mat'), 'per_replicate');
    G = load(fullfile(pub, 'group_two_level', 'group_fission_analysis_N1000.mat'), 'A');
    nv = [arrayfun(@(r) r.nVm, P.per_replicate(1:9)), sum(G.A.fission.occ.nVm), sum(G.A.fission_sham.occ.nVm)];
    nn = [arrayfun(@(r) r.ntot, P.per_replicate(1:9)), sum(G.A.fission.occ.nSampled), sum(G.A.fission_sham.occ.nSampled)];
    S = struct('family', {}, 'label', {}, 'tag', {}, 'pub_nvm', {}, 'pub_n', {}, 'pub_occ', {});
    add = @(k, f, l, t) struct('family', f, 'label', l, 'tag', t, ...
        'pub_nvm', nv(k), 'pub_n', nn(k), 'pub_occ', round(100 * nv(k) / nn(k), 1));
    S(1)  = add(1,  'signaling', 'graded separate',     'graded_separate');
    S(2)  = add(2,  'signaling', 'graded contrast',     'graded_contrast');
    S(3)  = add(3,  'signaling', 'saturating separate', 'saturating_separate');
    S(4)  = add(4,  'signaling', 'saturating contrast', 'saturating_contrast');
    S(5)  = add(5,  'spatial',   'well-mixed',          'wellmixed');
    S(6)  = add(6,  'spatial',   'ring k=1',            'ring_k1');
    S(7)  = add(7,  'spatial',   'ring k=2',            'ring_k2');
    S(8)  = add(8,  'spatial',   'ring k=5',            'ring_k5');
    S(9)  = add(9,  'spatial',   'ring k=10',           'ring_k10');
    S(10) = add(10, 'group',     'group selection',     'fission');
    S(11) = add(11, 'group',     'random splits',       'fission_sham');
end


% =========================================================================
function smp = draw_sample(pub, s)
% One setting's sample: M and B stacks (4x4xn), the replicate of each wiring, success per replicate
% (final cooperation >= 0.95), and the membership-test settings.
    switch s.family
        case 'signaling'
            pats = {'graded_separate_100reps.mat', 'graded_contrast_100reps.mat', ...
                    'saturating_separate_100reps.mat', 'saturating_contrast_100reps.mat'};
            rng(42);
            for m = 1:4
                d = load_latest(pub, pats{m});
                [Ms, Bs, rp] = subsample(d.results, 30);
                if strcmp(pats{m}, [s.tag '_100reps.mat']), break; end
            end
            cfg = d.config;
            succ = cellfun(@(x) double(x.historyCooperation(end)) >= 0.95, d.results(:));
        case 'spatial'
            d = load_latest(pub, sprintf('sweep_spatial_N50_%s_100reps_*.mat', s.tag));
            [Ms, Bs, rp] = subsample(d.results, Inf);
            cfg = d.config;
            succ = cellfun(@(x) double(x.historyCooperation(end)) >= 0.95, d.results(:));
        case 'group'
            f = fullfile(pub, 'group_two_level', ...
                sprintf('group_fission_%s_N1000_185reps_500gen.mat', strrep(s.tag, 'fission_sham', 'random_split')));
            d = load(f, 'cfg', 'res', 'finalCoop');
            rng(42);
            [Ms, Bs, rp] = subsample(d.res, 30);
            cfg = d.cfg;
            succ = double(d.finalCoop(:)) >= 0.95;
    end
    smp = struct('M', Ms, 'B', Bs, 'rep', rp, 'succ', succ, 'params', membership_settings(cfg));
end

function d = load_latest(pub, pat)
    f = dir(fullfile(pub, pat));
    assert(~isempty(f), 'No file matching %s', pat);
    [~, k] = max([f.datenum]);
    d = load(fullfile(pub, f(k).name));
end

function [Ms, Bs, rp] = subsample(results, nSub)
% One randperm per replicate in replicate order (nSub = Inf takes every wiring in index order).
    Ms = zeros(4, 4, 0); Bs = Ms; rp = zeros(0, 1);
    for r = 1:numel(results)
        M = results{r}.finalNetworkM; B = results{r}.finalNetworkB;
        nAg = size(B, 3);
        if isinf(nSub), pick = 1:nAg; else, pick = randperm(nAg, min(nSub, nAg)); end
        Ms = cat(3, Ms, double(M(:, :, pick)));
        Bs = cat(3, Bs, double(B(:, :, pick)));
        rp = [rp; r * ones(numel(pick), 1)]; %#ok<AGROW>
    end
end

% =========================================================================
function R = test_membership_sample(smp, nRepLimit)
% The membership test on every sampled wiring ((C1) and membership), toggle_stable on every member.
    keep = smp.rep <= nRepLimit;
    idx = find(keep);
    n = numel(idx);
    c1 = false(n, 1); is_member = false(n, 1); tog = true(n, 1); cls = zeros(n, 2);
    for q = 1:n
        k = idx(q);
        M = smp.M(:, :, k); B = smp.B(:, :, k);
        e = membership_test(M .* B, smp.params);
        c1(q) = e.is_weak; is_member(q) = e.is_css;
        if is_member(q)
            [tog(q), cls(q, :)] = toggle_stable(M, B, smp.params, 600);
        end
    end
    fail = is_member & ~tog;
    R = struct('rep', smp.rep(idx), 'is_weak', c1, 'is_css', is_member, 'toggle_ok', tog, ...
        'toggle_cls', cls, 'n', n, 'nMem', sum(is_member), 'nFail', sum(fail), ...
        'nC1', sum(c1), 'nPartnerFail', sum(fail & cls(:, 1) > 0), ...
        'nOtherOnlyFail', sum(fail & cls(:, 1) == 0), ...
        'invPartner', sum(cls(is_member, 1)), 'invOther', sum(cls(is_member, 2)), ...
        'scoredPartner', 4 * sum(is_member), 'scoredOther', 12 * sum(is_member));
end
