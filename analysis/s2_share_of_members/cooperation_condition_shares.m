function cooperation_condition_shares()
% COOPERATION_CONDITION_SHARES  How many wirings cooperate with their clone, and how many of those
% also retaliate? The share meeting the cooperation condition (C1) and the share of members (C1
% and the retaliation condition, C2), at generations 1 and 500, in the four signaling architectures.
%
% Produces : at generation 1, the (C1) shares and the share of (C1) wirings that fail (C2) (SI S2.3,
%            "Most random wirings that cooperate with their clone do not retaliate")
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures)
% Outputs  : publication/css_c1_c2_split.mat   (per architecture, at generations 1 and 500: the (C1)
%              share `c1_*`, the share of members `vm_*`, the share of (C1) wirings that fail (C2)
%              `lost_*`; success (%, final cooperation at least 0.95); sample sizes;
%              scope_checks_gen1 and scope_checks_summary read only the generation-1 fields); printed tables of
%              the same
% Run      : cooperation_condition_shares   (about 10 min)
% Seed     : generation 1: rng(2026) (the default generator) before each architecture's 10,000
%            random wirings, weights N(0, 0.1^2), so all four test the same wirings; generation
%            500: each architecture is sampled with its own seeded draw, rng(42) before its 30
%            wirings per replicate (randperm); membership_test seeds each wiring itself
%
% One call of the membership test (membership_test, 600 steps, the simulate_pair route) gives both
% results: field is_weak is (C1), is_css is membership. The generation-500 shares here are on
% this file's own sample: gen500_share_of_members (Table S2) sets rng(42) once before all four
% architectures, so the two samples agree for graded separate only. The generation-500 shares
% saved here are read by no other script; Table S2 takes its generation-500 values from
% gen500_share_of_members.

    pubDir = fullfile(grn_data_file('publication'));

    archs = {
        struct('name','graded-separate',     'sm','separate',     'hc',Inf, 'file','graded_separate_100reps.mat'), ...
        struct('name','graded-contrast',     'sm','contrast', 'hc',Inf, 'file','graded_contrast_100reps.mat'), ...
        struct('name','saturating-separate', 'sm','separate',     'hc',1,   'file','saturating_separate_100reps.mat'), ...
        struct('name','saturating-contrast', 'sm','contrast', 'hc',1,   'file','saturating_contrast_100reps.mat')
    };
    condition_names = cellfun(@(a) a.name, archs, 'UniformOutput', false)';

    nGen1     = 10000;   % random wirings per architecture at generation 1 (SI S2.3)
    nPerRep   = 30;      % wirings sampled per replicate at generation 500 (Methods)
    sigma_init = 0.1;

    fprintf('\n=== GENERATION 1 (N = %d random wirings, sigma_init = %.1f) ===\n', nGen1, sigma_init);
    [c1_gen1, vm_gen1] = measure_gen1(archs, nGen1, sigma_init);
    report(condition_names, c1_gen1, vm_gen1);

    fprintf('\n=== GENERATION 500 (%d per replicate) ===\n', nPerRep);
    [c1_gen500, vm_gen500, success_empirical] = measure_gen500(archs, pubDir, nPerRep);
    report(condition_names, c1_gen500, vm_gen500);

    % Share of the wirings satisfying (C1) that fail (C2)
    lost_gen1   = (c1_gen1   - vm_gen1)   ./ c1_gen1;
    lost_gen500 = (c1_gen500 - vm_gen500) ./ c1_gen500;

    outFile = fullfile(pubDir, 'css_c1_c2_split.mat');
    save(outFile, 'condition_names', 'success_empirical', ...
                  'c1_gen1', 'vm_gen1', 'lost_gen1', ...
                  'c1_gen500', 'vm_gen500', 'lost_gen500', ...
                  'nGen1', 'nPerRep', 'sigma_init');
    fprintf('\nSaved to: %s\n', outFile);
end


% =========================================================================
function [c1, is_member] = measure_gen1(archs, N, sigma_init)
% Random wirings, each weight drawn from N(0, sigma_init^2), the same N for every architecture,
% tested for membership on the simulate_pair route.
    n = numel(archs);
    c1 = nan(n,1); is_member = nan(n,1);
    for a = 1:n
        params = struct( ...
            'sensingMode',       archs{a}.sm, ...
            'hillCoefficient',   archs{a}.hc, ...
            'discreteOutput',    false, ...
            'michaelisConstant', 1.0, ...
            'degradationRate',   0.1, ...
            'timeStep',          1.0, ...
            'developmentSteps',  600, ...
            'backend',           'simulate_pair');
        rng(2026);                       % same wiring sample per architecture
        nWeak = 0; nCss = 0; t0 = tic;
        for i = 1:N
            r = membership_test(randn(4) * sigma_init, params);
            nWeak = nWeak + double(r.is_weak);
            nCss  = nCss  + double(r.is_css);
        end
        c1(a) = nWeak / N;
        is_member(a) = nCss  / N;
        fprintf('  %-22s C1 %.4f  members %.4f  (%.0fs)\n', archs{a}.name, c1(a), is_member(a), toc(t0));
    end
end


% =========================================================================
function [c1, is_member, success] = measure_gen500(archs, pubDir, nPerRep)
% Final-generation wirings, nPerRep per replicate by randperm, tested for membership on the simulate_pair route
% with the run's own settings; success (%), the share of replicates whose final cooperation
% is at least 0.95.
    n = numel(archs);
    c1 = nan(n,1); is_member = nan(n,1); success = nan(n,1);
    for a = 1:n
        data = load(fullfile(pubDir, archs{a}.file));
        cfg  = data.config;
        success(a) = 100 * sum(cellfun(@(r) r.historyCooperation(end) >= 0.95, data.results)) / numel(data.results);
        params = struct( ...
            'sensingMode',       cfg.sensingMode, ...
            'hillCoefficient',   cfg.hillCoefficient, ...
            'discreteOutput',    cfg.discreteOutput, ...
            'michaelisConstant', cfg.michaelisConstant, ...
            'degradationRate',   cfg.degradationRate, ...
            'timeStep',          cfg.timeStep, ...
            'developmentSteps',  600, ...
            ... % the simulate_pair route for all four architectures, as at
            ... % generation 1: 'auto' would send graded separate alone to the
            ... % linear route, which requires convergence and so judges it by
            ... % a stricter rule (membership_test, the comment on its routes)
            'backend',           'simulate_pair');
        rng(42);                          % subsample seed, reset for each architecture
        nWeak = 0; nCss = 0; nTot = 0; t0 = tic;
        for rep = 1:numel(data.results)
            B = data.results{rep}.finalNetworkB;
            M = data.results{rep}.finalNetworkM;
            nAg = size(B, 3);
            for k = randperm(nAg, min(nPerRep, nAg))
                W = double(M(:,:,k)) .* double(B(:,:,k));
                r = membership_test(W, params);
                nWeak = nWeak + double(r.is_weak);
                nCss  = nCss  + double(r.is_css);
                nTot  = nTot  + 1;
            end
        end
        c1(a) = nWeak / nTot;
        is_member(a) = nCss  / nTot;
        fprintf('  %-22s C1 %.4f  members %.4f  (n=%d, %.0fs)\n', ...
            archs{a}.name, c1(a), is_member(a), nTot, toc(t0));
    end
end


% =========================================================================
function report(names, c1, is_member)
    fprintf('\n  %-22s %8s %8s %10s %16s\n', ...
        'architecture', 'C1', 'members', 'C1 & ~C2', 'of C1, lost');
    for a = 1:numel(names)
        fprintf('  %-22s %8.4f %8.4f %10.4f %15.1f%%\n', ...
            names{a}, c1(a), is_member(a), c1(a)-is_member(a), 100*(c1(a)-is_member(a))/c1(a));
    end
end
