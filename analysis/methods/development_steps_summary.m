function T = development_steps_summary()
% DEVELOPMENT_STEPS_SUMMARY  How cooperation success at saturating contrast depends on the
% number of development steps allowed: 30, 50, 100, 200, 300 or 500.
%
% Produces : success per step limit and the test of 300 against 500 steps, printed: SI Extended
%            Materials and Methods, "Developmental integration and convergence"
% Inputs   : publication/devsteps_standard_steps<S>_100reps.mat   (S = 30, 50, 100, 200, 300,
%            500; coopValues; written by run_development_steps)
% Outputs  : printed per step limit: success k/n (last-generation cooperation at least 0.95) with
%            its Clopper-Pearson exact 95% interval and the mean last-generation cooperation, in
%            percent; two-sided Fisher exact p for adjacent limits; returns table T, writes no file
% Run      : development_steps_summary
% Seed     : none (no random draws)
%
% The step limits share the replicate seeds (replicate r under rng(12345 + r), set in
% evolve_population), so they are paired by replicate; the Fisher tests treat them as independent.

    dataDir = fullfile(grn_data_file('publication'));
    steps = [30 50 100 200 300 500];

    n = numel(steps);
    k = zeros(n, 1); N = k; mc = k; lo = k; hi = k;
    for i = 1:n
        d = load(fullfile(dataDir, sprintf('devsteps_standard_steps%d_100reps.mat', steps(i))), ...
                 'coopValues');
        v = double(d.coopValues(:));
        k(i) = sum(v >= 0.95); N(i) = numel(v); mc(i) = 100 * mean(v);
        [~, ci] = binofit(k(i), N(i));          % Clopper-Pearson exact
        lo(i) = 100 * ci(1); hi(i) = 100 * ci(2);
    end

    fprintf('\n=== Developmental duration, saturating contrast (T = 5) ===\n');
    fprintf('%6s  %11s  %17s  %13s\n', 'steps', 'success', '95% CI (CP), %', 'mean coop, %');
    for i = 1:n
        fprintf('%6d  %7d/%-3d  %7.1f - %-7.1f  %13.1f\n', steps(i), k(i), N(i), lo(i), hi(i), mc(i));
    end

    % Adjacent step counts; the 300-500 pair added here is one of them, so unique keeps it once
    pairs = [(1:n-1)' (2:n)'; find(steps == 300) find(steps == 500)];
    pairs = unique(pairs, 'rows', 'stable');
    fprintf('\nFisher exact (two-sided), success counts:\n');
    for j = 1:size(pairs, 1)
        a = pairs(j, 1); b = pairs(j, 2);
        [~, p] = fishertest([k(a) N(a)-k(a); k(b) N(b)-k(b)]);
        fprintf('  %4d vs %-4d  %d/%d vs %d/%d  p = %.3g\n', ...
            steps(a), steps(b), k(a), N(a), k(b), N(b), p);
    end

    T = table(steps(:), k, N, lo, hi, mc, ...
        'VariableNames', {'steps', 'success', 'n', 'ciLoPct', 'ciHiPct', 'meanCoopPct'});
end
