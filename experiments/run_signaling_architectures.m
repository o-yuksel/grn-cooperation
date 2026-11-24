function run_signaling_architectures(nReps, nWorkers)
% RUN_SIGNALING_ARCHITECTURES  Evolves 100 populations in each of the four signaling architectures,
% and runs the sweeps of the temptation T and the population size N at saturating contrast.
%
% Produces : the runs behind Fig. 1, Fig. 2, Fig. 3, Fig. 5, Table S2, Table S3 and Table S5;
%            Abstract and Significance; main text, "Mutual sensing lets gene networks escape
%            the dilemma, yet many populations fail to hold cooperation" (the T and N rows of Table S5 are read
%            by tableS5_parameter_sweeps)
% Inputs   : nReps      replicates per setting (default: the paper's 100, and 50 in the N sweep)
%            nWorkers   pool size if no pool is open (default 4)
%            (the standard setting of publication_config, the architecture, T or N changed)
% Outputs  : publication/<architecture>_<n>reps.mat, <architecture> one of graded_separate,
%              graded_contrast, saturating_separate, saturating_contrast (n = 100: the stored
%              files' names; an existing file stops the run)
%            publication/sweep_T_<T>_<n>reps_<time>.mat   (T = 3.5, 4.0, 4.5, 5.5, 5.75, the
%              last named 5.8)
%            publication/sweep_N_<N>_<n>reps_<time>.mat   (N = 50, 100, 200, 500)
% Run      : run_signaling_architectures   (several hours on 4 workers); run_signaling_architectures(3)
%            runs 3 replicates per setting
% Seed     : replicate r runs under rng(12345 + r, 'twister'), the same in every run

if nargin < 1, nReps = []; end
if nargin < 2 || isempty(nWorkers), nWorkers = 4; end
nMain = 100; nPop = 50;
if ~isempty(nReps), nMain = nReps; nPop = nReps; end

resultsDir = fullfile(grn_data_file('publication'));
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

if isempty(gcp('nocreate'))
    parpool('Processes', nWorkers);
end

fprintf('=== EVOLUTION BATCH ===\n');
fprintf('Config: mut=0.005, deg=0.10, s=1.0, T=5.0, N=1000, Hill n=1, contrast\n');
fprintf('Started: %s\n\n', datestr(now));

% The four signaling architectures: graded or saturating output (hillCoefficient Inf or 1),
% separate or contrast sensing (sensingMode 'separate' or 'contrast'); Methods, "Model"
fprintf('=== [1/3] NORMALIZATION FACTORIAL (4 settings) ===\n');

norm_settings = {
    'graded_separate',    struct('hillCoefficient', Inf, 'discreteOutput', false, 'sensingMode', 'separate');
    'graded_contrast',    struct('hillCoefficient', Inf, 'discreteOutput', false, 'sensingMode', 'contrast');
    'saturating_separate', struct('hillCoefficient', 1,  'discreteOutput', false, 'sensingMode', 'separate');
    'saturating_contrast', struct('hillCoefficient', 1,  'discreteOutput', false, 'sensingMode', 'contrast');
};

for i = 1:size(norm_settings, 1)
    label = norm_settings{i, 1};
    overrides = norm_settings{i, 2};
    fprintf('\n--- %s ---\n', label);

    config = publication_config('continuous_output', ...
        'hillCoefficient', overrides.hillCoefficient, ...
        'discreteOutput', overrides.discreteOutput, ...
        'sensingMode', overrides.sensingMode);
    config.numReplicates = nMain;
    filename = sprintf('%s_%dreps.mat', label, nMain);
    assert(~isfile(fullfile(resultsDir, filename)), 'run_signaling_architectures:Exists', ...
        '%s exists; not overwriting', filename);

    tic;
    results = evolve_population(config);
    elapsed = toc;

    coopValues = cellfun(@(r) r.historyCooperation(end), results);
    successRate = mean(coopValues >= 0.95) * 100;
    meanCoop = mean(coopValues) * 100;

    fprintf('  Success: %.0f%%, Mean coop: %.1f%%, Time: %.0f sec\n', successRate, meanCoop, elapsed);

    save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'successRate', 'meanCoop', 'elapsed', '-v7.3');
    fprintf('  Saved: %s\n', filename);
end

% Temptation T at saturating contrast (Table S5); T = 5.0 is saturating_contrast_100reps.mat,
% and T = 2.5 and 3.0 come from run_low_temptation
fprintf('\n=== [2/3] T-SENSITIVITY ===\n');

T_values = [3.5, 4.0, 4.5, 5.5, 5.75];

for i = 1:length(T_values)
    T = T_values(i);
    fprintf('\n--- T=%.1f ---\n', T);

    config = publication_config('continuous_output', 'payoffDC', T);
    config.numReplicates = nMain;

    tic;
    results = evolve_population(config);
    elapsed = toc;

    coopValues = cellfun(@(r) r.historyCooperation(end), results);
    successRate = mean(coopValues >= 0.95) * 100;
    meanCoop = mean(coopValues) * 100;

    fprintf('  Success: %.0f%%, Mean coop: %.1f%%, Time: %.0f sec\n', successRate, meanCoop, elapsed);

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('sweep_T_%.1f_%dreps_%s.mat', T, nMain, timestamp);
    save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'successRate', 'meanCoop', 'elapsed', '-v7.3');
    fprintf('  Saved: %s\n', filename);
end

% Population size N at saturating contrast, 50 replicates each (Table S5); N = 1,000 is
% saturating_contrast_100reps.mat
fprintf('\n=== [3/3] POPSIZE SWEEP ===\n');

N_values = [50, 100, 200, 500];

for i = 1:length(N_values)
    N = N_values(i);
    fprintf('\n--- N=%d ---\n', N);

    config = publication_config('continuous_output', 'populationSize', N);
    config.numReplicates = nPop;

    tic;
    results = evolve_population(config);
    elapsed = toc;

    coopValues = cellfun(@(r) r.historyCooperation(end), results);
    successRate = mean(coopValues >= 0.95) * 100;
    meanCoop = mean(coopValues) * 100;

    fprintf('  Success: %.0f%%, Mean coop: %.1f%%, Time: %.0f sec\n', successRate, meanCoop, elapsed);

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('sweep_N_%d_%dreps_%s.mat', N, nPop, timestamp);
    save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'successRate', 'meanCoop', 'elapsed', '-v7.3');
    fprintf('  Saved: %s\n', filename);
end

fprintf('\n=== ALL EXPERIMENTS COMPLETE: %s ===\n', datestr(now));
end
