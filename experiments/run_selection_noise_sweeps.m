function run_selection_noise_sweeps(nReps, nWorkers)
% RUN_SELECTION_NOISE_SWEEPS  Evolves 100 saturating-contrast populations at each value of the
% selection power zeta and of the noise in gene expression.
%
% Produces : the selection and noise rows of Table S5 (read by tableS5_parameter_sweeps)
% Inputs   : nReps      replicates per value (default 100, the paper's runs)
%            nWorkers   pool size if no pool is open (default 4)
%            (the standard setting of publication_config with one parameter changed)
% Outputs  : publication/sweep_selection_s<zeta>_<nReps>reps_<time>.mat (zeta = 0.5, 1, 2, 4, 8)
%            and publication/sweep_noise_<sigma>_<nReps>reps_<time>.mat (sigma, the noise's
%            standard deviation, = 0, 0.01, 0.05, 0.10)
% Run      : run_selection_noise_sweeps   (several hours on 4 workers); run_selection_noise_sweeps(3)
%            runs 3 replicates per value
% Seed     : replicate r runs under rng(12345 + r, 'twister') at every value

if nargin < 1 || isempty(nReps), nReps = 100; end
if nargin < 2 || isempty(nWorkers), nWorkers = 4; end

resultsDir = fullfile(grn_data_file('publication'));

if isempty(gcp('nocreate'))
    parpool('Processes', nWorkers);
end

fprintf('=== PARAMETER SWEEPS ===\n');
fprintf('Started: %s\n\n', datestr(now));

% Selection power zeta (Table S5)
fprintf('=== [1/2] SELECTION SWEEP ===\n');
s_values = [0.5, 1.0, 2.0, 4.0, 8.0];

for i = 1:length(s_values)
    s = s_values(i);
    fprintf('\n--- s=%.1f ---\n', s);
    config = publication_config('continuous_output', 'selectionPower', s);
    config.numReplicates = nReps;
    tic;
    results = evolve_population(config);
    elapsed = toc;
    coopValues = cellfun(@(r) r.historyCooperation(end), results);
    successRate = mean(coopValues >= 0.95) * 100;
    meanCoop = mean(coopValues) * 100;
    fprintf('  Success: %.0f%%, Mean coop: %.1f%%, Time: %.0f sec\n', successRate, meanCoop, elapsed);
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('sweep_selection_s%.1f_%dreps_%s.mat', s, nReps, timestamp);
    save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'successRate', 'meanCoop', 'elapsed', '-v7.3');
    fprintf('  Saved: %s\n', filename);
end

% Noise in gene expression of standard deviation sigma, added in evolve_population's Euler step
% (Table S5)
fprintf('\n=== [2/2] NOISE SWEEP ===\n');
sigma_values = [0.0, 0.01, 0.05, 0.10];

for i = 1:length(sigma_values)
    sigma = sigma_values(i);
    fprintf('\n--- noise sigma=%.2f ---\n', sigma);
    config = publication_config('continuous_output', 'noiseSigma', sigma);
    config.numReplicates = nReps;
    tic;
    results = evolve_population(config);
    elapsed = toc;
    coopValues = cellfun(@(r) r.historyCooperation(end), results);
    successRate = mean(coopValues >= 0.95) * 100;
    meanCoop = mean(coopValues) * 100;
    fprintf('  Success: %.0f%%, Mean coop: %.1f%%, Time: %.0f sec\n', successRate, meanCoop, elapsed);
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('sweep_noise_%.2f_%dreps_%s.mat', sigma, nReps, timestamp);
    save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'successRate', 'meanCoop', 'elapsed', '-v7.3');
    fprintf('  Saved: %s\n', filename);
end

fprintf('\n=== ALL SENSITIVITY SWEEPS COMPLETE: %s ===\n', datestr(now));
end
