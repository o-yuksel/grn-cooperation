function run_mutation_size_sweep(nReps)
% RUN_MUTATION_SIZE_SWEEP  Evolves 100 populations at each size sigma_mut of a weight mutation (its
% standard deviation, 0.01 to 1.0), at saturating contrast and at graded separate.
%
% Produces : the mutation-size rows of Table S5 (Part 1, saturating contrast; read by
%            tableS5_parameter_sweeps). Part 2 (graded separate) feeds only a range that
%            tableS5_parameter_sweeps prints and the SI does not state.
% Inputs   : none (the standard setting of publication_config, mutationStdWeight = sigma_mut)
% Outputs  : publication/sweep_mutation_size_saturating_contrast_<sigma_mut>_<nReps>reps_<time>.mat (Part 1) and
%            publication/sweep_mutation_size_graded_separate_<sigma_mut>_<nReps>reps_<time>.mat (Part 2),
%            sigma_mut written 0.010 to 1.000; a value whose file exists is skipped
% Run      : run_mutation_size_sweep   (nReps = 100; several hours on the parallel pool)
% Seed     : replicate r runs under rng(12345 + r, 'twister') at every value

    if nargin < 1, nReps = 100; end

    resultsDir = fullfile(grn_data_file('publication'));
    if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

    totalTimer = tic;

    % Part 1: weight-mutation size sigma_mut at saturating contrast (Table S5)
    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('  PART 1: sigma_mut sweep Sat+Con (%d reps)\n', nReps);
    fprintf('%s\n\n', repmat('=', 1, 70));

    mutValues = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
    nMut = length(mutValues);

    for c = 1:nMut
        mutStd = mutValues(c);
        existingFiles = dir(fullfile(resultsDir, ...
            sprintf('sweep_mutation_size_saturating_contrast_%.3f_%dreps_*.mat', mutStd, nReps)));
        if ~isempty(existingFiles)
            fprintf('  [%d/%d] sigma_mut=%.3f: exists, skipping\n', c, nMut, mutStd);
            continue;
        end
        fprintf('  [%d/%d] sigma_mut=%.3f (Sat+Con, %d reps)\n', c, nMut, mutStd, nReps);
        config = publication_config('continuous_output');
        config.mutationStdWeight = mutStd;
        config.numReplicates = nReps;
        condTimer = tic;
        results = evolve_population(config);
        elapsed = toc(condTimer);
        coopValues = zeros(nReps, 1);
        for i = 1:nReps
            coopValues(i) = results{i}.historyCooperation(end);
        end
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('sweep_mutation_size_saturating_contrast_%.3f_%dreps_%s.mat', mutStd, nReps, timestamp);
        save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'elapsed', '-v7.3');
        fprintf('    Success: %.0f%% | Mean coop: %.1f%% | Time: %.1f min\n', ...
            mean(coopValues >= 0.95)*100, mean(coopValues)*100, elapsed/60);
    end

    % Part 2: weight-mutation size sigma_mut at graded separate (not a row of Table S5; its range
    % is printed by tableS5_parameter_sweeps)
    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('  PART 2: sigma_mut sweep Grad+Sep (%d reps)\n', nReps);
    fprintf('%s\n\n', repmat('=', 1, 70));

    for c = 1:nMut
        mutStd = mutValues(c);
        existingFiles = dir(fullfile(resultsDir, ...
            sprintf('sweep_mutation_size_graded_separate_%.3f_%dreps_*.mat', mutStd, nReps)));
        if ~isempty(existingFiles)
            fprintf('  [%d/%d] sigma_mut=%.3f: exists, skipping\n', c, nMut, mutStd);
            continue;
        end
        fprintf('  [%d/%d] sigma_mut=%.3f (Grad+Sep, %d reps)\n', c, nMut, mutStd, nReps);
        config = publication_config('continuous_output');
        config.hillCoefficient = Inf;
        config.sensingMode = 'separate';
        config.mutationStdWeight = mutStd;
        config.numReplicates = nReps;
        condTimer = tic;
        results = evolve_population(config);
        elapsed = toc(condTimer);
        coopValues = zeros(nReps, 1);
        for i = 1:nReps
            coopValues(i) = results{i}.historyCooperation(end);
        end
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('sweep_mutation_size_graded_separate_%.3f_%dreps_%s.mat', mutStd, nReps, timestamp);
        save(fullfile(resultsDir, filename), 'results', 'config', 'coopValues', 'elapsed', '-v7.3');
        fprintf('    Success: %.0f%% | Mean coop: %.1f%% | Time: %.1f min\n', ...
            mean(coopValues >= 0.95)*100, mean(coopValues)*100, elapsed/60);
    end

    totalElapsed = toc(totalTimer);
    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('  ALL PARTS COMPLETE. Total: %.1f hours\n', totalElapsed/3600);
    fprintf('%s\n', repmat('=', 1, 70));
end
