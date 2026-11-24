function run_spatial_n50(nReps)
% RUN_SPATIAL_N50  Evolves saturating-contrast populations of N = 50, well mixed or on a ring with
% k = 1, 2, 5 or 10 neighbors on each side, 100 populations per setting.
%
% Produces : the runs behind SI "Extended Materials and Methods", "Spatial reproduction", Table S2
%            and Fig. 3 (network reciprocity); Methods, "Game and evolution"
% Inputs   : none (the standard setting of publication_config with N and the ring changed)
% Outputs  : publication/sweep_spatial_N50_<setting>_<nReps>reps_<time>.mat   (setting
%              wellmixed, ring_k1, ring_k2, ring_k5 or ring_k10; a setting whose file exists is
%              read back) and publication/sweep_spatial_N50_summary_<nReps>reps_<time>.mat
% Run      : run_spatial_n50   (nReps = 100, the paper's runs)
% Seed     : replicate r runs under rng(12345 + r, 'twister') in every setting
%
% On the ring (in evolve_population) the partner is one of the 2k ring neighbors, and each
% position is refilled by a parent drawn fitness-proportionally from them and its own occupant.

    if nargin < 1, nReps = 100; end

    resultsDir = fullfile(grn_data_file('publication'));
    if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end

    spatialSettings = {
        struct('name', 'wellmixed', 'spatialStructure', 'wellmixed', 'k', 0);
        struct('name', 'ring_k1',   'spatialStructure', 'ring',      'k', 1);
        struct('name', 'ring_k2',   'spatialStructure', 'ring',      'k', 2);
        struct('name', 'ring_k5',   'spatialStructure', 'ring',      'k', 5);
        struct('name', 'ring_k10',  'spatialStructure', 'ring',      'k', 10);
    };

    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('  N=50 SPATIAL SWEEP\n');
    fprintf('  %d settings x %d reps\n', numel(spatialSettings), nReps);
    fprintf('%s\n\n', repmat('=', 1, 70));

    summary = table();

    for c = 1:numel(spatialSettings)
        setting = spatialSettings{c};

        existing = dir(fullfile(resultsDir, ...
            sprintf('sweep_spatial_N50_%s_%dreps_*.mat', setting.name, nReps)));
        if ~isempty(existing)
            fprintf('  [%d/%d] %s: exists, skipping\n', c, numel(spatialSettings), setting.name);
            d = load(fullfile(resultsDir, existing(end).name), 'finalCoop');
            successRate = mean(d.finalCoop >= 0.95) * 100;
            meanCoop    = mean(d.finalCoop) * 100;
            row = table(string(setting.name), 50, successRate, meanCoop, 0, ...
                'VariableNames', {'Setting', 'N', 'SuccessRate', 'MeanCoop', 'ElapsedSec'});
            summary = [summary; row]; %#ok<AGROW>
            continue;
        end

        fprintf('  [%d/%d] %s (k=%d)\n', c, numel(spatialSettings), setting.name, setting.k);

        config = publication_config('continuous_output', ...
            'populationSize',   50, ...
            'spatialStructure', setting.spatialStructure, ...
            'neighborhoodSize', setting.k);
        config.numReplicates = nReps;
        config.conditionName = ['spatial_N50_' setting.name];

        condTimer = tic;
        results = evolve_population(config);
        elapsed = toc(condTimer);

        finalCoop = zeros(nReps, 1);
        for r = 1:nReps
            finalCoop(r) = results{r}.historyCooperation(end);
        end
        successRate = mean(finalCoop >= 0.95) * 100;
        meanCoop    = mean(finalCoop) * 100;

        fprintf('    Success: %.0f%%  MeanCoop: %.1f%%  Time: %.1f min\n', ...
            successRate, meanCoop, elapsed/60);

        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        fname = sprintf('sweep_spatial_N50_%s_%dreps_%s.mat', setting.name, nReps, timestamp);
        save(fullfile(resultsDir, fname), ...
            'results', 'config', 'finalCoop', 'successRate', 'meanCoop', 'elapsed', '-v7.3');

        row = table(string(setting.name), 50, successRate, meanCoop, elapsed, ...
            'VariableNames', {'Setting', 'N', 'SuccessRate', 'MeanCoop', 'ElapsedSec'});
        summary = [summary; row]; %#ok<AGROW>
    end

    fprintf('\n========== N=50 SPATIAL SWEEP SUMMARY ==========\n');
    disp(summary);

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save(fullfile(resultsDir, sprintf('sweep_spatial_N50_summary_%dreps_%s.mat', nReps, timestamp)), ...
        'summary');
    fprintf('Summary saved.\n');
end
