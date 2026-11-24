function run_development_steps(stepsList, nReps)
% RUN_DEVELOPMENT_STEPS  How does success depend on the limit on development steps? Evolves 100
% saturating-contrast populations at each limit (30, 50, 100, 200, 300, 500 steps).
%
% Produces : SI "Extended Materials and Methods", "Developmental integration and convergence"
%            (read by development_steps_summary)
% Inputs   : publication/saturating_contrast_100reps.mat   (its config; run_signaling_architectures)
% Outputs  : publication/devsteps_standard_steps<S>_<nReps>reps.mat, S the step limit
% Run      : run_development_steps   (all six limits, 300 first, 100 replicates; several hours),
%            then development_steps_summary
% Seed     : replicate r runs under rng(12345 + r, 'twister') at every limit (paired by replicate)
%
% Only developmentSteps (and numReplicates) differ from the stored config, so the 300-step run
% repeats the stored run; controlCheck stops the script if their success counts differ.

    if nargin < 1 || isempty(stepsList), stepsList = [300 30 50 100 200 500]; end
    if nargin < 2 || isempty(nReps), nReps = 100; end

    outDir  = fullfile(grn_data_file('publication'));
    srcFile = fullfile(outDir, 'saturating_contrast_100reps.mat');

    src = load(srcFile, 'config');
    base = src.config;
    base.numReplicates = nReps;

    tAll = tic;
    for S = stepsList
        config = base;
        config.developmentSteps = S;

        t0 = tic;
        results = evolve_population(config);
        elapsed = toc(t0);

        coopValues  = cellfun(@(r) r.historyCooperation(end), results);
        successRate = mean(coopValues >= 0.95) * 100;
        meanCoop    = mean(coopValues) * 100;
        configSource = srcFile; %#ok<NASGU>

        f = fullfile(outDir, sprintf('devsteps_standard_steps%d_%dreps.mat', S, nReps));
        save(f, 'results', 'config', 'coopValues', 'successRate', 'meanCoop', ...
             'elapsed', 'configSource', '-v7.3');
        fprintf('CELL steps=%d  success=%d/%d  meancoop=%.1f%%  (%.0f s)  -> %s\n', ...
            S, sum(coopValues >= 0.95), nReps, meanCoop, elapsed, f);

        if S == 300 && nReps == 100
            controlCheck(results, srcFile);
        end
    end
    fprintf('TOTAL %.1f min\n', toc(tAll)/60);
end


function controlCheck(results, srcFile)
% The 300-step cell against the stored saturating-contrast run.
    pub = load(srcFile, 'results');
    n = numel(results);
    same = false(n, 1);
    for r = 1:n
        same(r) = isequal(results{r}.historyCooperation, pub.results{r}.historyCooperation);
    end
    vNew = cellfun(@(x) x.historyCooperation(end), results);
    vPub = cellfun(@(x) x.historyCooperation(end), pub.results);
    kNew = sum(vNew >= 0.95); kPub = sum(vPub >= 0.95);

    fprintf('\nCONTROL (300 steps vs saturating_contrast_100reps.mat)\n');
    fprintf('  success      %d/%d  vs published %d/%d\n', kNew, n, kPub, n);
    fprintf('  mean coop    %.4f%% vs published %.4f%%\n', 100*mean(vNew), 100*mean(vPub));
    fprintf('  historyCooperation bit-identical in %d/%d replicates\n', sum(same), n);
    if ~all(same)
        fprintf('  differing replicates: %s\n', mat2str(find(~same)'));
        fprintf('  max |final coop diff| = %.3g\n', max(abs(double(vNew) - double(vPub))));
    end
    if kNew ~= kPub
        error('run_development_steps:controlFailed', ...
              'Control cell gives %d/%d, published %d/%d. Stop before other cells.', ...
              kNew, n, kPub, n);
    end
end
