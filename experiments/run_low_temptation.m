function run_low_temptation(nReps)
% RUN_LOW_TEMPTATION  Evolves 100 saturating-contrast populations at two temptations, T = 2.5 and
% T = 3.0, with R = 3, P = 1 and S = 0. At T = 2.5 each player does best to match its partner (a
% coordination game, the stag hunt); at T = 3.0 defecting on a cooperator earns the same as
% cooperating with it.
%
% Produces : the T = 2.5 and T = 3.0 rows of Table S5 (read by tableS5_parameter_sweeps)
% Inputs   : none (the standard setting of publication_config with T changed)
% Outputs  : publication/sweep_T_<T>_<nReps>reps_<time>.mat   (results, config, elapsed)
% Run      : run_low_temptation   (nReps = 100, the paper's runs; run_low_temptation(n) runs n)
% Seed     : replicate r runs under rng(12345 + r, 'twister'), at both T

    if nargin < 1, nReps = 100; end

    tValues = [2.5, 3.0];
    nSettings = length(tValues);

    resultsDir = fullfile(grn_data_file('publication'));
    if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('  T-BOUNDARY SWEEP: Where does the capacity constraint bind?\n');
    fprintf('  T values: %s\n', mat2str(tValues));
    fprintf('  Config: publication_config (the paper defaults, Hill n=1, contrast)\n');
    fprintf('  %d reps per setting\n', nReps);
    fprintf('%s\n\n', repmat('=', 1, 70));

    % Game type of each T against the reward R
    R = 3; P = 1; S = 0;
    for i = 1:nSettings
        T = tValues(i);
        [gameType, ~] = game_name(T, R, P, S);
        fprintf('  T=%.1f: %s\n', T, gameType);
    end
    fprintf('\n');

    summarySuccess = zeros(nSettings, 1);
    summaryMeanCoop = zeros(nSettings, 1);
    summaryElapsed = zeros(nSettings, 1);
    totalTimer = tic;

    for c = 1:nSettings
        T = tValues(c);
        fprintf('\n%s\n', repmat('=', 1, 70));
        fprintf('  [%d/%d] T = %.1f\n', c, nSettings, T);
        fprintf('%s\n', repmat('=', 1, 70));

        config = publication_config('continuous_output', 'payoffDC', T);
        config.numReplicates = nReps;

        settingTimer = tic;
        results = evolve_population(config);
        elapsed = toc(settingTimer);
        summaryElapsed(c) = elapsed;

        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('sweep_T_%.1f_%dreps_%s.mat', T, nReps, timestamp);
        save(fullfile(resultsDir, filename), 'results', 'config', 'elapsed', '-v7.3');

        % Success: cooperation at least 0.95 in the last generation
        coop = zeros(nReps, 1);
        for i = 1:nReps
            coop(i) = results{i}.historyCooperation(end);
        end
        summarySuccess(c) = mean(coop >= 0.95) * 100;
        summaryMeanCoop(c) = mean(coop) * 100;

        % Clopper-Pearson 95% interval of the success share
        nSucc = sum(coop >= 0.95);
        [~, ci] = binofit(nSucc, nReps, 0.05);

        fprintf('\n  T=%.1f | Success: %.0f%% [%.1f, %.1f] | Mean coop: %.1f%% | Time: %.1f min\n', ...
            T, summarySuccess(c), ci(1)*100, ci(2)*100, summaryMeanCoop(c), elapsed/60);
        fprintf('  Saved: %s\n', filename);
    end

    totalElapsed = toc(totalTimer);

    fprintf('\n\n%s\n', repmat('=', 1, 70));
    fprintf('  T-BOUNDARY RESULTS\n');
    fprintf('%s\n', repmat('=', 1, 70));
    fprintf('  %-6s  %-20s  %-10s  %-10s\n', 'T', 'Game type', 'Success', 'Mean Coop');
    fprintf('  %s\n', repmat('-', 1, 60));

    for c = 1:nSettings
        T = tValues(c);
        [~, gt] = game_name(T, R, P, S);
        fprintf('  T=%.1f   %-20s  %5.0f%%     %5.1f%%\n', ...
            T, gt, summarySuccess(c), summaryMeanCoop(c));
    end

    fprintf('\n  Total time: %.1f minutes\n', totalElapsed/60);
    fprintf('%s\n', repmat('=', 1, 70));
end


function [long, short] = game_name(T, R, P, S)
% The symmetric 2 x 2 game of the payoff order (T temptation, R reward, P punishment, S sucker's).
    if T > R && R > P && P > S
        long = 'Prisoner''s Dilemma (T>R>P>S)'; short = 'PD';
    elseif T == R && R > P && P > S
        long = 'boundary of Prisoner''s Dilemma and Stag Hunt (T=R>P>S)'; short = 'PD/SH boundary';
    elseif R > T && T > P && P > S
        long = 'Stag Hunt (R>T>P>S)'; short = 'Stag Hunt';
    else
        long = sprintf('other order (T=%g, R=%g, P=%g, S=%g)', T, R, P, S); short = 'other';
    end
end
