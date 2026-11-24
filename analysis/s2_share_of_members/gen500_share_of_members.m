function gen500_share_of_members()
% GEN500_SHARE_OF_MEMBERS  What share of each signaling architecture's final populations are
% members? The membership test on 30 sampled wirings per replicate at generation 500, and the
% correlation of that share with success over the four architectures.
%
% Produces : the four architectures' shares of members in Fig. 3 and Table S2
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures; the
%              newest file matching the name)
% Outputs  : publication/css_predicts_sensing_gen500.mat   (`reach_occupancy`, the share of members;
%              `n_in_Vm`, the member count; `n_subsample`, the wirings tested; success (%, final
%              cooperation at least 0.95); `r_sensing_reach`, the correlation), read by fig3_data;
%              a printed table
% Run      : gen500_share_of_members   (about 5 min)
% Seed     : rng(42) (the default generator) once before the four architectures; randperm then
%            draws 30 wirings per replicate, architectures in order, in one stream;
%            membership_test seeds each wiring itself
%
% Each architecture's wirings get the membership test with that run's settings, 600 steps, the
% simulate_pair route: 3,000 wirings per architecture from 100 replicates.

    pubDir = fullfile(grn_data_file('publication'));

    modes = {
        struct('name','graded-separate',     'pat','graded_separate_100reps.mat'), ...
        struct('name','graded-contrast',     'pat','graded_contrast_100reps.mat'), ...
        struct('name','saturating-separate', 'pat','saturating_separate_100reps.mat'), ...
        struct('name','saturating-contrast', 'pat','saturating_contrast_100reps.mat')
    };
    n_mode  = numel(modes);
    success_empirical = nan(n_mode, 1);    % success (%) of the four architectures (Table S2)
    n_sub   = 30;                       % wirings sampled per replicate (Methods)
    rng(42);                            % subsample seed, once for all four architectures
    reach_occupancy = nan(n_mode, 1);
    n_subsample     = zeros(n_mode, 1);
    n_in_Vm         = zeros(n_mode, 1);
    condition_names = cell(n_mode, 1);

    for m = 1:n_mode
        sf = modes{m};
        files = dir(fullfile(pubDir, sf.pat));
        assert(~isempty(files), sprintf('No files matching %s', sf.pat));
        [~, idx] = max([files.datenum]);
        fpath = fullfile(pubDir, files(idx).name);
        fprintf('\n=== Mode %d/%d: %s (%s) ===\n', m, n_mode, sf.name, files(idx).name);

        data = load(fpath);
        cfg  = data.config;
        success_empirical(m) = 100 * sum(cellfun(@(r) r.historyCooperation(end) >= 0.95, data.results)) / numel(data.results);
        fprintf('  config: sensingMode=%s, hillCoef=%g, discreteOutput=%d\n', ...
                cfg.sensingMode, cfg.hillCoefficient, cfg.discreteOutput);

        params = struct( ...
            'sensingMode',       cfg.sensingMode, ...
            'hillCoefficient',   cfg.hillCoefficient, ...
            'discreteOutput',    cfg.discreteOutput, ...
            'michaelisConstant', cfg.michaelisConstant, ...
            'degradationRate',   cfg.degradationRate, ...
            'timeStep',          cfg.timeStep, ...
            'developmentSteps',  600, ...
            ... % sigma_mut is membership_test's default, 0.2, the standard deviation of
            ... % a weight mutation. The route is simulate_pair for all four: 'auto'
            ... % would send graded separate alone to the linear route, which
            ... % requires convergence and counts a mutant pair that does not
            ... % converge as failing the retaliation condition, a stricter rule than
            ... % the other three get (membership_test, the comment on its routes).
            'backend',           'simulate_pair');

        % n_sub final-generation wirings per replicate, W = M .* B, one row per wiring
        Wstack = zeros(0, 16);
        for r = 1:numel(data.results)
            B = data.results{r}.finalNetworkB;
            M = data.results{r}.finalNetworkM;
            n_agents = size(B, 3);
            pick = randperm(n_agents, min(n_sub, n_agents));
            for a = pick
                W = double(M(:,:,a)) .* double(B(:,:,a));
                Wstack(end+1, :) = W(:)';                                  %#ok<AGROW>
            end
        end
        n_subsample(m) = size(Wstack, 1);
        fprintf('  Subsampled %d wirings (%d per rep x %d reps)\n', ...
                n_subsample(m), n_sub, numel(data.results));

        n_pass = 0;
        t0 = tic;
        for k = 1:n_subsample(m)
            W = reshape(Wstack(k, :), 4, 4);
            res = membership_test(W, params);
            if res.is_css, n_pass = n_pass + 1; end
            if mod(k, 500) == 0
                fprintf('    %d/%d (%.1f%%, members=%d, elapsed %.1fs)\n', ...
                        k, n_subsample(m), 100*k/n_subsample(m), n_pass, toc(t0));
            end
        end
        n_in_Vm(m)         = n_pass;
        reach_occupancy(m) = n_pass / n_subsample(m);
        condition_names{m} = sf.name;
        fprintf('  Done in %.1fs. share of members = %.4f (%d/%d)\n', ...
                toc(t0), reach_occupancy(m), n_pass, n_subsample(m));
    end

    r_sensing_reach = corr(reach_occupancy, success_empirical, 'Type', 'Pearson');

    outFile = fullfile(pubDir, 'css_predicts_sensing_gen500.mat');
    save(outFile, 'condition_names', 'success_empirical', ...
                  'reach_occupancy', 'n_subsample', 'n_in_Vm', ...
                  'r_sensing_reach');

    fprintf('\n=================================================================\n');
    fprintf('Generation-500 share of members against success\n');
    fprintf('=================================================================\n');
    fprintf('%-22s | %-10s | %-12s | %s\n', 'architecture', 'success(%)', 'members', 'members/sampled');
    for m = 1:n_mode
        fprintf('%-22s | %5.0f      | %.4f      | %d/%d\n', ...
                condition_names{m}, success_empirical(m), ...
                reach_occupancy(m), n_in_Vm(m), n_subsample(m));
    end
    fprintf('\nr (share of members vs success): %+.3f\n', r_sensing_reach);
    fprintf('\nSaved to: %s\n', outFile);
end
