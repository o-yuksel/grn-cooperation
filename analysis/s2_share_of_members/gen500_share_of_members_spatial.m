function gen500_share_of_members_spatial()
% GEN500_SHARE_OF_MEMBERS_SPATIAL  What share of the final populations are members in the five
% spatial settings (N = 50, well mixed and rings with neighborhood size k = 1, 2, 5, 10)? The
% membership test on every generation-500 wiring, and the share's correlation with success.
%
% Produces : the five spatial shares of members in Fig. 3 and Table S2
% Inputs   : publication/sweep_spatial_N50_<setting>_100reps_<time>.mat   (written by
%              run_spatial_n50; the newest file of each setting)
% Outputs  : publication/css_predicts_spatial.mat   (`reach_occupancy`, the share of members;
%              `n_in_Vm`, members; `n_evolved_total`, wirings tested; success (%, final cooperation
%              at least 0.95); `r_spatial`, the share's correlation with success; the runs' config)
% Run      : gen500_share_of_members_spatial   (about 8 min)
% Seed     : none (every wiring is tested in order; membership_test seeds each wiring itself)
%
% The membership test uses the runs' own settings (all five are saturating contrast; the run stops
% if their sensing differs), 600 steps, the simulate_pair route.

    pubDir = fullfile(grn_data_file('publication'));

    spatial_files = {
        struct('cond','wellmixed', 'pat','sweep_spatial_N50_wellmixed_100reps_*.mat'), ...
        struct('cond','ring_k1',   'pat','sweep_spatial_N50_ring_k1_100reps_*.mat'), ...
        struct('cond','ring_k2',   'pat','sweep_spatial_N50_ring_k2_100reps_*.mat'), ...
        struct('cond','ring_k5',   'pat','sweep_spatial_N50_ring_k5_100reps_*.mat'), ...
        struct('cond','ring_k10',  'pat','sweep_spatial_N50_ring_k10_100reps_*.mat')
    };

    n_settings = length(spatial_files);

    success_empirical = nan(n_settings, 1);                        % success (%) per setting (Table S2)
    reach_occupancy = nan(n_settings, 1);
    n_evolved_total = zeros(n_settings, 1);
    n_in_Vm         = zeros(n_settings, 1);
    spatial_sensing_config = [];

    for c = 1:n_settings
        sf = spatial_files{c};
        files = dir(fullfile(pubDir, sf.pat));
        assert(~isempty(files), sprintf('No files matching %s', sf.pat));
        [~, idx] = max([files.datenum]);
        fpath = fullfile(pubDir, files(idx).name);
        fprintf('\n=== Setting %d/%d: %s (loading %s) ===\n', c, n_settings, sf.cond, files(idx).name);

        data = load(fpath);
        success_empirical(c) = 100 * sum(cellfun(@(r) r.historyCooperation(end) >= 0.95, data.results)) / numel(data.results);

        if c == 1
            spatial_sensing_config = data.config;
            fprintf('Spatial runs use: sensingMode=%s, hillCoefficient=%g, discreteOutput=%d\n', ...
                data.config.sensingMode, data.config.hillCoefficient, data.config.discreteOutput);
        else
            assert(strcmp(stored_value('sensingMode', data.config.sensingMode), stored_value('sensingMode', spatial_sensing_config.sensingMode)), ...
                'Sensing mode mismatch across spatial files');
        end

        params = struct( ...
            'sensingMode',       data.config.sensingMode, ...
            'hillCoefficient',   data.config.hillCoefficient, ...
            'discreteOutput',    data.config.discreteOutput, ...
            'michaelisConstant', data.config.michaelisConstant, ...
            'degradationRate',   data.config.degradationRate, ...
            'timeStep',          data.config.timeStep, ...
            'developmentSteps',  600, ...
            ... % sigma_mut is membership_test's default, 0.2, the standard deviation of
            ... % a weight mutation; the simulate_pair route, as for the signaling
            ... % architectures.
            'backend',           'simulate_pair');

        % Extract every final-generation wiring across all replicates
        Wstack = zeros(0, 16);
        for r = 1:numel(data.results)
            B_final = data.results{r}.finalNetworkB;  % 4x4xN
            M_final = data.results{r}.finalNetworkM;  % 4x4xN
            n_agents = size(B_final, 3);
            for a = 1:n_agents
                W_a = double(M_final(:,:,a)) .* double(B_final(:,:,a));
                Wstack(end+1, :) = W_a(:)';                                 %#ok<AGROW>
            end
        end
        n_evolved_total(c) = size(Wstack, 1);
        fprintf('  Loaded %d evolved wirings across %d replicates.\n', n_evolved_total(c), numel(data.results));

        member_W = zeros(0, 16);
        t0 = tic;
        for k = 1:n_evolved_total(c)
            W = reshape(Wstack(k, :), 4, 4);
            res = membership_test(W, params);
            if res.is_css
                member_W(end+1, :) = Wstack(k, :);                              %#ok<AGROW>
            end
            if mod(k, 500) == 0
                fprintf('    %d/%d (%.1f%%, members=%d, elapsed %.1fs)\n', ...
                    k, n_evolved_total(c), 100*k/n_evolved_total(c), ...
                    size(member_W, 1), toc(t0));
            end
        end

        n_in_Vm(c) = size(member_W, 1);
        reach_occupancy(c) = n_in_Vm(c) / n_evolved_total(c);
        fprintf('  Done in %.1fs. share=%.4f\n', toc(t0), reach_occupancy(c));
    end

    % Correlation r of the share of members with success over the five settings
    r_spatial = corr(reach_occupancy, success_empirical, 'Type', 'Pearson');

    condition_names = cellfun(@(s) s.cond, spatial_files, 'UniformOutput', false);
    outFile = fullfile(pubDir, 'css_predicts_spatial.mat');
    save(outFile, 'condition_names', 'success_empirical', ...
                  'reach_occupancy', 'r_spatial', ...
                  'n_evolved_total', 'n_in_Vm', 'spatial_sensing_config');

    fprintf('\n=================================================================\n');
    fprintf('Summary: generation-500 share of members against success, spatial settings\n');
    fprintf('=================================================================\n');
    fprintf('%-12s | %-13s | %-9s\n', 'setting', 'success(%)', 'members');
    for c = 1:n_settings
        fprintf('%-12s | %5.0f         | %.4f\n', ...
            condition_names{c}, success_empirical(c), reach_occupancy(c));
    end
    fprintf('\nr with success:\n');
    fprintf('  share of members vs success: %+.3f\n', r_spatial(1));
    fprintf('\nSaved to: %s\n', outFile);
end
