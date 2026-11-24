function gen1_share_of_members()
% GEN1_SHARE_OF_MEMBERS  How many random starting wirings are members? The generation-1 share of
% members in each signaling architecture, from the membership test (600 steps, the simulate_pair
% route) on 10,000 random wirings.
%
% Produces : the generation-1 share of members (Table S3; SI S2.3); the stored member wirings are
%            read by scope_checks_gen1 and silent_and_moving_members
% Inputs   : publication/<architecture>_100reps.mat   (run_signaling_architectures; success only)
% Outputs  : publication/css_predicts_sensing.mat   (per architecture, named in `mode_names`: the
%              member count `n_css_stable`; the member wirings `css_stack_per_mode`, one row each;
%              `descriptor_table`, the share of members; success (%, final cooperation at least
%              0.95); `r_sensing`, the share's correlation with success over the four
%              architectures); a printed table of the same
% Run      : gen1_share_of_members   (about 10 min)
% Seed     : rng(2026) (the default generator) before each architecture's 10,000 random wirings,
%            weights N(0, 0.1^2), so all four test the same ones; membership_test seeds each wiring

    outDir = fullfile(grn_data_file('publication'));
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    N_per_mode = 10000;
    sigma_init = 0.1;
    sigma_mut  = 0.20;   % standard deviation of a weight mutation, variance 0.04 (Methods,
                         % "Game and evolution"); sets the test grid, -3 to +3 sigma_mut

    modes = {
        struct('name','graded-separate',     'sensingMode','separate',     'hillCoefficient',Inf, 'discreteOutput',false, 'file','graded_separate_100reps.mat'), ...
        struct('name','graded-contrast',     'sensingMode','contrast', 'hillCoefficient',Inf, 'discreteOutput',false, 'file','graded_contrast_100reps.mat'), ...
        struct('name','saturating-separate', 'sensingMode','separate',     'hillCoefficient',1,   'discreteOutput',false, 'file','saturating_separate_100reps.mat'), ...
        struct('name','saturating-contrast', 'sensingMode','contrast', 'hillCoefficient',1,   'discreteOutput',false, 'file','saturating_contrast_100reps.mat')
    };

    n_modes = length(modes);
    % success (%) of the four architectures (Table S2): replicates with final cooperation >= 0.95
    success_empirical = nan(n_modes, 1);
    for m = 1:n_modes
        d = load(grn_data_file(['publication/' modes{m}.file]), 'results');
        success_empirical(m) = 100 * sum(cellfun(@(r) r.historyCooperation(end) >= 0.95, d.results)) / numel(d.results);
    end
    clear d

    % The share of members per architecture
    descriptor_names = {'occupancy'};
    descriptor_table = nan(n_modes, 1);
    n_css_stable = zeros(n_modes, 1);
    css_stack_per_mode = cell(n_modes, 1);

    for m = 1:n_modes
        mode = modes{m};
        fprintf('\n=== Mode %d/%d: %s (N=%d) ===\n', m, n_modes, mode.name, N_per_mode);

        params = struct( ...
            'sensingMode',       mode.sensingMode, ...
            'hillCoefficient',   mode.hillCoefficient, ...
            'discreteOutput',    mode.discreteOutput, ...
            'michaelisConstant', 1.0, ...
            'degradationRate',   0.1, ...
            'timeStep',          1.0, ...
            'developmentSteps',  600, ...
            'sigma_mut',         sigma_mut, ...
            'backend',           'simulate_pair');

        rng(2026);  % reset for each architecture, so all four test the same wirings
        css_W = zeros(0, 16);
        n_css = 0;

        t0 = tic;
        for i = 1:N_per_mode
            W = randn(4) * sigma_init;
            r = membership_test(W, params);
            if r.is_css
                n_css = n_css + 1;
                css_W(end+1, :) = W(:)';                      %#ok<AGROW>
            end
            if mod(i, 500) == 0
                fprintf('  %d/%d (%.1f%%, members=%d, elapsed %.1fs)\n', ...
                    i, N_per_mode, 100*i/N_per_mode, n_css, toc(t0));
            end
        end

        n_css_stable(m) = n_css;
        css_stack_per_mode{m} = css_W;

        descriptor_table(m, 1) = n_css / N_per_mode;                                 % share of members

        fprintf('  Done in %.1fs. share=%.4f\n', toc(t0), descriptor_table(m, 1));
    end

    % Correlation r of the share of members with success over the four architectures
    r_sensing = corr(descriptor_table(:, 1), success_empirical, 'Type', 'Pearson');

    % Save
    outFile = fullfile(outDir, 'css_predicts_sensing.mat');
    mode_names = cellfun(@(m) m.name, modes, 'UniformOutput', false);
    save(outFile, 'mode_names', 'success_empirical', 'descriptor_table', ...
                  'descriptor_names', 'r_sensing', 'n_css_stable', ...
                  'css_stack_per_mode', 'sigma_init', 'sigma_mut', ...
                  'N_per_mode');

    % Print summary
    fprintf('\n=================================================================\n');
    fprintf('Summary: generation-1 share of members against success\n');
    fprintf('=================================================================\n');
    fprintf('%-22s | %-7s\n', 'architecture', 'members');
    for m = 1:n_modes
        fprintf('%-22s | %.4f\n', modes{m}.name, descriptor_table(m, 1));
    end
    fprintf('\nSuccess rates: %s\n', mat2str(success_empirical(:)'));
    fprintf('\nr (share of members vs success): %+.3f\n', r_sensing);
    fprintf('\nSaved to: %s\n', outFile);
end
