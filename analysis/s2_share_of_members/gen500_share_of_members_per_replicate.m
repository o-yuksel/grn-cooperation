function gen500_share_of_members_per_replicate()
% GEN500_SHARE_OF_MEMBERS_PER_REPLICATE  How much does the share of members vary between replicate
% populations? Per signaling architecture, each replicate's generation-500 share of members on 50
% sampled wirings, and ||x*|| (the size of the settled expression state) of every member found.
%
% Produces : the input that share_interval_effective_n reads for the width of the interval on the
%            share of members (SI S2.3)
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures; the
%              newest matching file)
% Outputs  : publication/share_of_members_cloud_sensing.mat   (`N_SUB`, wirings sampled per replicate; `conds`,
%              per architecture: name; success (%, Table S2, written in); member_x, ||x*|| of each
%              member; member_reach, the share of members of each member's replicate; rep_reach,
%              each replicate's share of members; reach_mean, reach_std; xnorm_median); prints per
%              architecture the mean share, the members found and their median ||x*||
% Run      : gen500_share_of_members_per_replicate   (not timed; it tests 20,000 wirings, where
%              gen500_share_of_members tests 12,000 in about 5 min)
% Seed     : rng(42) (the default generator) once before the four architectures; randperm then
%            draws 50 wirings per replicate, architectures in order, in one stream;
%            membership_test seeds each wiring itself
%
% The membership test uses each run's settings, 600 steps and the simulate_pair route for all
% four. Stored names say "reach" for the share of members.

    pub = fullfile(grn_data_file('publication'));

    N_SUB = 50;
    rng(42);

    modes = {
        struct('name', 'graded-separate',     'succ', 7,  'pat', 'graded_separate_100reps.mat')
        struct('name', 'graded-contrast',     'succ', 34, 'pat', 'graded_contrast_100reps.mat')
        struct('name', 'saturating-separate', 'succ', 58, 'pat', 'saturating_separate_100reps.mat')
        struct('name', 'saturating-contrast', 'succ', 76, 'pat', 'saturating_contrast_100reps.mat')
    };
    modes = vertcat(modes{:});

    conds = struct('name', {}, 'success', {}, 'member_x', {}, ...
                   'member_reach', {}, 'rep_reach', {}, 'reach_mean', {}, ...
                   'reach_std', {}, 'xnorm_median', {});

    for m = 1:numel(modes)
        s = modes(m);
        files = dir(fullfile(pub, s.pat));
        assert(~isempty(files), 'No file matching %s', s.pat);
        [~, mx] = max([files.datenum]);
        d = load(fullfile(pub, files(mx).name));
        cfg = d.config;
        % sigma_mut is membership_test's default, 0.2, the standard deviation of a weight mutation
        params = struct('sensingMode', cfg.sensingMode, ...
            'hillCoefficient', cfg.hillCoefficient, ...
            'discreteOutput', cfg.discreteOutput, ...
            'michaelisConstant', cfg.michaelisConstant, ...
            'degradationRate', cfg.degradationRate, ...
            'timeStep', cfg.timeStep, 'developmentSteps', 600, ...
            'backend', 'simulate_pair');

        nrep = numel(d.results);
        member_x     = [];
        member_reach = [];
        rep_reach    = nan(nrep, 1);
        t0 = tic;
        for r = 1:nrep
            B = d.results{r}.finalNetworkB;
            M = d.results{r}.finalNetworkM;
            na = size(B, 3);
            pick = randperm(na, min(N_SUB, na));
            xnorms_r = zeros(0, 1);
            for a = pick
                W = double(M(:, :, a)) .* double(B(:, :, a));
                res = membership_test(reshape(W, 4, 4), params);
                if res.is_css
                    xnorms_r(end+1, 1) = res.x_star_norm;  %#ok<AGROW>
                end
            end
            % the replicate's share of members: members among its sampled wirings
            rr = numel(xnorms_r) / numel(pick);
            rep_reach(r) = rr;
            member_x     = [member_x;     xnorms_r];             %#ok<AGROW>
            member_reach = [member_reach; repmat(rr, numel(xnorms_r), 1)]; %#ok<AGROW>
        end

        conds(m).name         = s.name;
        conds(m).success      = s.succ;
        conds(m).member_x     = member_x;
        conds(m).member_reach = member_reach;
        conds(m).rep_reach    = rep_reach;           % per-replicate share of members (n = nrep)
        conds(m).reach_mean   = mean(rep_reach, 'omitnan');
        conds(m).reach_std    = std(rep_reach, 'omitnan');
        conds(m).xnorm_median = median(member_x);
        fprintf('%-22s share=%.3f  %d members  medianX=%.2f (%.1fs)\n', ...
            s.name, conds(m).reach_mean, numel(member_x), ...
            conds(m).xnorm_median, toc(t0));
    end

    outFile = fullfile(pub, 'share_of_members_cloud_sensing.mat');
    save(outFile, 'conds', 'N_SUB');
    fprintf('\nSaved %s\n', outFile);
end
