function out = silent_and_moving_members()
% SILENT_AND_MOVING_MEMBERS  How many members have a clone state that is silent or still moving,
% and do the shares of members change when those members are set aside? Generation 1 in the four
% architectures; generation 500 in the nine signaling and spatial settings.
%
% Produces : SI S2.3 (members whose clone state is silent or still moving cooperate in the game and
%            count as members)
% Inputs   : stored members, shares and success, for the checks and the "as stored" columns:
%              publication/css_predicts_sensing.mat (gen1_share_of_members),
%              css_predicts_sensing_gen500.mat (gen500_share_of_members), css_predicts_spatial.mat
%              (gen500_share_of_members_spatial), css_radius_sensitivity.mat (part of the data
%              archive)
%            final wirings: publication/<architecture>_100reps.mat (run_signaling_architectures),
%              sweep_spatial_N50_<setting>_100reps_<time>.mat (run_spatial_n50)
% Outputs  : publication/member_degeneracy.mat   (per setting the flags and counts; the shares,
%              correlations with success and spreads, as stored and with flagged members set
%              aside); console tables
% Run      : silent_and_moving_members   (parfor on the current pool; about 4 min on 3 workers)
% Seed     : rng(2026) (MATLAB's default generator) before the 10,000 generation-1 draws
%            randn(4) * 0.1, one set for all four architectures; rng(42) once before the four
%            signaling architectures; membership_test seeds each wiring itself
%
% Flags, on the clone state the membership test reads (600 steps from the fixed start): silent,
% max(x_C, x_D) below 1e-3 (and, separately, 1e-2); still moving, rerun for 20,000 steps without the
% stop on convergence, the state moves by more than 5% or the readout is not positive. The run stops
% unless it reproduces the stored members (55 / 219 / 109 / 301 at generation 1, and the stored
% membership of every generation-500 wiring).

    pub = fullfile(grn_data_file('publication'));
    tStart = tic;

    FLOORS = [1e-3, 1e-2];
    g1 = gen1_pass(pub, FLOORS);
    t_gen1 = toc(tStart);
    fprintf('\nGen-1 pass done in %.1f min.\n', t_gen1 / 60);
    g500 = gen500_pass(pub, FLOORS);
    runtime_s = toc(tStart);

    out = summarize(g1, g500, FLOORS);
    out.runtime_s = runtime_s; out.t_gen1 = t_gen1;
    outFile = fullfile(pub, 'member_degeneracy.mat');
    save(outFile, 'g1', 'g500', 'out', 'FLOORS', 'runtime_s', 't_gen1');
    fprintf('\nRuntime %.1f min (gen 1 %.1f). Saved %s\n', runtime_s / 60, t_gen1 / 60, outFile);
end


% =========================================================================
function g1 = gen1_pass(pub, FLOORS)
    D = load(fullfile(pub, 'css_predicts_sensing.mat'), 'css_stack_per_mode');
    names = {'graded separate', 'graded contrast', 'saturating separate', 'saturating contrast'};
    cfgs  = {{'separate', Inf}, {'contrast', Inf}, {'separate', 1}, {'contrast', 1}};
    stored_n = [55 219 109 301];
    nW = 10000;
    rng(2026);                       % same stack for every architecture
    Ws = zeros(nW, 16);
    for i = 1:nW, W = randn(4) * 0.1; Ws(i, :) = W(:)'; end
    for m = 1:4
        cfg = struct('sensingMode', cfgs{m}{1}, 'hillCoefficient', cfgs{m}{2}, ...
            'discreteOutput', false, 'michaelisConstant', 1.0, 'degradationRate', 0.1, ...
            'timeStep', 1.0);
        params = membership_settings(cfg);
        [isc, xs] = test_membership(Ws, params);
        memW = Ws(isc, :);
        ok = sum(isc) == stored_n(m) && isequal(memW, D.css_stack_per_mode{m});
        fprintf('CONTROL gen 1 %-20s members %d (stored %d), stack identical %d\n', ...
            names{m}, sum(isc), stored_n(m), isequal(memW, D.css_stack_per_mode{m}));
        if ~ok, error('silent_and_moving_members:gen1', 'Gen-1 control fails for %s.', names{m}); end
        F = flag_members(memW, xs(isc, :), params, FLOORS);
        g1(m) = pack(names{m}, 'gen1', NaN, nW, isc, F); %#ok<AGROW>
    end
end


% =========================================================================
function g500 = gen500_pass(pub, FLOORS)
    conds = [load_sensing(pub), load_spatial(pub)];
    R = load(fullfile(pub, 'css_radius_sensitivity.mat'), 'is_css_3', 'condition_names');
    for c = 1:numel(conds)
        params = conds(c).params;
        [isc, xs] = test_membership(conds(c).Wstack, params);
        n = size(conds(c).Wstack, 1);
        same = isequal(isc(:), R.is_css_3{c}(:)) && strcmp(R.condition_names{c}, conds(c).name);
        fprintf('CONTROL gen 500 %-20s members %4d (stored %4d), n %4d (stored %4d), per-wiring = radius file %d\n', ...
            conds(c).name, sum(isc), conds(c).stored_n_members, n, conds(c).stored_n, same);
        if sum(isc) ~= conds(c).stored_n_members || n ~= conds(c).stored_n || ~same
            error('silent_and_moving_members:gen500', 'Gen-500 control fails for %s.', conds(c).name);
        end
        F = flag_members(conds(c).Wstack(isc, :), xs(isc, :), params, FLOORS);
        g = pack(conds(c).name, conds(c).family, conds(c).success, n, isc, F);
        g.stored_occ = conds(c).stored_occ;
        g500(c) = g; %#ok<AGROW>
    end
end


% =========================================================================
function [isc, xs] = test_membership(Ws, params)
    n = size(Ws, 1); isc = false(n, 1); xs = nan(n, 4);
    parfor k = 1:n
        r = membership_test(reshape(Ws(k, :), 4, 4), params);
        isc(k) = r.is_css;
        xs(k, :) = r.x_star(:)';
    end
end


% =========================================================================
function F = flag_members(Wm, xs, params, FLOORS)
% The member's 600-step clone state, re-simulated with the membership test's pair settings and
% asserted equal to x_star; then the 20,000-step run without the stop on convergence.
    sp = sim_params_of(params);
    spL = sp; spL.developmentSteps = 20000; spL.stabilityThreshold = 0;
    n = size(Wm, 1); M = ones(4);
    x6 = nan(n, 4); xl = nan(n, 4);
    parfor k = 1:n
        W = reshape(Wm(k, :), 4, 4);
        [~, f6] = simulate_pair(M, W, [], M, W, [], sp);
        [~, fl] = simulate_pair(M, W, [], M, W, [], spL);
        x6(k, :) = f6(1:4)'; xl(k, :) = fl(1:4)';
    end
    assert(isequal(x6, xs), 'Re-simulated state differs from membership_test x_star.');
    F.x = x6; F.xlong = xl;
    F.xnorm  = sqrt(sum(x6.^2, 2));
    F.outmax = max(x6(:, 1:2), [], 2);
    F.rel    = sqrt(sum((xl - x6).^2, 2)) ./ F.xnorm;
    F.flip   = (xl(:, 1) - xl(:, 2)) <= 0;
    F.unconv = F.rel > 0.05 | F.flip;
    F.silent = F.outmax < FLOORS;                    % n x numel(FLOORS)
    F.either = F.silent | F.unconv;
end


% =========================================================================
function g = pack(name, family, success, n, isc, F)
    g = struct('name', name, 'family', family, 'success', success, 'n', n, ...
        'n_mem', sum(isc), 'is_css', isc, 'F', F, ...
        'n_silent', sum(F.silent, 1), 'n_unconv', sum(F.unconv), ...
        'n_either', sum(F.either, 1), 'n_flip', sum(F.flip), ...
        'med_outmax', median(F.outmax), 'stored_occ', NaN);
    g.occ = g.n_mem / n;
    g.occ_clean = (g.n_mem - g.n_either) / n;        % 1 x numel(FLOORS)
end


% =========================================================================
function out = summarize(g1, g500, FLOORS)
    print_table('GENERATION 1 (10,000 draws per architecture)', g1, FLOORS);
    print_table('GENERATION 500 (9 settings)', g500, FLOORS);
    out = standard(g1, g500, FLOORS);
end

function print_table(label, G, FLOORS)
    fprintf('\n=== %s ===\n', label);
    fprintf(['%-20s %5s %5s %6s %6s %6s %6s %6s %5s %9s %9s %9s %10s\n'], 'setting', 'n', ...
        'mem', 'SIL-3', 'SIL-2', 'UNCON', 'EIT-3', 'EIT-2', 'flip', 'members %', ...
        'clean-3 %', 'clean-2 %', 'med out');
    for c = 1:numel(G)
        g = G(c);
        fprintf('%-20s %5d %5d %6d %6d %6d %6d %6d %5d %9.2f %9.2f %9.2f %10.3g\n', g.name, ...
            g.n, g.n_mem, g.n_silent(1), g.n_silent(2), g.n_unconv, g.n_either(1), ...
            g.n_either(2), g.n_flip, 100 * g.occ, 100 * g.occ_clean(1), ...
            100 * g.occ_clean(2), g.med_outmax);
    end
    fprintf('(SIL-3/-2: max(x_C,x_D) < %g / %g; EIT = SILENT or UNCONVERGED; med out = median member max(x_C,x_D))\n', ...
        FLOORS(1), FLOORS(2));
end

function out = standard(g1, g500, FLOORS)
% Ranges of the share of members, its correlation with success (pooled and per family) and the
% architecture spread at generations 1 and 500, as stored and with flagged members removed.
    fprintf('\n=== Shares of members and r with success: as stored, and with flagged members removed ===\n');
    share1 = [g1.occ]'; share1c = reshape([g1.occ_clean], numel(FLOORS), [])';
    share5 = [g500.occ]'; share5c = reshape([g500.occ_clean], numel(FLOORS), [])';
    succ = [g500.success]'; fam = {g500.family}';
    fams = {'signaling', 'spatial'};
    cols = [share5, share5c]; tags = [{'published'}, arrayfun(@(f) sprintf('clean %g', f), FLOORS, 'UniformOutput', false)];
    o1 = [share1, share1c]; o5s = cols(strcmp(fam, 'signaling'), :);
    for j = 1:size(cols, 2)
        rp = corr(cols(:, j), succ);
        rf = cellfun(@(f) corr(cols(strcmp(fam, f), j), succ(strcmp(fam, f))), fams);
        s1 = 100 * (max(o1(:, j)) - min(o1(:, j))); s5 = 100 * (max(o5s(:, j)) - min(o5s(:, j)));
        e1 = 100 * (o1(4, j) - o1(1, j));  e5 = 100 * (o5s(4, j) - o5s(1, j));
        fprintf(['%-11s gen-1 %.2f-%.2f%%  gen-500 %.2f-%.2f%%  r pooled %+.4f  r sig/spa %+.4f %+.4f  ', ...
            'spread gen-1 %.2f pp gen-500 %.2f pp evolved %.1f%%  (SC-GS: %.2f / %.2f pp, %.1f%%)\n'], ...
            tags{j}, 100 * min(o1(:, j)), 100 * max(o1(:, j)), 100 * min(cols(:, j)), ...
            100 * max(cols(:, j)), rp, rf, s1, s5, 100 * (s5 - s1) / s5, e1, e5, 100 * (e5 - e1) / e5);
        out.tags{j} = tags{j}; out.r_pooled(j) = rp; out.r_family(:, j) = rf(:);
        out.spread1(j) = s1; out.spread500(j) = s5; out.evolved(j) = (s5 - s1) / s5;
    end
    out.occ1 = o1; out.occ500 = cols; out.success = succ; out.family = fam;
end


function c = make_setting(family, name, success, params, Wstack, stored_share, stored_n_members, stored_n)
    c = struct('family', family, 'name', name, 'success', success, ...
               'params', params, 'Wstack', Wstack, 'stored_occ', stored_share, ...
               'stored_n_members', stored_n_members, 'stored_n', stored_n);
end

function Wstack = subsample(results, n_sub)
    Wstack = zeros(0, 16);
    for r = 1:numel(results)
        B = results{r}.finalNetworkB;
        M = results{r}.finalNetworkM;
        nAg = size(B, 3);
        if isinf(n_sub), pick = 1:nAg; else, pick = randperm(nAg, min(n_sub, nAg)); end
        for a = pick
            W = double(M(:,:,a)) .* double(B(:,:,a));
            Wstack(end+1, :) = W(:)';                                       %#ok<AGROW>
        end
    end
end

function conds = load_sensing(pubDir)
    S = load(fullfile(pubDir, 'css_predicts_sensing_gen500.mat'));
    pats = {'graded_separate_100reps.mat', 'graded_contrast_100reps.mat', ...
            'saturating_separate_100reps.mat', 'saturating_contrast_100reps.mat'};
    rng(42);
    conds = [];
    for m = 1:numel(pats)
        files = dir(fullfile(pubDir, pats{m}));
        [~, idx] = max([files.datenum]);
        data = load(fullfile(pubDir, files(idx).name));
        Ws = subsample(data.results, 30);
        conds = [conds, make_setting('signaling', S.condition_names{m}, ...
                 S.success_empirical(m), membership_settings(data.config), Ws, ...
                 S.reach_occupancy(m), S.n_in_Vm(m), S.n_subsample(m))];  %#ok<AGROW>
    end
end

function conds = load_spatial(pubDir)
    S = load(fullfile(pubDir, 'css_predicts_spatial.mat'));
    tags = {'wellmixed', 'ring_k1', 'ring_k2', 'ring_k5', 'ring_k10'};
    conds = [];
    for c = 1:numel(tags)
        files = dir(fullfile(pubDir, sprintf('sweep_spatial_N50_%s_100reps_*.mat', tags{c})));
        [~, idx] = max([files.datenum]);
        data = load(fullfile(pubDir, files(idx).name));
        Ws = subsample(data.results, Inf);
        conds = [conds, make_setting('spatial', S.condition_names{c}, ...
                 S.success_empirical(c), membership_settings(data.config), Ws, ...
                 S.reach_occupancy(c), S.n_in_Vm(c), S.n_evolved_total(c))];  %#ok<AGROW>
    end
end
