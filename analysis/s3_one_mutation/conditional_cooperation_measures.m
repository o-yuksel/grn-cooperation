function varargout = conditional_cooperation_measures(mode, varargin)
% CONDITIONAL_COOPERATION_MEASURES  Builds the networks the response dial starts from (SI S3.5):
% 500 saturating-contrast networks from the stored runs, 30 random mutations of each, and how each
% network reacts to 80 partners from the same populations. Also holds the constants K and the
% helpers that the response-dial scripts share.
%
% Produces : as input only, for main text, "Conditional cooperation stays within a mutation's
%            reach of the threshold"; Methods, "Mutation tests in the game"; Fig. 4; SI S3.4;
%            SI S3.5; Table S8. The 'np_*' phases (np: natural partners) run response_dial_networks, whose stored
%            networks and reactions response_dial reads; 'np_handles' hands K and the helpers to
%            response_dial and dial_loss_shares, and through them to fig4b_retaliation_from_rest
% Inputs   : mode, the phase (for 'np_analyze' also 'full' or 'timing_run', also accepted as 'pilot')
%            publication/<architecture>_100reps.mat   (the run configs; run_signaling_architectures);
%            the other stored inputs of the 'np_*' phases are listed in response_dial_networks
% Outputs  : see response_dial_networks; 'np_handles' returns {K, helpers, its helpers}
% Run      : conditional_cooperation_measures('np_pool')   draws the networks, their mutations and
%              the partners (about 2 min); then, in order, 'np_check1' (two of
%              mutate_with_refusal's encounter calls agree on 200 network-partner pairs),
%              'np_check2' (the copied mutation draws reproduce 200 stored census draws) and
%              'np_timing_run' (a trial on networks 1 to 50 that projects the full time), these three
%              about 3 min;
%              then 'np_score' (all 500; about 28 min on 3 workers; runs only if both checks
%              passed). 'np_analyze' prints checks and measures no paper item uses.
% Seed     : K.seed = 20261028; each draw set uses rng(K.seed + its own fixed offset, 'twister')

    K = consts();
    if ~isfolder(K.out), mkdir(K.out); end
    switch mode
        case {'np_pool', 'np_check1', 'np_check2', 'np_timing_run', 'np_score', 'np_analyze'}
            response_dial_networks(mode, K, helpers(), varargin{:});
        case 'np_handles'
            [Kn, NP] = response_dial_networks(mode, K, helpers()); varargout = {Kn, helpers(), NP};
        otherwise, error('unknown mode %s', mode);
    end
end


function H = helpers()
% Handles to the local helpers that response_dial_networks and response_dial call.
    H = struct('load_configs', @load_configs, 'boot_fit', @boot_fit, 'plant', @plant);
end


function K = consts()
    % K.seed and K.nBoot are read; K.filt, K.chunk, K.nWorkers, K.budgetH, K.eChunk, K.eBudgetH and
    % K.eDir are read by no script in the release
    % (response_dial_networks takes its pool size as an argument, default 3).
    K.pub = grn_data_file('publication');
    K.out = fullfile(K.pub, 'conditionality_score');
    K.archs = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'};
    K.cls = {'member', 'cooperate_only'};
    K.seed = 20261028; K.filt = 0.99; K.nBoot = 2000; K.chunk = 400; K.nWorkers = 4; K.budgetH = 3;
    K.eChunk = 5000; K.eBudgetH = 2; K.eDir = fullfile(K.out, 'encounter_rule');
end


function cfg = load_configs(K)
    cfg = cell(4, 1);
    for m = 1:4
        d = load(fullfile(K.pub, sprintf('%s_100reps.mat', K.archs{m})), 'config'); cfg{m} = d.config;
    end
end


% =========================================================================
% FITTING (boot_fit, plant: used by response_dial_networks' checks)
% =========================================================================
function v = wm(x, w)
    v = sum(x .* w) / sum(w);
end


function F = boot_fit(T, Xf, K, seed)
% Weighted logistic fits (binomial k of 30 per row, clone-count weights) for each design in Xf,
% point estimates and K.nBoot resamples of replicates (cluster = architecture x replicate),
% resampled within architecture; all designs share each resample.
    k = sum(T.lost, 2); n = size(T.lost, 2) * ones(size(k));
    [keys, ~, loc] = unique([T.arch T.rep], 'rows'); nJ = numel(Xf);
    X = cellfun(@(f) f(T), Xf, 'uni', 0);
    F.b = cellfun(@(x) logit_fit(x, k, n, T.count), X, 'uni', 0);
    F.B = cellfun(@(x) nan(K.nBoot, size(x, 2)), X, 'uni', 0);
    rs = rng; rng(seed, 'twister');
    for b = 1:K.nBoot
        mult = accumarray(strat_idx(keys(:, 1)), 1, [size(keys, 1) 1]);
        w = T.count .* mult(loc);
        for j = 1:nJ, F.B{j}(b, :) = logit_fit(X{j}, k, n, w)'; end
    end
    rng(rs);
    F.nFail = sum(any(isnan([F.B{:}]), 2));
end


function b = logit_fit(X, k, n, w)
% Newton-Raphson for the weighted binomial logit; NaN on non-convergence or separation.
    keep = w > 0; X = X(keep, :); k = k(keep); n = n(keep); w = w(keep);
    p = size(X, 2); b = zeros(p, 1);
    r0 = sum(w .* k) / sum(w .* n);
    if r0 <= 0 || r0 >= 1 || rank(X) < p, b = nan(p, 1); return; end
    b(1) = log(r0 / (1 - r0));
    for it = 1:100
        mu = 1 ./ (1 + exp(-X * b));
        g = X' * (w .* (k - n .* mu)); H = X' * (X .* (w .* n .* mu .* (1 - mu)));
        step = H \ g; b = b + step;
        if ~all(isfinite(b)) || max(abs(b)) > 1e3, b = nan(p, 1); return; end
        if max(abs(step)) < 1e-9, return; end
    end
    b = nan(p, 1);
end


function T = plant(T, orPer01, seed, unit)
% Replace the flags with draws from logit p = b0 + log(orPer01)/unit * score; b0 matches the
% clone-weighted observed rate of the same rows. unit (default 0.1): the score step the planted
% odds ratio is per (response_dial_networks passes 1 for a score in SD units).
    if nargin < 4, unit = 0.1; end
    b1 = log(orPer01) / unit; target = wm(mean(T.lost, 2), T.count);
    f = @(b0) wm(1 ./ (1 + exp(-(b0 + b1 * T.score))), T.count) - target;
    b0 = fzero(f, [-60 60]);
    p = 1 ./ (1 + exp(-(b0 + b1 * T.score)));
    rs = rng; rng(seed, 'twister'); T.lost = rand(numel(p), size(T.lost, 2)) < p; rng(rs);
end
