function scope_checks_gen1(nW, withDenominator)
% SCOPE_CHECKS_GEN1  On the 10,000 random wirings per architecture behind the generation-1 share of
% members: which members also resist every single topology toggle, and which wirings have a clone
% state that is silent or still moving?
%
% Produces : the input that scope_checks_summary reports (Methods, "Success, share of members and
%            statistics"; SI S2.3; Tables S2 and S3)
% Inputs   : publication/css_predicts_sensing.mat   (gen1_share_of_members; the stored members)
%            publication/css_c1_c2_split.mat   (cooperation_condition_shares; the stored shares
%              meeting the cooperation condition (C1)); both for the check
%            nW               wirings per architecture (default 10,000, the paper's run; a smaller
%                             nW takes the first nW of the same draws and skips the check)
%            withDenominator  flag every wiring, not only those meeting (C1) (default true)
% Outputs  : publication/scope_checks/gen1_toggle_silent.mat   (G, per architecture and wiring:
%              (C1) `is_weak`, membership `is_css`, the toggle result and the flags F); one console
%              line per architecture. The folder (the scope checks' output folder) must exist
%              (scope_checks_gen500 creates it).
% Run      : scope_checks_gen1   (parfor on the current pool; nW = 3 takes about 1 min; the full
%              run is not timed)
% Seed     : rng(2026) (MATLAB's default generator) before the 10,000 draws randn(4) * 0.1, one
%            set for all four architectures; membership_test seeds each wiring itself
%
% A full run stops unless its (C1) count, member count and member wirings equal the stored ones
% (the members of css_predicts_sensing.mat: 55 / 219 / 109 / 301). Toggles: toggle_stable on each member with all 16 connections present, so
% every toggle removes one. Flags, on the clone state the membership test reads (600 steps from the
% fixed start): silent, max(x_C, x_D) < 1e-3; still moving, rerun from the same start for 20,000
% steps without the stop on convergence, the state moves by more than 5% or the readout x_C - x_D
% changes sign.

    if nargin < 1 || isempty(nW), nW = 10000; end
    if nargin < 2 || isempty(withDenominator), withDenominator = true; end
    pub = fullfile(grn_data_file('publication'));
    D = load(fullfile(pub, 'css_predicts_sensing.mat'), 'css_stack_per_mode');
    C = load(fullfile(pub, 'css_c1_c2_split.mat'), 'c1_gen1');
    names = {'graded separate', 'graded contrast', 'saturating separate', 'saturating contrast'};
    cfgs  = {{'separate', Inf}, {'contrast', Inf}, {'separate', 1}, {'contrast', 1}};
    full = nW == 10000;

    rng(2026);
    Ws = zeros(nW, 16);
    for i = 1:nW, W = randn(4) * 0.1; Ws(i, :) = W(:)'; end
    for m = 1:4
        t0 = tic;
        params = struct('sensingMode', cfgs{m}{1}, 'hillCoefficient', cfgs{m}{2}, ...
            'discreteOutput', false, 'michaelisConstant', 1.0, 'degradationRate', 0.1, ...
            'timeStep', 1.0, 'developmentSteps', 600, 'backend', 'simulate_pair', ...
            'sigma_mut', 0.2, 'cssMagGrid', linspace(-3*0.2, 3*0.2, 11));
        [c1, is_member, xs] = test_membership(Ws, params);
        if full, control(m, names{m}, c1, is_member, Ws, D, C); end
        [tog, cls] = toggle_members(Ws, is_member, params);
        F = flags(Ws, c1, xs, params, withDenominator);
        G(m) = struct('name', names{m}, 'n', nW, 'is_weak', c1, 'is_css', is_member, ...
            'toggle_ok', tog, 'toggle_cls', cls, 'F', F, 'seconds', toc(t0)); %#ok<AGROW>
        fprintf('%-20s C1 %d  members %d  toggle-fail %d  C1 silent %d unconv %d either %d  %.0fs\n', ...
            names{m}, sum(c1), sum(is_member), sum(is_member & ~tog), sum(c1 & F.silent), sum(c1 & F.unconv), ...
            sum(c1 & (F.silent | F.unconv)), G(m).seconds);
    end
    save(fullfile(pub, 'scope_checks', 'gen1_toggle_silent.mat'), 'G', 'nW', 'withDenominator');
end


% =========================================================================
function [c1, is_member, xs] = test_membership(Ws, params)
    n = size(Ws, 1); c1 = false(n, 1); is_member = c1; xs = nan(n, 4);
    parfor k = 1:n
        r = membership_test(reshape(Ws(k, :), 4, 4), params);
        c1(k) = r.is_weak; is_member(k) = r.is_css; xs(k, :) = r.x_star(:)';
    end
end

function control(m, name, c1, is_member, Ws, D, C)
% The (C1) count against css_c1_c2_split.mat, the members against the stored member wirings.
    nStored = size(D.css_stack_per_mode{m}, 1);
    okC1 = sum(c1) == round(C.c1_gen1(m) * 10000);
    okMember = sum(is_member) == nStored && isequal(Ws(is_member, :), D.css_stack_per_mode{m});
    fprintf('CONTROL %-20s C1 %d (stored %d)  members %d (stored %d), stack identical %d\n', ...
        name, sum(c1), round(C.c1_gen1(m) * 10000), sum(is_member), nStored, okMember);
    if ~(okC1 && okMember), error('scope_checks_gen1:control', 'Control fails for %s.', name); end
end

function [tog, cls] = toggle_members(Ws, is_member, params)
    n = size(Ws, 1); tog = true(n, 1); cls = zeros(n, 2);
    idx = find(is_member);
    t = true(numel(idx), 1); c = zeros(numel(idx), 2);
    parfor q = 1:numel(idx)
        W = reshape(Ws(idx(q), :), 4, 4);
        [t(q), c(q, :)] = toggle_stable(ones(4), W, params, 600);
    end
    tog(idx) = t; cls(idx, :) = c;
end


% =========================================================================
function F = flags(Ws, c1, xs, params, withDenominator)
% The silent and still-moving flags (silent_and_moving_members's), on the wirings that satisfy (C1),
% or on every wiring with withDenominator. The 600-step clone state is re-simulated here and must
% equal membership_test's x_star where (C1) holds.
    sp = sim_params_of(params);
    spL = sp; spL.developmentSteps = 20000; spL.stabilityThreshold = 0;
    n = size(Ws, 1);
    todo = c1;
    if withDenominator, todo = true(n, 1); end
    idx = find(todo);
    x6 = nan(numel(idx), 4); xl = x6;
    parfor q = 1:numel(idx)
        W = reshape(Ws(idx(q), :), 4, 4); M = ones(4);
        [~, f6] = simulate_pair(M, W, [], M, W, [], sp);
        [~, fl] = simulate_pair(M, W, [], M, W, [], spL);
        x6(q, :) = f6(1:4)'; xl(q, :) = fl(1:4)';
    end
    assert(isequal(x6(c1(idx), :), xs(idx(c1(idx)), :)), 'Re-simulated state differs from x_star.');
    X6 = nan(n, 4); XL = X6; X6(idx, :) = x6; XL(idx, :) = xl;
    F.flagged = todo;
    F.x = X6; F.xlong = XL;
    F.outmax = max(X6(:, 1:2), [], 2);
    F.rel = sqrt(sum((XL - X6).^2, 2)) ./ sqrt(sum(X6.^2, 2));
    r6 = X6(:, 1) - X6(:, 2); rl = XL(:, 1) - XL(:, 2);
    F.flip = (r6 > 0 & rl <= 0) | (r6 <= 0 & rl > 0);
    F.unconv = todo & (F.rel > 0.05 | F.flip);
    F.silent = todo & F.outmax < 1e-3;
end
