function result = membership_test(W, params)
% MEMBERSHIP_TEST  The membership test: does a wiring belong to the conditional genotype network?
%
% Implements : SI S2.3, "Deterministic-CSS definition of the conditional genotype network", and
%              Methods, "Success, share of members and statistics": (C1) the cooperation
%              condition, the clone pair cooperates; (C2) the retaliation condition, no
%              single-weight mutant on the 16 x 11 grid (Delta from -3 sigma_mut to +3 sigma_mut,
%              Delta = 0 skipped, 160 mutants) earns more than the clone payoff R against the
%              wiring. Each pair develops with simulate_pair from x_j = 0.005 for up to 600 steps
%              (SI S2.3, "Developmental initial condition").
% Inputs     : W        4-by-4 wiring
%              params   struct: sensingMode ('separate' or 'contrast'; a data file's older spelling is read through stored_value),
%                       hillCoefficient (1 saturating, Inf graded), discreteOutput (false in the
%                       paper); optional fields and their defaults are set at the top of the code
% Outputs    : result   struct: .is_weak (C1 holds); .is_css (a member: C1 and C2); .max_pi_M
%                       (highest payoff of a grid mutant against the wiring); .x_star, .x_star_norm
%                       (the clone pair's settled state x*, ||x*||_2); .alpha, .alpha_p, .beta at
%                       x* (NaN unless C1 holds); .membership_fraction, .is_strong, .min_r_M,
%                       .worst_retal, .worst_retal_slot, .backend_used (explained where set)
% Seed       : sets its own seed per wiring and draw (css_draw_seed) and restores the caller's
%              random stream on return, so the verdict is a function of the wiring
%
% (C1) is decided by the sign of the readout x_C - x_D alone, with no floor on expression, so a
% clone pair whose outputs have both decayed towards zero is classed by the sign of what remains
% (SI S2.3, "Most random wirings that cooperate with their clone do not retaliate"). With T <= R,
% (C2) cannot fail (a mutant earns more than R only by collecting T) and membership is (C1) alone.

    % Defaults
    if ~isfield(params, 'michaelisConstant'), params.michaelisConstant = 1.0; end
    if ~isfield(params, 'degradationRate'),   params.degradationRate   = 0.1; end
    if ~isfield(params, 'timeStep'),          params.timeStep          = 1.0; end
    if ~isfield(params, 'developmentSteps'),  params.developmentSteps  = 600; end
    % sigma_mut is a standard deviation: the run's own mutationStdWeight where given, so the grid
    % matches the mutations the wirings evolved under; else 0.2, the paper's value (Table S1)
    if ~isfield(params, 'sigma_mut')
        if isfield(params, 'mutationStdWeight')
            params.sigma_mut = params.mutationStdWeight;
        else
            params.sigma_mut = 0.2;
        end
    end
    if ~isfield(params, 'cssMagGrid')
        params.cssMagGrid = linspace(-3*params.sigma_mut, 3*params.sigma_mut, 11);
    end
    % Route: 'simulate_pair' (the paper's) develops each pair with simulate_pair from the fixed
    % start x_j = 0.005. 'linear' integrates the joint pair with solve_fixed_point, from a random
    % start and to convergence, and is valid only for separate sensing with graded output: there
    % (C1) requires convergence and a mutant pair that does not converge counts as breaking (C2).
    % 'auto' takes 'linear' where it is valid and 'simulate_pair' elsewhere, so it mixes the two
    % rules across architectures.
    if ~isfield(params, 'backend'),           params.backend = 'simulate_pair'; end
    % Draws per wiring, each with its own seed (css_draw_seed); a draw matters only where the
    % start is random (the 'linear' route)
    if ~isfield(params, 'nDraws'),            params.nDraws  = 1; end
    % gradedC2: test the whole grid instead of stopping at the first mutant that earns more than R
    % (the scores are explained where they are computed); is_css is decided the same way either way
    if ~isfield(params, 'gradedC2'),          params.gradedC2 = false; end
    % Sensor gain g (Eq. 4), 1.0 in the paper
    if ~isfield(params, 'sensorIntegrationRate'), params.sensorIntegrationRate = 1.0; end
    % A constant input carried by the genotype: added to every network that carries the wiring
    % (the wiring, its clone, its grid mutants), not to the zero-wiring partner of is_strong;
    % 'simulate_pair' route only
    if ~isfield(params, 'geneInputBias'),     params.geneInputBias = zeros(4,1); end
    % A random start for the 'simulate_pair' route is not used in the paper and stops the run
    if isfield(params, 'randomStart') && params.randomStart
        error('membership_test:randomStart', 'params.randomStart is not supported; the test starts from x_j = 0.005');
    end

    % Route: 'auto' takes 'linear' only for separate sensing with graded output (see the defaults)
    if strcmp(params.backend, 'auto')
        is_linear_compat = isinf(params.hillCoefficient) && ...
                           ~params.discreteOutput && ...
                           strcmp(stored_value('sensingMode', params.sensingMode), 'separate');
        if is_linear_compat
            params.backend = 'linear';
        else
            params.backend = 'simulate_pair';
        end
    end

    % Result, filled in by the route below
    result.is_weak      = false;
    result.is_strong    = false;
    result.is_css       = false;
    result.membership_fraction = 0;
    result.alpha        = NaN;
    result.alpha_p      = NaN;
    result.beta         = NaN;
    result.x_star       = nan(4, 1);
    result.x_star_norm  = NaN;
    result.max_pi_M     = NaN;
    result.min_r_M      = NaN;
    result.worst_retal  = NaN;
    result.worst_retal_slot = NaN;
    result.backend_used = params.backend;

    % Mutant's payoff payoff_I_table(a_R + 1, a_M + 1), a_R the wiring's action and a_M the
    % mutant's (0 = D, 1 = C); params.payoffs = [T R P S] replaces the paper's T = 5, R = 3,
    % P = 1, S = 0 (Methods, "Game and evolution"). A mutant breaks (C2) when its payoff exceeds R.
    if isfield(params, 'payoffs') && ~isempty(params.payoffs)
        p4 = params.payoffs;
        T_p = p4(1); R_p = p4(2); P_p = p4(3); S_p = p4(4);
    else
        T_p = 5; R_p = 3; P_p = 1; S_p = 0;
    end
    payoff_I_table = [P_p, S_p; T_p, R_p];

    % Each draw is seeded from the wiring itself, so a wiring always gets the same verdict; the
    % caller's stream is saved here and restored after the last draw
    rng_state = rng();
    n_pass    = 0;
    first     = result;
    for draw = 1:params.nDraws
        rng(css_draw_seed(W, draw));
        r_draw = result;
        if strcmp(params.backend, 'linear')
            r_draw = css_linear_backend(W, params, r_draw, payoff_I_table);
        else
            r_draw = css_simulate_pair_backend(W, params, r_draw, payoff_I_table);
        end
        if draw == 1, first = r_draw; end
        n_pass = n_pass + double(r_draw.is_css);
    end
    rng(rng_state);

    % x_star, alpha and beta from the first draw; membership over all draws: with nDraws > 1 a
    % member (is_css) passed in more than half of the draws (membership_fraction)
    result = first;
    result.membership_fraction = n_pass / params.nDraws;
    result.is_css              = result.membership_fraction > 0.5;
end


% =========================================================================
function s = css_draw_seed(W, draw)
% Seed for one wiring and one draw, a deterministic function of the 16 weights and the draw
% number; the Mersenne Twister's initialization separates nearby seeds, so a wiring and its
% single-weight mutants get unrelated draws.
    mixer = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53];
    s = mod(round(abs(W(:)' * mixer(:)) * 1e6) + 7919 * draw, 2^31 - 2) + 1;
end


% =========================================================================
function result = css_linear_backend(W, params, result, payoff_I_table)
% 'linear' route: the joint pair integrated by solve_fixed_point; separate sensing with graded output only.

    mu = params.degradationRate;
    km = params.michaelisConstant;
    dt = params.timeStep;
    maxSteps = 1000;
    tol = 1e-10;

    % Coupling E: gene 3 (S_C) receives the partner's gene 1 (C), gene 4 (S_D) its gene 2 (D)
    E = zeros(4); E(3,1) = 1; E(4,2) = 1;

    s_zero4 = zeros(4, 1);
    s_zero8 = zeros(8, 1);
    s_def   = [0; 0; 0; 1];  % constant input on the zero-wiring partner's S_D (is_strong)

    grad_r_self_4 = [1; -1; 0; 0];       % readout x_C - x_D of the clone pair (4 genes)
    grad_r_V_8    = [1; -1; 0; 0; 0; 0; 0; 0];  % the wiring's readout in the joint pair (8 genes)
    grad_r_I_8    = [0;  0; 0; 0; 1; -1; 0; 0];  % the mutant's readout in the joint pair

    W_def = zeros(4);  % the zero-wiring partner of is_strong
    R_p   = payoff_I_table(2, 2);  % R, the clone payoff

    % (C1): the clone pair converges with readout > 0
    W_joint = W + E;
    [xJ, infoJ] = solve_fixed_point(W_joint, s_zero4, mu, km, dt, maxSteps, tol);
    if ~(infoJ.converged && (grad_r_self_4' * xJ > 0))
        return;  % fails (C1)
    end
    result.is_weak     = true;
    result.x_star      = xJ;
    result.x_star_norm = norm(xJ);

    % is_strong: the wiring's readout <= 0 against the zero-wiring partner
    W_block = [W, E; E, W_def];
    s_block = [s_zero4; s_def];
    [x_pair, infoP] = solve_fixed_point(W_block, s_block, mu, km, dt, maxSteps, tol);
    if infoP.converged && (grad_r_V_8' * x_pair <= 0)
        result.is_strong = true;
    end

    % (C2): every single-weight mutant W + Delta e_ij, 16 slots x the grid values, paired with
    % the wiring; a member if no mutant's payoff exceeds R
    max_pi_M   = -inf;
    any_invasion = false;

    for ii_pos = 1:4
        for jj_pos = 1:4
            for delta = params.cssMagGrid
                if delta == 0, continue; end

                W_mut = W;
                W_mut(ii_pos, jj_pos) = W_mut(ii_pos, jj_pos) + delta;

                W_block_mut = [W, E; E, W_mut];
                [x_mut, infoM] = solve_fixed_point(W_block_mut, s_zero8, mu, km, dt, maxSteps, tol);

                if ~infoM.converged
                    % A mutant pair that does not converge counts as a mutant earning 5, more than R
                    max_pi_M = max(max_pi_M, 5);
                    any_invasion = true;
                    break;
                end

                r_R_mut   = grad_r_V_8' * x_mut;
                r_M_mut   = grad_r_I_8' * x_mut;
                a_R_mut   = double(r_R_mut > 0);
                a_M_mut   = double(r_M_mut > 0);
                pi_M      = payoff_I_table(a_R_mut+1, a_M_mut+1);
                max_pi_M  = max(max_pi_M, pi_M);

                if pi_M > R_p
                    any_invasion = true;
                    break;
                end
            end
            if any_invasion, break; end
        end
        if any_invasion, break; end
    end

    result.max_pi_M = max_pi_M;
    if ~any_invasion
        result.is_css = true;
    end

    % alpha and beta at the clone pair's state, where (C1) holds
    if result.is_weak
        [result.alpha, result.beta] = compute_alpha_beta_linear(W, xJ, params);
    end
end


% =========================================================================
function result = css_simulate_pair_backend(W, params, result, payoff_I_table)
% 'simulate_pair' route: every pair developed by simulate_pair (see the defaults, "Route").

    M     = ones(4);
    K_vec = [];

    % Settings of every pair of the test; geneInputBias on both networks
    sim_params = struct( ...
        'sensingMode',       params.sensingMode, ...
        'hillCoefficient',   params.hillCoefficient, ...
        'discreteOutput',    params.discreteOutput, ...
        'michaelisConstant', params.michaelisConstant, ...
        'degradationRate',   params.degradationRate, ...
        'timeStep',          params.timeStep, ...
        'developmentSteps',  params.developmentSteps, ...
        'payoffDC',          0, ...
        'sensorIntegrationRate', params.sensorIntegrationRate, ...
        'externalSignal1',   params.geneInputBias(:), ...
        'externalSignal2',   params.geneInputBias(:));

    R_p = payoff_I_table(2, 2);  % R, the clone payoff

    % (C1): the wiring cooperates with its clone
    [outcome_clone, final_clone, ~] = simulate_pair(M, W, K_vec, M, W, K_vec, sim_params);
    is_weak = strcmp(outcome_clone, 'CC');
    if ~is_weak
        return;
    end
    result.is_weak     = true;
    result.x_star      = final_clone(1:4);
    result.x_star_norm = norm(result.x_star);

    % is_strong: whether the wiring defects against a partner with zero wiring and a constant
    % input of 1 on its S_D gene, which reaches no other gene, so that partner's outputs decay
    % towards zero; not condition (C2)
    sim_params_def = sim_params;
    sim_params_def.externalSignal2 = [0; 0; 0; 1];
    [outcome_def, ~, ~] = simulate_pair(M, W, K_vec, M, zeros(4), K_vec, sim_params_def);
    if outcome_def(1) ~= 'C'
        result.is_strong = true;
    end

    % (C2): every single-weight mutant on the grid, paired with the wiring
    max_pi_M   = -inf;
    min_r_M    = inf;      % lowest readout of a mutant against the wiring
    worst_retal = -inf;    % highest readout of the wiring against a defecting mutant
    worst_retal_slot = NaN;
    any_invasion = false;

    for ii_pos = 1:4
        for jj_pos = 1:4
            for delta = params.cssMagGrid
                if delta == 0, continue; end

                W_mut = W;
                W_mut(ii_pos, jj_pos) = W_mut(ii_pos, jj_pos) + delta;

                [outcome_mut, final_mut, ~] = simulate_pair(M, W, K_vec, M, W_mut, K_vec, sim_params);
                a_R_mut = double(outcome_mut(1) == 'C');
                a_M_mut = double(outcome_mut(2) == 'C');
                pi_M    = payoff_I_table(a_R_mut+1, a_M_mut+1);
                max_pi_M = max(max_pi_M, pi_M);

                % Continuous (C2) scores (gradedC2 true; the grid is then tested whole).
                % min_r_M is the lowest readout any grid mutant reaches against the wiring;
                % worst_retal is, over the mutants that defect, the highest readout of the
                % wiring itself (r_R < 0 is retaliation), its negative the retaliation margin
                % (SI S2.3, "The retaliation margin"); worst_retal_slot is the weight that
                % mutant changed, as the linear index i + 4(j - 1). The mutant is organism 2,
                % its readout final(5) - final(6); it earns more than R by defecting
                % (readout < 0) while the wiring cooperates.
                if params.gradedC2
                    r_M = final_mut(5) - final_mut(6);
                    if isfinite(r_M), min_r_M = min(min_r_M, r_M); end
                    % Over the mutants that defect: the wiring's own readout r_R (< 0 is
                    % retaliation); the highest is the closest it comes to being exploited
                    if a_M_mut == 0
                        r_R = final_mut(1) - final_mut(2);
                        if isfinite(r_R) && r_R > worst_retal
                            worst_retal = r_R;
                            worst_retal_slot = ii_pos + 4*(jj_pos - 1);
                        end
                    end
                end

                if pi_M > R_p
                    any_invasion = true;
                    if ~params.gradedC2, break; end
                end
            end
            if any_invasion && ~params.gradedC2, break; end
        end
        if any_invasion && ~params.gradedC2, break; end
    end

    result.max_pi_M = max_pi_M;
    if params.gradedC2 && isfinite(min_r_M), result.min_r_M = min_r_M; end
    if params.gradedC2 && isfinite(worst_retal)
        result.worst_retal = worst_retal;
        result.worst_retal_slot = worst_retal_slot;
    end
    if ~any_invasion
        result.is_css = true;
    end

    % alpha and beta by finite differences at the clone pair's state (compute_alpha_beta_simpair)
    if result.is_weak
        [result.alpha, result.beta, result.alpha_p] = compute_alpha_beta_simpair(W, result.x_star, sim_params);
    end
end


% =========================================================================
function [alpha, beta] = compute_alpha_beta_linear(W, x_star, params)
% 'linear' route: alpha by implicit differentiation of the fixed-point condition, beta by finite
% differences of the fixed point over the 16 weights.

    mu = params.degradationRate;
    km = params.michaelisConstant;
    E  = zeros(4); E(3,1) = 1; E(4,2) = 1;

    % At the fixed point mu x = f((W + E) x + s), s = 0, so
    % dx/ds = (mu I - diag(f') (W + E))^{-1} diag(f')
    u     = (W + E) * x_star;
    u_pos = max(u, 0);
    % f'(u) = K_M / (K_M + u)^2 for u >= 0, else 0
    f_prime = (km ./ (km + u_pos).^2) .* double(u >= 0);
    A = mu * eye(4) - diag(f_prime) * (W + E);
    if rcond(A) < 1e-10    % numerically singular: alpha and beta undefined
        alpha = NaN; beta = NaN;
        return;
    end
    J_alpha = A \ diag(f_prime);
    alpha   = norm(J_alpha, 'fro');

    % beta: finite differences of the clone fixed point over the 16 weights
    eps_fd = 1e-5;
    s_zero4 = zeros(4, 1);

    [x_base, info_base] = solve_fixed_point(W + E, s_zero4, mu, km, params.timeStep, 1000, 1e-10);
    if ~info_base.converged
        beta = NaN;
        return;
    end

    dx_dW_sq = 0;
    for ii = 1:4
        for jj = 1:4
            W_p = W; W_p(ii,jj) = W_p(ii,jj) + eps_fd;
            [x_p, info_p] = solve_fixed_point(W_p + E, s_zero4, mu, km, params.timeStep, 1000, 1e-10);
            if ~info_p.converged
                beta = NaN;
                return;
            end
            % A step that flips the sign of the readout leaves beta undefined
            r_base = x_base(1) - x_base(2);
            r_p    = x_p(1)    - x_p(2);
            if sign(r_base) ~= sign(r_p)
                beta = NaN;
                return;
            end
            dx_ij    = (x_p - x_base) / eps_fd;
            dx_dW_sq = dx_dW_sq + sum(dx_ij.^2);
        end
    end
    beta = sqrt(dx_dW_sq);
end


% =========================================================================
function [alpha, beta, alpha_p] = compute_alpha_beta_simpair(W, x_star, sim_params)
% 'simulate_pair' route: alpha, beta and alpha_p by finite differences of the clone pair's
% settled state. beta is the Frobenius norm of dx*/dW over the 16 weights (finite differences of
% 1e-4); alpha the Frobenius norm of dx*/ds over an input at each of the four genes (1e-3), and
% alpha_p the same over the two sensor genes only. Both alphas norm the whole state x*, not the
% readout-projected alpha of SI Eq. S6, and no value the paper prints comes from them.

    eps_s = 1e-3;  % perturbation of the input to each gene (alpha)
    eps_W = 1e-4;  % perturbation for weight sensitivity (beta)
    M     = ones(4);
    K_vec = [];

    base = x_star;

    % alpha: an input at each of the four genes, since beta is taken over all 16 weights and
    % beta = alpha ||x*|| (Eq. 2) holds when both are taken over every gene; alpha_p over the
    % two sensor genes (3, 4) only, the responsiveness to the partner's signal (SI S1.5)
    dx_ds_sq = 0; dx_ds_sq_p = 0;
    for k = 1:4
        sp = sim_params;
        sp.externalSignal1    = sim_params.externalSignal1;          % geneInputBias
        sp.externalSignal1(k) = sp.externalSignal1(k) + eps_s;
        [~, fs_p, ~] = simulate_pair(M, W, K_vec, M, W, K_vec, sp);
        dx_k     = (fs_p(1:4) - base) / eps_s;
        dx_ds_sq = dx_ds_sq + sum(dx_k.^2);
        if k >= 3, dx_ds_sq_p = dx_ds_sq_p + sum(dx_k.^2); end
    end
    alpha   = sqrt(dx_ds_sq);
    alpha_p = sqrt(dx_ds_sq_p);

    % beta: finite differences over the 16 weights, both networks of the clone pair changed
    dx_dW_sq = 0;
    for ii = 1:4
        for jj = 1:4
            W_p = W; W_p(ii,jj) = W_p(ii,jj) + eps_W;
            [~, fs_p, ~] = simulate_pair(M, W_p, K_vec, M, W_p, K_vec, sim_params);
            % A step that flips the action leaves beta undefined
            base_coop = (base(1) > base(2));
            new_coop  = (fs_p(1) > fs_p(2));
            if base_coop ~= new_coop
                beta = NaN;
                return;
            end
            dx_ij    = (fs_p(1:4) - base) / eps_W;
            dx_dW_sq = dx_dW_sq + sum(dx_ij.^2);
        end
    end
    beta = sqrt(dx_dW_sq);
end
