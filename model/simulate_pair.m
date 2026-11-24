function [outcome, final_state, coop_trace] = simulate_pair(M1, B1, K1, M2, B2, K2, params)
% SIMULATE_PAIR  Develops two networks that meet as one coupled system and reads each one's action.
%
% Implements : Eq. 3 and Eq. 4 (Methods, "Model") for one encounter, with the stop on convergence
%              of SI Extended Materials and Methods, "Developmental integration and convergence";
%              the membership test (membership_test) and most analyses develop pairs here
% Inputs     : M1, M2   4-by-4 connectivity matrices M (M(i,j) = 1 when gene j regulates gene i)
%              B1, B2   4-by-4 signed weights B, so that W = M .* B
%              K1, K2   per-gene K_M (4 values), or [] for params.michaelisConstant on every gene
%              params   struct: sensingMode ('separate' or 'contrast' sensing; older spellings via stored_value),
%                       hillCoefficient (1 saturating output; Inf or absent graded output),
%                       discreteOutput (true passes the partner's action C > D as 0/1),
%                       michaelisConstant, degradationRate mu, timeStep Delta t, developmentSteps
%                       (the step limit, set by the caller: 600 in the membership test, 300 in
%                       the runs); optional fields, with their defaults, where they are read
% Outputs    : outcome      'CC', 'CD', 'DC' or 'DD': the action of organism 1, then of organism 2
%              final_state  expression [x_C x_D x_SC x_SD] each action is read from, organism 1
%                           then organism 2 (8 values)
%              coop_trace   x_C of organism 1 and of organism 2 at each step (steps-by-2)
% Seed       : none (deterministic from the starting state)
%
% Both organisms start by default at x_j = 0.005, the membership test's start (SI S2.3,
% "Developmental initial condition"); the runs start from x_j ~ U(0, 0.01) (Methods, "Model"),
% which a caller passes as initialState1/2. An organism's action is read at the first
% convergence check it passes alone, or where both pass (the rule is beside the check).

    % Parameters
    nGenes   = size(M1, 1);
    maxSteps = params.developmentSteps;
    dt       = params.timeStep;
    mu       = params.degradationRate;

    % Optional fields, with defaults: sensor gain g, convergence threshold and check interval,
    % positions of C, D, S_C and S_D (genes 1 to 4)
    intRate   = getf(params, 'sensorIntegrationRate', 1);
    threshold = getf(params, 'stabilityThreshold',    0.01);
    interval  = getf(params, 'stabilityCheckInterval', 40);
    gC = getf(params, 'geneCooperator', 1);
    gD = getf(params, 'geneDefector',   2);
    gS1 = getf(params, 'geneSensor1',   3);
    gS2 = getf(params, 'geneSensor2',   4);

    % Constant input added to every gene's input at each step (externalSignal1 for organism 1,
    % externalSignal2 for organism 2; default 0)
    extSig1 = getf(params, 'externalSignal1', zeros(nGenes, 1));
    extSig2 = getf(params, 'externalSignal2', zeros(nGenes, 1));
    extSig1 = extSig1(:);
    extSig2 = extSig2(:);
    if length(extSig1) ~= nGenes
        error('externalSignal1 must be %dx1; got %dx1', nGenes, length(extSig1));
    end
    if length(extSig2) ~= nGenes
        error('externalSignal2 must be %dx1; got %dx1', nGenes, length(extSig2));
    end

    useSeparate = strcmp(stored_value('sensingMode', params.sensingMode), 'separate');
    useDiscrete = isfield(params, 'discreteOutput') && params.discreteOutput;

    hillCoeff = Inf;
    if isfield(params, 'hillCoefficient'), hillCoeff = params.hillCoefficient; end
    useHill = isfinite(hillCoeff);

    % K_M per gene: the vector supplied, else params.michaelisConstant on every gene
    if nargin >= 3 && ~isempty(K1)
        km1 = double(K1(:));
    else
        km1 = ones(nGenes, 1) * params.michaelisConstant;
    end
    if nargin >= 6 && ~isempty(K2)
        km2 = double(K2(:));
    else
        km2 = ones(nGenes, 1) * params.michaelisConstant;
    end

    % Regulatory matrix W = M .* B (Methods, "Model")
    W1 = double(M1) .* double(B1);
    W2 = double(M2) .* double(B2);

    thresholdSq = threshold^2;

    % Starting state (initialState1, initialState2): default 0.005 on every gene (header)
    x1 = getf(params, 'initialState1', 0.005 * ones(nGenes, 1));
    x2 = getf(params, 'initialState2', 0.005 * ones(nGenes, 1));
    x1 = x1(:); x2 = x2(:);
    if numel(x1) ~= nGenes || numel(x2) ~= nGenes
        error('initialState1/2 must be %dx1', nGenes);
    end

    x1_old = x1;          % state at the previous convergence check
    x2_old = x2;
    stable1 = false;
    stable2 = false;
    x1_snap = x1;         % state each organism's action is read from
    x2_snap = x2;
    lastCheck = 0;

    coop_trace = zeros(maxSteps, 2);

    active = true;

    for t = 1:maxSteps
        coop_trace(t, 1) = x1(gC);
        coop_trace(t, 2) = x2(gC);

        if ~active
            continue;
        end

        % Regulatory input sum_j W_ij x_j
        inputF = W1 * x1;
        inputO = W2 * x2;

        % What each organism's sensors receive of its partner's outputs (Methods, "Model"):
        % 0/1 action (discreteOutput), x_C^n / (x_C^n + x_D^n) and its complement (saturating
        % output, n = hillCoefficient), or x_C and x_D as they are (graded output); 1e-12 keeps
        % the saturating output defined when x_C = x_D = 0
        if useDiscrete
            cur_opp_CF = double(x2(gC) > x2(gD));
            cur_opp_DF = 1.0 - cur_opp_CF;
            cur_opp_CO = double(x1(gC) > x1(gD));
            cur_opp_DO = 1.0 - cur_opp_CO;
        elseif useHill
            oC_n = x2(gC)^hillCoeff;
            oD_n = x2(gD)^hillCoeff;
            denom_o = oC_n + oD_n + 1e-12;
            cur_opp_CF = oC_n / denom_o;
            cur_opp_DF = oD_n / denom_o;
            fC_n = x1(gC)^hillCoeff;
            fD_n = x1(gD)^hillCoeff;
            denom_f = fC_n + fD_n + 1e-12;
            cur_opp_CO = fC_n / denom_f;
            cur_opp_DO = fD_n / denom_f;
        else
            cur_opp_CF = x2(gC);
            cur_opp_DF = x2(gD);
            cur_opp_CO = x1(gC);
            cur_opp_DO = x1(gD);
        end

        % Sensor input (Eq. 4): separate sensing S_C = g x_C, S_D = g x_D; contrast sensing
        % S_C = -S_D = g (x_C - x_D)
        if useSeparate
            inputF(gS1) = inputF(gS1) + intRate * cur_opp_CF;
            inputF(gS2) = inputF(gS2) + intRate * cur_opp_DF;
            inputO(gS1) = inputO(gS1) + intRate * cur_opp_CO;
            inputO(gS2) = inputO(gS2) + intRate * cur_opp_DO;
        else
            sigF = intRate * (cur_opp_CF - cur_opp_DF);
            sigO = intRate * (cur_opp_CO - cur_opp_DO);
            inputF(gS1) = inputF(gS1) + sigF;
            inputF(gS2) = inputF(gS2) - sigF;
            inputO(gS1) = inputO(gS1) + sigO;
            inputO(gS2) = inputO(gS2) - sigO;
        end

        % Constant input
        inputF = inputF + extSig1;
        inputO = inputO + extSig2;

        % Production f(u) = max(u, 0) / (K_M + max(u, 0)) (Eq. 3)
        inputF = max(inputF, 0);
        inputO = max(inputO, 0);
        prodF = inputF ./ (km1 + inputF);
        prodO = inputO ./ (km2 + inputO);

        % Euler step of dx_i/dt = f(u_i) - mu x_i, expression kept non-negative
        x1 = x1 + dt * (prodF - mu * x1);
        x2 = x2 + dt * (prodO - mu * x2);
        x1 = max(x1, 0);
        x2 = max(x2, 0);

        % Stop on convergence. Every 40 steps each organism is tested against
        % ||x(t) - x(t-40)||_2 / max_i x_i(t) < 0.01, in squares; 1e-9 keeps it defined when every
        % x_i = 0. An organism's action is read from its state at the first check it passes on
        % its own; it keeps developing while its partner has not converged, and when both pass
        % the same check, both are read there instead and development stops. An organism that
        % never passes is read at the step limit.
        if t - lastCheck >= interval
            maxF = max(max(x1), 1e-9);
            maxO = max(max(x2), 1e-9);
            diffF = sum((x1 - x1_old).^2) / maxF^2;
            diffO = sum((x2 - x2_old).^2) / maxO^2;
            convF = diffF < thresholdSq;
            convO = diffO < thresholdSq;

            if convF && convO
                % Both converged: read both here, stop developing, hold the trace
                x1_snap = x1;
                x2_snap = x2;
                stable1 = true;
                stable2 = true;
                active = false;
                for tt = t+1:maxSteps
                    coop_trace(tt, 1) = x1(gC);
                    coop_trace(tt, 2) = x2(gC);
                end
            else
                % One converged: read it at this first check it passes
                if convF && ~stable1
                    x1_snap = x1;
                    stable1 = true;
                end
                if convO && ~stable2
                    x2_snap = x2;
                    stable2 = true;
                end
            end
            x1_old = x1;
            x2_old = x2;
            lastCheck = t;
        end
    end

    % An organism that never converged is read at the step limit
    if ~stable1, x1_snap = x1; end
    if ~stable2, x2_snap = x2; end

    % Action: cooperate when x_C > x_D, otherwise defect (a tie is defection)
    coop1 = x1_snap(gC) > x1_snap(gD);
    coop2 = x2_snap(gC) > x2_snap(gD);

    if coop1 && coop2
        outcome = 'CC';
    elseif coop1 && ~coop2
        outcome = 'CD';
    elseif ~coop1 && coop2
        outcome = 'DC';
    else
        outcome = 'DD';
    end

    final_state = [x1_snap(:); x2_snap(:)];
end

function v = getf(s, name, default)
% A field of s, or the default when the field is absent or empty
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
