function [a1, a2] = encounter_batch(W1, W2, X1, X2, sp, steps, snapMode)
% ENCOUNTER_BATCH  Plays many encounters at once: in each pair the two networks develop together,
% each sensing the other's outputs, and a network cooperates if x_C > x_D once settled.
%
% Implements : one encounter of Methods, "Model" (dynamics and sensing equations, Eq. 3 and
%              Eq. 4, Euler steps), for both ways the paper develops a pair: the runs'
%              encounter rule, as evolve_population plays it, and the membership test (SI S2.3,
%              "Developmental initial condition")
% Inputs     : W1, W2     4 x 4 x P wirings, first and second network of each pair
%              X1, X2     4 x P start states
%              sp         struct: degradationRate, timeStep, michaelisConstant,
%                         sensorIntegrationRate (the gain g), hillCoefficient (1 saturating,
%                         Inf graded), sensingMode ('separate' or 'contrast'; older spellings
%                         via stored_value); discreteOutput must be false
%              steps      development steps (the step limit)
%              snapMode   'encounter_rule' (the runs' rule): a network's action is read at the last
%                         convergence check it passes; 'simpair' (the membership test's rule,
%                         as in simulate_pair): at the first it passes alone, or where both pass
% Outputs    : a1, a2   1 x P logical, true = cooperate
% Seed       : none (no random draws; the caller supplies the start states)
%
% The caller sets the rule. Runs: each network starts from X ~ U(0, 0.01), 300 steps, 'encounter_rule'.
% Membership test: X = 0.005, 600 steps, 'simpair'.

    P = size(X1, 2);
    mu = sp.degradationRate; dt = sp.timeStep; km = sp.michaelisConstant;
    ir = sp.sensorIntegrationRate; hill = sp.hillCoefficient; useHill = isfinite(hill);
    useSeparate = strcmp(stored_value('sensingMode', sp.sensingMode), 'separate');
    snapMode = stored_value('snapMode', snapMode);
    assert(any(strcmp(snapMode, {'encounter_rule', 'simpair'})), 'snapMode must be encounter_rule or simpair');
    enc_rule = strcmp(snapMode, 'encounter_rule');
    assert(~sp.discreteOutput, 'discreteOutput not supported');
    xs1 = X1; xs2 = X2; st1 = false(1, P); st2 = false(1, P);
    aIdx = 1:P; x1 = X1; x2 = X2; x1o = x1; x2o = x2; last = 0;
    for t = 1:steps
        % One Euler step of dx/dt = f(W x + S) - mu x with f(u) = u / (K_M + u), x kept >= 0
        [in1, in2] = pair_inputs(W1, W2, x1, x2, useHill, hill, useSeparate, ir);
        x1 = max(x1 + dt * (in1 ./ (km + in1) - mu * x1), 0);
        x2 = max(x2 + dt * (in2 ./ (km + in2) - mu * x2), 0);
        % Every 40 steps a network has converged when ||x(t) - x(t-40)||_2 / max_i x_i < 0.01
        % (compared squared); a pair leaves the batch once both networks have converged
        if t - last < 40, continue; end
        c1 = sum((x1 - x1o).^2, 1) ./ max(max(x1, [], 1), 1e-9).^2 < 1e-4;
        c2 = sum((x2 - x2o).^2, 1) ./ max(max(x2, [], 1), 1e-9).^2 < 1e-4;
        b = c1 & c2; n1 = c1 & ~b; n2 = c2 & ~b;
        if ~enc_rule, n1 = n1 & ~st1(aIdx); n2 = n2 & ~st2(aIdx); end
        s1 = b | n1; s2 = b | n2;
        xs1(:, aIdx(s1)) = x1(:, s1); st1(aIdx(s1)) = true;
        xs2(:, aIdx(s2)) = x2(:, s2); st2(aIdx(s2)) = true;
        if any(b)
            k = ~b; aIdx = aIdx(k); W1 = W1(:,:,k); W2 = W2(:,:,k); x1 = x1(:, k); x2 = x2(:, k);
            if isempty(aIdx), break; end
        end
        x1o = x1; x2o = x2; last = t;
    end
    % A network that never converged takes its state at the last step
    u1 = ~st1(aIdx); u2 = ~st2(aIdx);
    xs1(:, aIdx(u1)) = x1(:, u1); xs2(:, aIdx(u2)) = x2(:, u2);
    a1 = xs1(1,:) > xs1(2,:); a2 = xs2(1,:) > xs2(2,:);
end


function [in1, in2] = pair_inputs(W1, W2, x1, x2, useHill, hill, useSeparate, ir)
% Each network's rectified input W x + S. The partner's output is x_C, x_D (graded) or
% x_C^n / (x_C^n + x_D^n) and its complement (saturating, n = 1); separate sensing adds it to
% S_C and S_D on two channels, contrast sensing adds g (C - D) to S_C and its negative to S_D.
    in1 = reshape(pagemtimes(W1, reshape(x1, 4, 1, [])), 4, []);
    in2 = reshape(pagemtimes(W2, reshape(x2, 4, 1, [])), 4, []);
    if useHill
        oC = x2(1,:).^hill; oD = x2(2,:).^hill; d = oC + oD + 1e-12; cF = oC ./ d; dF = oD ./ d;
        fC = x1(1,:).^hill; fD = x1(2,:).^hill; d = fC + fD + 1e-12; cO = fC ./ d; dO = fD ./ d;
    else
        cF = x2(1,:); dF = x2(2,:); cO = x1(1,:); dO = x1(2,:);
    end
    if useSeparate
        in1(3,:) = in1(3,:) + ir * cF; in1(4,:) = in1(4,:) + ir * dF;
        in2(3,:) = in2(3,:) + ir * cO; in2(4,:) = in2(4,:) + ir * dO;
    else
        s1 = ir * (cF - dF); s2 = ir * (cO - dO);
        in1(3,:) = in1(3,:) + s1; in1(4,:) = in1(4,:) - s1;
        in2(3,:) = in2(3,:) + s2; in2(4,:) = in2(4,:) - s2;
    end
    in1 = max(in1, 0); in2 = max(in2, 0);
end
