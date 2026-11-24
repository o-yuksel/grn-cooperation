function [x1, x2, it, ok] = rest_pair(W1, b1, x1, W2, b2, x2, sp, tol, cap)
%REST_PAIR  Develops two networks together from given states until both are at rest.
%
% Implements : simulate_pair's Euler step for both networks of an encounter (saturating contrast sensing), with
%              a constant input added to each; stops when both change by less than tol per step (tol < 0: never
%              stop); used by response_dial, response_dial_as_input, response_dial_from_rest and
%              fig4b_retaliation_from_rest
% Inputs     : W1, b1, x1 and W2, b2, x2 (wiring, constant input and start state of each network), sp
%              (simulate_pair settings), tol, cap (step limit)
% Outputs    : x1, x2 (final states), it (steps taken), ok (true if both came to rest)
% Seed       : none

    km = sp.michaelisConstant; mu = sp.degradationRate; dt = sp.timeStep; ir = sp.sensorIntegrationRate;
    n = sp.hillCoefficient; ok = false;
    for it = 1:cap
        u1 = W1 * x1; u2 = W2 * x2;
        oC = x2(1)^n; oD = x2(2)^n; dO = oC + oD + 1e-12; fC = x1(1)^n; fD = x1(2)^n; dF = fC + fD + 1e-12;
        sF = ir * (oC / dO - oD / dO); sO = ir * (fC / dF - fD / dF);
        u1(3) = u1(3) + sF; u1(4) = u1(4) - sF; u2(3) = u2(3) + sO; u2(4) = u2(4) - sO;
        u1 = max(u1 + b1, 0); u2 = max(u2 + b2, 0);
        n1 = max(x1 + dt * (u1 ./ (km + u1) - mu * x1), 0); n2 = max(x2 + dt * (u2 ./ (km + u2) - mu * x2), 0);
        ch = max(max(abs(n1 - x1)), max(abs(n2 - x2))); x1 = n1; x2 = n2;
        if ch < tol, ok = true; return; end
    end
end
