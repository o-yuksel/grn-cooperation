function [ok, cls] = toggle_stable(M, B, P, toggleSteps)
% TOGGLE_STABLE  Does a member also resist its single topology toggles? True when no mutant that
% switches one connection on or off, M(i,j) -> ~M(i,j), earns more against the wiring than the
% wiring's own clone does.
%
% Implements : the membership test extended to topology toggles (SI S2.3, "Scope: mutations and
%              wirings outside the test"), for scope_checks_gen1 and scope_checks_gen500
% Inputs     : M, B          4-by-4 connections (1 = present) and signed weights of the wiring
%                            (W = M .* B), kept apart so that a removed connection keeps its weight
%              P             membership-test settings (sensingMode, hillCoefficient, discreteOutput,
%                            michaelisConstant, degradationRate, timeStep)
%              toggleSteps   development step limit of each pair (the callers pass 600, as in the
%                            membership test)
% Outputs    : ok    true when none of the 16 toggle mutants invades
%              cls   counts of invading toggles: [among the four sensor-to-output connections
%                    M(1:2,3:4), among the other twelve]
% Seed       : none (no random draws)
%
% A toggle mutant invades when the wiring cooperates against it and it defects, so it earns T = 5
% against the clone's R = 3: the retaliation condition (C2) at the paper's payoffs. Each pair
% develops with simulate_pair from the fixed start, with the membership test's pair settings.
    Mask = ones(4); Kv = [];
    sim_params = struct('sensingMode', P.sensingMode, ...
        'hillCoefficient', P.hillCoefficient, 'discreteOutput', P.discreteOutput, ...
        'michaelisConstant', P.michaelisConstant, 'degradationRate', P.degradationRate, ...
        'timeStep', P.timeStep, 'developmentSteps', toggleSteps, ...
        'payoffDC', 0, 'externalSignal1', zeros(4,1), 'externalSignal2', zeros(4,1));
    W = M .* B;
    isPartner = false(4); isPartner(1:2,3:4) = true;
    ok = true; cls = [0 0];
    for i = 1:4
        for j = 1:4
            Mm = M; Mm(i,j) = 1 - Mm(i,j);
            Wm = Mm .* B;
            out = simulate_pair(Mask, W, Kv, Mask, Wm, Kv, sim_params);
            % Invades when the wiring cooperates and the mutant defects (it collects T > R)
            if out(1) == 'C' && out(2) ~= 'C'
                ok = false;
                if isPartner(i,j), cls(1) = cls(1)+1; else, cls(2) = cls(2)+1; end
            end
        end
    end
end
