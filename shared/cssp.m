function P = cssp()
%CSSP  The membership-test settings of the saturating-contrast architecture.
%
% Implements : contrast sensing, Hill 1, K_M 1, mu 0.1, dt 1, 600 steps, the simulate_pair route; used by
%              run_safety_search and run_safety_search_without_rules
% Inputs     : none
% Outputs    : P (membership_test parameters)
% Seed       : none

    P = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'backend', 'simulate_pair');
end
