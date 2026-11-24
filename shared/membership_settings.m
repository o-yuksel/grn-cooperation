function params = membership_settings(cfg)
%MEMBERSHIP_SETTINGS  Membership-test settings from a run's config.
%
% Implements : 600 steps, the simulate_pair route, and membership_test's defaults written out (sigma_mut = 0.2,
%              the 11 grid values from -0.6 to 0.6); used by scope_checks_gen500 and silent_and_moving_members
% Inputs     : cfg (a run's config)
% Outputs    : params (membership_test parameters)
% Seed       : none

    params = struct('sensingMode', cfg.sensingMode, 'hillCoefficient', cfg.hillCoefficient, ...
        'discreteOutput', cfg.discreteOutput, 'michaelisConstant', cfg.michaelisConstant, ...
        'degradationRate', cfg.degradationRate, 'timeStep', cfg.timeStep, ...
        'developmentSteps', 600, 'backend', 'simulate_pair', 'sigma_mut', 0.2, ...
        'cssMagGrid', linspace(-3*0.2, 3*0.2, 11));
end
