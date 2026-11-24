function c = coop_now(M, B, config, n)
%COOP_NOW  Cooperation of a population under the runs' encounter rule.
%
% Implements : the mean of n one-generation evaluations of evolve_population (its call that returns the
%              individuals' cooperation without evolving the population); used by
%              run_repair_test, run_repair_test_observed and run_repair_restart_points
% Inputs     : M, B (the population), config (the run's config), n (evaluations)
% Outputs    : c (mean cooperation)
% Seed       : none set here: the evaluations draw from the caller's stream

    K = ones(4, size(M, 3), 'single') * config.michaelisConstant;
    c = 0;
    for q = 1:n
        x = evolve_population('__test_evaluateFitness__', M, B, K, config);
        c = c + x{2} / n;
    end
end
