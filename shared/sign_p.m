function p = sign_p(d)
%SIGN_P  Two-sided exact sign-test p of paired differences.
%
% Implements : binomial with p = 0.5 on the count of positive differences; used by run_tracking_held_still and
%              run_tracking_moving_target
% Inputs     : d (differences)
% Outputs    : p
% Seed       : none

    n = numel(d); k = sum(d > 0);
    p = min(1, 2 * min(binocdf(k, n, 0.5), 1 - binocdf(k - 1, n, 0.5)));
end
