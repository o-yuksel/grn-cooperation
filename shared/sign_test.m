function p = sign_test(d)
%SIGN_TEST  Two-sided exact sign-test p of the nonzero, non-NaN differences.
%
% Implements : binomial with p = 0.5 on the count of positive differences; p = 1 when none is left; used by
%              tracking_load_vs_identity and tracking_responsiveness
% Inputs     : d (differences)
% Outputs    : p
% Seed       : none

    d = d(d ~= 0 & ~isnan(d)); n = numel(d);
    if n == 0, p = 1; return; end
    k = sum(d > 0);
    p = min(1, 2 * min(binocdf(k, n, 0.5), 1 - binocdf(k-1, n, 0.5)));
end
