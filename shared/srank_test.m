function p = srank_test(d)
%SRANK_TEST  Wilcoxon signed-rank p of the non-NaN values, NaN where signrank fails.
%
% Implements : used by tracking_load_vs_identity and tracking_responsiveness
% Inputs     : d (differences)
% Outputs    : p
% Seed       : none

    d = d(~isnan(d));
    try, p = signrank(d); catch, p = NaN; end
end
