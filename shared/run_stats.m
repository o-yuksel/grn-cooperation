function S = run_stats(o)
%RUN_STATS  The cut in mutation load of one tracking run, with its interval.
%
% Implements : the cut is minus the mean paired difference from CONTROL over CONTROL's median mutation load, in
%              percent, with the paired bootstrap interval; the ALPHA/XNORM ratio of the two differences; the
%              second control's difference; used by run_tracking_cut_vs_shift and
%              run_tracking_cut_vs_mutation_size
% Inputs     : o (a run's saved output, with o.analysis.table and o.analysis.Y)
% Outputs    : S (c0, cutA, ciA, nA, cutX, ciX, nX, ratio, ratioCi, d2, ci2, pct2)
% Seed       : none

    T = o.analysis.table; Y = o.analysis.Y; c0 = median(Y(:, 1, 1), 'omitnan');
    S.c0 = c0;
    S.cutA = -100 * T(1, 3).d / c0; S.ciA = -100 * fliplr(T(1, 3).ci) / c0; S.nA = T(1, 3).n;
    S.cutX = -100 * T(1, 4).d / c0; S.ciX = -100 * fliplr(T(1, 4).ci) / c0; S.nX = T(1, 4).n;
    rb = T(1, 3).boot ./ T(1, 4).boot;
    S.ratio = T(1, 3).d / T(1, 4).d; S.ratioCi = prctile(rb, [2.5 97.5]);
    S.d2 = T(1, 2).d; S.ci2 = T(1, 2).ci; S.pct2 = 100 * T(1, 2).d / c0;
end
