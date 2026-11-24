function [bA, bX, kw] = best_series(J, rs, a)
%BEST_SERIES  The window-best genotype's alpha and ||x*|| over the generations after the split.
%
% Implements : relative to the CONTROL split ThBest; used by report_added_input
% Inputs     : J (the arms' runs), rs (replicates), a (arm)
% Outputs    : bA, bX (replicate x generation), kw (generations)
% Seed       : none

    c0 = J([J.r] == rs(1) & [J.arm] == 1); ok = isfinite([c0.track.wBest]);
    kw = [c0.track(ok).k]; bA = nan(numel(rs), numel(kw)); bX = bA;
    for q = 1:numel(rs)
        c = J([J.r] == rs(q) & [J.arm] == 1); z = J([J.r] == rs(q) & [J.arm] == a);
        if isempty(z) || z.void, continue; end
        t = z.track(isfinite([z.track.wBest]));
        bA(q, :) = [t.alphaBest] / c.track(1).alphaBest;
        bX(q, :) = [t.xnBest] / c.track(1).xnBest;
    end
end
