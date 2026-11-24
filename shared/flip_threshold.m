function thr = flip_threshold(flips)
%FLIP_THRESHOLD  The smallest input size at which a flip test turns true (the tipping input).
%
% Implements : a scan of 41 sizes from 1e-4 to 10 (log spaced), then bisection to 0.1% of the size; Inf if no
%              size flips; used by safety_search_losses, tipping_input_flag and worked_example
% Inputs     : flips (handle: input size -> true if the network flips)
% Outputs    : thr
% Seed       : none

    g = logspace(-4, 1, 41); lo = 0; hi = NaN;
    for k = 1:numel(g)
        if flips(g(k)), hi = g(k); break; end
        lo = g(k);
    end
    if isnan(hi), thr = Inf; return; end
    for it = 1:60
        if hi - lo <= 1e-3 * hi, break; end
        mid = (lo + hi) / 2;
        if flips(mid), hi = mid; else, lo = mid; end
    end
    thr = hi;
end
