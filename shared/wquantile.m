function q = wquantile(v, w, p)
%WQUANTILE  Weighted quantile by linear interpolation of the weighted midpoints.
%
% Implements : the clone-weighted quantiles of the census tables; used by per_mutation_census and
%              cooperate_only_losses
% Inputs     : v (values), w (weights), p (probability)
% Outputs    : q
% Seed       : none

    [v, i] = sort(v(:)); w = w(i); w = w(:);
    cw = cumsum(w) - 0.5*w;
    cw = cw / sum(w);
    if p <= cw(1), q = v(1); return; end
    if p >= cw(end), q = v(end); return; end
    q = interp1(cw, v, p, 'linear');
end
