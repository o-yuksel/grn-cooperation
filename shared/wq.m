function q = wq(v, w, p)
%WQ  Weighted quantile over the finite values with positive weight.
%
% Implements : linear interpolation of the weighted midpoints, NaN when nothing is left; used by
%              retaliation_margin_vs_rho and two_margin_plane
% Inputs     : v (values), w (weights), p (probability)
% Outputs    : q
% Seed       : none

    v = v(:); w = double(w(:));
    ok = isfinite(v) & isfinite(w) & w > 0;
    v = v(ok); w = w(ok);
    if isempty(v), q = NaN; return; end
    [v, i] = sort(v); w = w(i);
    cw = cumsum(w) - 0.5*w; cw = cw / sum(w);
    if p <= cw(1), q = v(1); return; end
    if p >= cw(end), q = v(end); return; end
    q = interp1(cw, v, p, 'linear');
end
