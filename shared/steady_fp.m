function x = steady_fp(Th, u, P)
%STEADY_FP  A tracking network's settled state by tracking.steady's Euler steps, to a tight stop.
%
% Implements : stopped when the squared residual of Eq. 3 is below 1e-18 (at most 50,000 steps); used by the
%              tracking analyses and run_tracking_single_genotypes
% Inputs     : Th (wiring), u (input), P (tracking parameters)
% Outputs    : x (settled state)
% Seed       : none

    x = P.x0; K = P.K; mu = P.mu; dt = P.dev_dt;
    for s = 1:50000
        b = Th*x + u; b(b < 0) = 0;
        dx = (b ./ (K + b)) - mu*x;
        x = x + dt*dx;
        if sum(dx.^2) < 1e-18, break; end
    end
end
