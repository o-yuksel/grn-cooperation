function X = integrate_fp(Th, U, P)
%INTEGRATE_FP  Settled states of one tracking network in several environments, by Euler steps.
%
% Implements : Eq. 3 integrated from x0 until the squared step is below 1e-12 (at most 200,000 steps), all
%              environments together; used by the tracking runners and tracking_held_still_checks
% Inputs     : Th (wiring), U (one input column per environment), P (tracking parameters)
% Outputs    : X (one settled state per column)
% Seed       : none

    M = size(U, 2); X = repmat(P.x0, 1, M); live = true(1, M);
    for s = 1:200000
        B = Th*X(:, live) + U(:, live); B(B < 0) = 0;
        dX = B ./ (P.K + B) - P.mu*X(:, live);
        X(:, live) = X(:, live) + P.dev_dt*dX;
        done = sum(dX.^2, 1) < 1e-12; idx = find(live); live(idx(done)) = false;
        if ~any(live), break; end
    end
end
