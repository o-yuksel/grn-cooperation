function [X, stUsed] = steady_batch(TA, u, P, tol, maxSteps)
%STEADY_BATCH  tracking.steady for a stack of wirings at once.
%
% Implements : the same Euler steps of Eq. 3 and the same stopping rule, applied to each wiring until it
%              settles; the stopping tolerance and step limit default to P.dev_tolerance and P.max_dev_steps;
%              used by the tracking runners and analyses
% Inputs     : TA (n x n x M wirings), u (input), P (tracking parameters), tol and maxSteps (optional)
% Outputs    : X (n x M settled states), stUsed (steps used)
% Seed       : none

    if nargin < 4, tol = P.dev_tolerance; end
    if nargin < 5, maxSteps = P.max_dev_steps; end
    M = size(TA,3); n = size(TA,1); stUsed = 0;
    X = repmat(P.x0, 1, M);
    cols = cell(1,n);
    for j = 1:n, cols{j} = reshape(TA(:,j,:), n, M); end
    live = true(1,M); K = P.K; mu = P.mu; dt = P.dev_dt;
    for st = 1:maxSteps
        Xl = X(:,live);
        B = zeros(n, size(Xl,2));
        for j = 1:n, B = B + cols{j}(:,live) .* Xl(j,:); end
        B = B + u;
        B(B < 0) = 0;
        Xn = Xl + dt*((B ./ (K + B)) - mu*Xl);
        done = sum((Xn - Xl).^2, 1) < tol;
        X(:,live) = Xn;
        idx = find(live); live(idx(done)) = false;
        stUsed = st;
        if ~any(live), break; end
    end
end
