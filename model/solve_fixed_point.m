function [x_star, info] = solve_fixed_point(W, s, mu, km, dt, maxSteps, tol)
% SOLVE_FIXED_POINT  Integrates one network with a constant input to its fixed point.
%
% Implements : Eq. 3, dx/dt = f(W x + s) - mu x with f(u) = max(u, 0) / (K_M + max(u, 0)),
%              integrated by the Euler method; used by the 'linear' route of membership_test
% Inputs     : W          n-by-n regulatory matrix (for a pair, the joint 8-by-8 matrix)
%              s          constant input to each gene (n values)
%              mu, km, dt degradation rate mu, Michaelis constant K_M, Euler step Delta t
%              maxSteps   step limit (membership_test: 1000)
%              tol        convergence tolerance (membership_test: 1e-10; a numerical tolerance, not
%                         a model parameter)
% Outputs    : x_star     state at the last step
%              info       .converged true when ||x(t) - x(t-1)||_2 < tol on 10 consecutive steps;
%                         .steps, the step at which that held (maxSteps otherwise)
% Seed       : uses the caller's random stream for the starting state
%
% The start is drawn as in the runs, x_j ~ U(0, 0.01) (Methods, "Model"), so a wiring with more
% than one fixed point can reach different ones from different draws; a caller that needs the
% same result for the same wiring seeds the stream first (membership_test does, per wiring).

    n = length(s);
    x = rand(n, 1) * 0.01;
    consecutive_small = 0;
    info.converged = false;
    info.steps = maxSteps;

    for t = 1:maxSteps
        x_prev = x;
        u = W * x + s;
        u_pos = max(u, 0);
        production = u_pos ./ (km + u_pos);
        x = x + dt * (production - mu * x);     % Euler step of Eq. 3
        x = max(x, 0);                          % expression kept non-negative

        delta = norm(x - x_prev, 2);
        if delta < tol
            consecutive_small = consecutive_small + 1;
            if consecutive_small >= 10
                info.converged = true;
                info.steps = t;
                break;
            end
        else
            consecutive_small = 0;
        end
    end

    x_star = x;
end
