function out = alpha_and_beta(W, x_star, s, mu, km, grad_r, kink_tol)
% ALPHA_AND_BETA  Responsiveness alpha and mutation-sensitivity beta of the readout at a fixed
% point x*, through the response operator Phi shared by the partner's signal and the wiring.
%
% Implements : SI S1.2 (the response operator Phi) and S1.3 (alpha, beta, Theorem S1.3);
%              Methods, "Congruence identity": u = W x* + s; D_f = diag(f'(u)) with
%              f'(u) = K_M / (K_M + u)^2 where u > kink_tol and 0 where u < -kink_tol;
%              A = mu I - D_f W; Phi = A^-1 D_f; v = Phi' grad_r; alpha = ||v||_2;
%              beta = alpha ||x*||_2
% Inputs     : W 4-by-4 wiring; x_star the fixed point; s the signal held fixed at each gene's
%              input; mu the degradation rate; km K_M; grad_r the readout gradient ([1; -1; 0; 0]
%              for R = x_C - x_D); kink_tol (optional, 1e-6) the band around u = 0 in which a
%              gene counts as at the rectifier's corner
% Outputs    : out with u, active / inactive / near_kink (logical per gene), and, when valid,
%              D_f, Phi, v, alpha, beta; valid is false (alpha, beta NaN) when any gene is at
%              the corner or |det A| < 1e-12. Stability of x* is not tested here.
% Seed       : none (no random draws)

    % kink_tol is a numerical tolerance, not a value the paper states
    if nargin < 7 || isempty(kink_tol), kink_tol = 1e-6; end

    n = length(x_star);
    u = W * x_star + s;

    out.active    = (u > kink_tol);
    out.inactive  = (u < -kink_tol);
    out.near_kink = (abs(u) <= kink_tol);

    out.u = u;
    out.valid = false;
    out.Phi = NaN(n);
    out.v = NaN(n, 1);
    out.alpha = NaN;
    out.beta = NaN;

    if any(out.near_kink)
        return;
    end

    fprime = zeros(n, 1);
    fprime(out.active) = km ./ (km + u(out.active)).^2;
    out.D_f = diag(fprime);

    A = mu * eye(n) - out.D_f * W;
    if abs(det(A)) < 1e-12
        return;
    end

    out.Phi = A \ out.D_f;
    out.v = (grad_r(:).' * out.Phi).';
    out.alpha = norm(out.v, 2);
    out.beta  = out.alpha * norm(x_star, 2);
    out.valid = true;
end
