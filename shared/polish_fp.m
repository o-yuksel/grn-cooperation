function [X, conv] = polish_fp(Th, U, P, X)
%POLISH_FP  Newton polish of settled states of a tracking network.
%
% Implements : at most 50 Newton steps per environment on Eq. 3 at rest; conv is true when every residual is
%              below 1e-10; used by the tracking runners and tracking_held_still_checks
% Inputs     : Th (wiring), U (inputs, one column per environment), P (tracking parameters), X (start states)
% Outputs    : X (polished states), conv
% Seed       : none

    n = P.total; conv = true;
    for e = 1:size(U, 2)
        x = X(:, e);
        for it = 1:50
            h = Th*x + U(:, e); hp = max(h, 0);
            F = hp ./ (P.K + hp) - P.mu*x;
            if norm(F) < 1e-13, break; end
            fp = zeros(n, 1); fp(h > 0) = P.K ./ (P.K + h(h > 0)).^2;
            x = x - (diag(fp)*Th - P.mu*eye(n)) \ F;
        end
        h = Th*x + U(:, e); hp = max(h, 0);
        conv = conv && norm(hp ./ (P.K + hp) - P.mu*x) < 1e-10;
        X(:, e) = x;
    end
end
