function d = best_check(M, ThBest, E, P)
%BEST_CHECK  Whether the best-individual search finds the saved best genotype in a rebuilt population.
%
% Implements : the original run's best-individual search on the rebuilt population; used by
%              run_tracking_held_still and run_tracking_moving_target
% Inputs     : M (rebuilt population), ThBest (saved best), E (environments), P (tracking parameters)
% Outputs    : d (largest weight difference; 0 when the search finds ThBest)
% Seed       : none

    n = size(E.U, 2); np = size(M, 3); d2 = zeros(np, n);
    for g = 1:n
        X = steady_batch(M, E.U(:, g), P);
        d2(:, g) = sum((X([3 4], :) - E.O(:, g)).^2, 1)';
    end
    fb = mean(exp(-P.selection_strength * d2), 2);
    [~, ib] = max(fb);
    d = max(abs(reshape(M(:, :, ib) - ThBest, [], 1)));
end
