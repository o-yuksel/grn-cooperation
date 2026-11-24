function [wBar, wBest, ib] = window_eval(M, knob, shift, E, P)
%WINDOW_EVAL  Mean and best window fitness of a population under an arm's dials.
%
% Implements : the original run's best-individual search, under the arm's dials: every individual's fitness
%              averaged over the whole scoring window; used by run_tracking_held_still and
%              run_tracking_moving_target
% Inputs     : M (population), knob, shift, E (environments), P (tracking parameters)
% Outputs    : wBar, wBest, ib (the best individual)
% Seed       : none

    n = size(E.U, 2); np = size(M, 3); d2 = zeros(np, n); bb = [0; 0; knob(1); knob(1)];
    for g = 1:n
        X = steady_batch(M, E.U(:, g) + bb, P);
        d2(:, g) = sum((X([3 4], :) - (E.O(:, g) + shift)).^2, 1)';
    end
    f = exp(-P.selection_strength * d2) * E.wt(:);
    wBar = mean(f); [wBest, ib] = max(f);
end
