function p = window_pheno(Th, knob, E, P)
%WINDOW_PHENO  The outputs of one genotype over the scoring window.
%
% Implements : on the evolution's developmental map (steady_batch); used by phenotype_shift and w_window
% Inputs     : Th (wiring), knob (the two dial values), E (environments), P (tracking parameters)
% Outputs    : p (2 x environments)
% Seed       : none

    n = size(E.U, 2); p = zeros(2, n); bb = [0; 0; knob(1); knob(1)];
    for g = 1:n
        x = steady_batch(Th, E.U(:, g) + bb, P); p(:, g) = x([3 4]);
    end
end
