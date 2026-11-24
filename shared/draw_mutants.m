function TA = draw_mutants(Th, P, M)
%DRAW_MUTANTS  M draws of the tracking task's mutation distribution.
%
% Implements : made as tracking.fragility makes them (W33 and W44 held), so the same seed gives the same draws;
%              used by tracking_load_excess, tracking_load_vs_identity and run_tracking_single_genotypes
% Inputs     : Th (wiring), P (tracking parameters), M (number of mutants)
% Outputs    : TA (n x n x M mutants)
% Seed       : none set here: the draws come from the caller's stream

    n = size(Th,1); TA = zeros(n,n,M);
    for k = 1:M
        Tm = Th; mask = rand(size(Tm)) < P.mutation_prob;
        Tm(mask) = Tm(mask) + randn(nnz(mask),1) * P.mutation_sd;
        Tm(3,3) = Th(3,3); Tm(4,4) = Th(4,4);
        TA(:,:,k) = Tm;
    end
end
