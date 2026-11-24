function beta = fragility_fast(Th, u, P, ndraws)
%FRAGILITY_FAST  The mean output shift per mutation of a tracking network, its mutants settled together.
%
% Implements : tracking.fragility's measure (W33 and W44 held), its mutants drawn in the same order and their
%              settled states solved together by steady_batch; used by run_tracking_task and
%              run_tracking_cut_vs_shift
% Inputs     : Th (wiring), u (input), P (tracking parameters), ndraws
% Outputs    : beta (the mean over the draws of the size of the output shift per mutation: a sampled
%              estimate, not the root of summed squares over weights of the README's mutation-sensitivity)
% Seed       : none set here: the mutants are drawn from the caller's stream

    p0 = tracking.steady(Th, u, P); p0 = p0([3 4]);
    n = size(Th,1); TA = zeros(n,n,ndraws);
    for k = 1:ndraws
        Tm = Th; mask = rand(size(Tm)) < P.mutation_prob;
        Tm(mask) = Tm(mask) + randn(nnz(mask),1) * P.mutation_sd;
        Tm(3,3) = Th(3,3); Tm(4,4) = Th(4,4);
        TA(:,:,k) = Tm;
    end
    X = steady_batch(TA, u, P);
    acc = 0;
    for k = 1:ndraws
        acc = acc + norm(X([3 4],k) - p0);
    end
    beta = acc / ndraws;
end
