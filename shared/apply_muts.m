function [Mm, Bm] = apply_muts(M, B, dr)
%APPLY_MUTS  Applies drawn mutations to one network, one mutant per draw.
%
% Implements : a toggle flips M(g1, g2); a weight change adds 0.2 * z to B(g1, g2); used by run_safety_search,
%              run_safety_search_without_rules and safety_search_losses
% Inputs     : M, B (the network), dr (from draw_muts)
% Outputs    : Mm, Bm (4 x 4 x n mutants)
% Seed       : none

    n = numel(dr.g1); Mm = repmat(M, 1, 1, n); Bm = repmat(B, 1, 1, n);
    for q = 1:n
        if dr.tog(q), Mm(dr.g1(q), dr.g2(q), q) = ~Mm(dr.g1(q), dr.g2(q), q);
        else, Bm(dr.g1(q), dr.g2(q), q) = Bm(dr.g1(q), dr.g2(q), q) + dr.z(q) * 0.2; end
    end
end
