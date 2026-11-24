function G = grid_scan(W, sp, deltas)
%GRID_SCAN  The membership test's retaliation condition on its full mutant grid, with the worst-case mutant.
%
% Implements : membership_test's simulate_pair loop for the retaliation condition (C2) with its option gradedC2
%              on (the whole grid is tested, not stopped at the first mutant that earns more than R; full grid,
%              same order, strict >), plus the move of the worst-case mutant, which membership_test does
%              not return; inv: some mutant defects on a cooperating W; wr: W's highest readout against a
%              defecting mutant; used by measure_retaliation_fixed_opponent and safety_search_losses
% Inputs     : W (wiring), sp (simulate_pair settings), deltas (the weight changes)
% Outputs    : G with fields inv, wr, slot, delta, nDef
% Seed       : none

    E = ones(4); G = struct('inv', false, 'wr', NaN, 'slot', NaN, 'delta', NaN, 'nDef', 0);
    best = -inf;
    for ii = 1:4
        for jj = 1:4
            for d = deltas
                if d == 0, continue; end
                Wn = W; Wn(ii, jj) = Wn(ii, jj) + d;
                [oc, f] = simulate_pair(E, W, [], E, Wn, [], sp);
                aR = oc(1) == 'C'; aM = oc(2) == 'C';
                if aR && ~aM, G.inv = true; end
                if ~aM
                    G.nDef = G.nDef + 1; rR = f(1) - f(2);
                    if isfinite(rR) && rR > best
                        best = rR; G.slot = ii + 4 * (jj - 1); G.delta = d;
                    end
                end
            end
        end
    end
    if isfinite(best), G.wr = best; end
end
