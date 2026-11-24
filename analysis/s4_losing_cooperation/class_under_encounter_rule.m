function r = class_under_encounter_rule(W, sp, seed, needMember)
% CLASS_UNDER_ENCOUNTER_RULE  The class of one wiring (member, cooperator that does not retaliate,
% or neither) when the cooperation (C1) and retaliation (C2) conditions are read under the runs'
% encounter rule instead of the membership test.
%
% Implements : the rule of SI S4.5 as membership_test_vs_encounter_rule applies it (Methods,
%              "Failure route"); called by repair_vs_sham_checks.
%              Encounters: K = 20 start pairs x ~ U(0, 0.01) per gene, the same in every
%              pairing; up to 300 steps, settled when checked every 40 steps; a network
%              cooperates if x_C > x_D.
%              C1: the clone pairing cooperates in at least half of its 2K actions.
%              C2: from the single-weight mutants (each weight shifted by one of the membership
%              test's 10 steps in [-0.6, 0.6]) that defect on the wiring from the test's start
%              (0.005, 600 steps) or in either of 2 screening encounters, up to 20 are sampled
%              and play K encounters each. C2 fails if one that defects in at least half of them
%              meets defection back in under half of its defections, or if one's mean payoff
%              exceeds the clone's by a one-sided paired z above 1.645.
% Inputs     : W (4 x 4 wiring), sp (development parameters of the architecture), seed,
%              needMember (false: only C1 is read)
% Outputs    : r.cls = 1 member, 2 cooperates but does not retaliate, 3 fails C1 (NaN when it
%              passes C1 and needMember is false); r.pC1, the C1 share; r.pool, the number of
%              mutants in the pool; nothing printed or saved
% Seed       : rng(seed, 'twister'), the caller's seed: the start states and the mutant sample

    K = 20; nNeigh = 20; nScreen = 2;
    rng(seed, 'twister');
    X1 = rand(4, K) * 0.01; X2 = rand(4, K) * 0.01;
    [c3a, c3b] = encounter_batch(repmat(W, 1, 1, K), repmat(W, 1, 1, K), X1, X2, sp, 300, 'encounter_rule');
    r.pC1 = mean([c3a c3b]); r.cls = 3; r.pool = 0;
    if r.pC1 < 0.5, return; end
    r.cls = NaN;
    if ~needMember, return; end
    D = grid_deltas();
    X0 = 0.005 * ones(4, 1);
    F = grid_pairs(W, sp, 600, X0, X0, D, 'simpair');
    S = grid_pairs(W, sp, 300, X1(:, 1:nScreen), X2(:, 1:nScreen), D, 'encounter_rule');
    pool = find(~F.aM(:)' | S.pDM > 0);
    r.pool = numel(pool);
    nbr = sort(pool(randperm(numel(pool), min(nNeigh, numel(pool)))));
    member = true;
    if ~isempty(nbr)
        g = grid_pairs(W, sp, 300, X1, X2, D(nbr, :), 'encounter_rule');
        def = g.pDM >= 0.5;
        db = double(g.dbar); sd = double(g.sd);
        ex = db > 0 & db > 1.6449 * sd / sqrt(K);
        member = all(g.retal(def) >= 0.5) && ~any(ex);
    end
    r.cls = 2 - member;
end


% =========================================================================
% The membership test's mutant grid: every entry (i, j) of W shifted by one of 10 steps.
function D = grid_deltas()
    d = linspace(-0.6, 0.6, 11); d = d(d ~= 0);
    [dd, jj, ii] = ndgrid(d, 1:4, 1:4);
    D = [ii(:), jj(:), dd(:)];
end


function g = grid_pairs(W, sp, steps, X1, X2, D, snapMode)
    K = size(X1, 2); nN = size(D, 1);
    Wn = repmat(W, 1, 1, nN);
    for q = 1:nN, Wn(D(q,1), D(q,2), q) = Wn(D(q,1), D(q,2), q) + D(q,3); end
    W2 = cat(3, repmat(W, 1, 1, K), Wn(:,:,repelem(1:nN, K)));
    [a1, a2] = encounter_batch(repmat(W, 1, 1, K * (nN + 1)), W2, repmat(X1, 1, nN + 1), ...
        repmat(X2, 1, nN + 1), sp, steps, snapMode);
    a1 = reshape(a1, K, nN + 1); a2 = reshape(a2, K, nN + 1);
    aR = a1(:, 2:end); aM = a2(:, 2:end);
    % d: the mutant's payoff minus the clone's, per draw (T, R, P, S = 5, 3, 1, 0)
    pay = @(r, m) 1 * (~r & ~m) + 0 * (~r & m) + 5 * (r & ~m) + 3 * (r & m);
    d = pay(aR, aM) - pay(a1(:,1), a2(:,1));
    nD = sum(~aM, 1);
    g.aM = aM;
    g.pDM = single(mean(~aM, 1));
    g.retal = single(sum(~aR & ~aM, 1) ./ max(nD, 1)); g.retal(nD == 0) = NaN;
    g.dbar = single(mean(d, 1)); g.sd = single(std(d, 0, 1));
end
