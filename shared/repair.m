function r = repair(M, B, sp)
%REPAIR  The repair of the repair test: cooperators that do not retaliate move to the nearest member genotype.
%
% Implements : classes by the membership test (1 member, 2 cooperates but does not retaliate, 3 other); each
%              individual of class 2 takes the nearest member genotype (Euclidean distance of W, ties to the
%              most common); used by run_repair_test_observed and run_repair_restart_points
% Inputs     : M, B (the population), sp (membership_test parameters)
% Outputs    : r with fields shares (the three classes), nMember, M, B (the repaired population), nReplaced,
%              meanDist
% Seed       : none

    N = size(M, 3);
    G = [reshape(double(M), 16, N); reshape(double(B), 16, N)]';
    [Gu, ~, ig] = unique(G, 'rows'); cnt = accumarray(ig, 1);
    Wg = Gu(:, 1:16) .* Gu(:, 17:32);
    [Wu, ~, iw] = unique(Wg, 'rows');
    cw = zeros(size(Wu, 1), 1);
    for u = 1:size(Wu, 1)
        rr = membership_test(reshape(Wu(u, :), 4, 4), sp);
        cw(u) = 3 - rr.is_weak - rr.is_css;
    end
    cls = cw(iw);
    r.shares = [sum(cnt(cls == 1)), sum(cnt(cls == 2)), sum(cnt(cls == 3))] / N;
    mem = find(cls == 1); r.nMember = numel(mem);
    r.M = M; r.B = B; r.nReplaced = 0; r.meanDist = NaN;
    if isempty(mem), return; end
    cooperate_only = find(cls == 2);
    target = zeros(numel(Gu(:, 1)), 1); dist = nan(size(target));
    for q = cooperate_only'
        dd = sqrt(sum((Wg(mem, :) - Wg(q, :)).^2, 2));
        best = find(dd == min(dd));
        [~, b] = max(cnt(mem(best)));
        target(q) = mem(best(b)); dist(q) = min(dd);
    end
    ind = find(cls(ig) == 2);
    for i = ind'
        t = target(ig(i));
        r.M(:,:,i) = reshape(logical(Gu(t, 1:16)), 4, 4);
        r.B(:,:,i) = reshape(single(Gu(t, 17:32)), 4, 4);
    end
    r.nReplaced = numel(ind);
    r.meanDist = mean(dist(ig(ind)));
end
