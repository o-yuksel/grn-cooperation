function [x1s, x2s] = simpair_batch(W1, W2)
%SIMPAIR_BATCH  simulate_pair for many pairs of saturating-contrast wirings at once.
%
% Implements : the membership test's encounter rule: start 0.005, 600 steps, a convergence check every 40 steps,
%              each network read where it first settles or where both do; used by member_batch's callers
%              (run_safety_search, run_safety_search_without_rules, safety_search_losses)
% Inputs     : W1, W2 (4 x 4 x n wirings, first and second network of each pair)
% Outputs    : x1s, x2s (4 x n settled states)
% Seed       : none

    n = size(W1, 3); x1 = 0.005 * ones(4, n); x2 = x1; x1o = x1; x2o = x2; s1 = x1; s2 = x2;
    st1 = false(1, n); st2 = false(1, n); act = true(1, n); thr = 0.01^2; mu = 0.1; last = 0;
    for t = 1:600
        a = find(act); if isempty(a), break; end
        y1 = x1(:, a); y2 = x2(:, a);
        in1 = reshape(pagemtimes(W1(:, :, a), reshape(y1, 4, 1, [])), 4, []);
        in2 = reshape(pagemtimes(W2(:, :, a), reshape(y2, 4, 1, [])), 4, []);
        d2 = y2(1, :) + y2(2, :) + 1e-12; d1 = y1(1, :) + y1(2, :) + 1e-12;
        sF = y2(1, :) ./ d2 - y2(2, :) ./ d2; sO = y1(1, :) ./ d1 - y1(2, :) ./ d1;
        in1(3, :) = in1(3, :) + sF; in1(4, :) = in1(4, :) - sF;
        in2(3, :) = in2(3, :) + sO; in2(4, :) = in2(4, :) - sO;
        in1 = max(in1, 0); in2 = max(in2, 0);
        y1 = max(y1 + (in1 ./ (1 + in1) - mu * y1), 0); y2 = max(y2 + (in2 ./ (1 + in2) - mu * y2), 0);
        x1(:, a) = y1; x2(:, a) = y2;
        if t - last >= 40
            m1 = max(max(y1, [], 1), 1e-9); m2 = max(max(y2, [], 1), 1e-9);
            c1 = sum((y1 - x1o(:, a)).^2, 1) ./ m1.^2 < thr; c2 = sum((y2 - x2o(:, a)).^2, 1) ./ m2.^2 < thr;
            both = c1 & c2; ab = a(both);
            s1(:, ab) = y1(:, both); s2(:, ab) = y2(:, both); st1(ab) = true; st2(ab) = true; act(ab) = false;
            n1 = a(~both & c1 & ~st1(a)); s1(:, n1) = x1(:, n1); st1(n1) = true;
            n2 = a(~both & c2 & ~st2(a)); s2(:, n2) = x2(:, n2); st2(n2) = true;
            x1o(:, a) = x1(:, a); x2o(:, a) = x2(:, a); last = t;
        end
    end
    s1(:, ~st1) = x1(:, ~st1); s2(:, ~st2) = x2(:, ~st2); x1s = s1; x2s = s2;
end
