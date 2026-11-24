function [mem, info] = member_batch(Wm, full)
%MEMBER_BATCH  The membership test for many saturating-contrast wirings at once.
%
% Implements : cooperation with the clone, then (if full) every single-weight mutant of the grid membership_test
%              uses (11 values from -0.6 to 0.6 without 0) against the clone; used by run_safety_search,
%              run_safety_search_without_rules and safety_search_losses
% Inputs     : Wm (4 x 4 x n wirings), full (false: the clone test only)
% Outputs    : mem (n x 1, true for a member), info (per wiring: weak, true when it cooperates with its clone;
%              rClone, its readout
%              against the clone; worstRet, the highest readout of the wiring against a mutant that defects on it;
%              nDefMut, the number of such mutants; minRM, the lowest readout among its mutants)
% Seed       : none

    n = size(Wm, 3); grid = linspace(-3 * 0.2, 3 * 0.2, 11); grid = grid(grid ~= 0);
    [xa, ~] = simpair_batch(Wm, Wm); rC = xa(1, :) - xa(2, :);
    info.weak = (xa(1, :) > xa(2, :))'; info.rClone = rC';
    mem = false(n, 1); info.worstRet = nan(n, 1); info.nDefMut = nan(n, 1); info.minRM = nan(n, 1);
    if ~full, return; end
    w = find(info.weak); nw = numel(w); if nw == 0, return; end
    nm = 16 * numel(grid); Wr = zeros(4, 4, nw * nm); Wx = Wr; c = 0;
    for a = 1:nw
        W = Wm(:, :, w(a));
        for slot = 1:16
            for dl = grid
                c = c + 1; Wr(:, :, c) = W; Wq = W; Wq(slot) = Wq(slot) + dl; Wx(:, :, c) = Wq;
            end
        end
    end
    [xr, xm] = simpair_batch(Wr, Wx);
    resC = reshape(xr(1, :) > xr(2, :), nm, nw); mutC = reshape(xm(1, :) > xm(2, :), nm, nw);
    rR = reshape(xr(1, :) - xr(2, :), nm, nw); rM = reshape(xm(1, :) - xm(2, :), nm, nw);
    inv = resC & ~mutC; mem(w) = ~any(inv, 1)';
    for a = 1:nw
        dm = ~mutC(:, a); info.nDefMut(w(a)) = sum(dm);
        if any(dm), info.worstRet(w(a)) = max(rR(dm, a)); end
        info.minRM(w(a)) = min(rM(:, a));
    end
end
