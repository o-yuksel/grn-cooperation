function [M, nd, top] = rebuild_population(popTh, popFreq, N, P)
%REBUILD_POPULATION  Rebuilds a population from its distinct wirings and their frequencies.
%
% Implements : distinct wirings times clone counts, in first-occurrence order; used by run_tracking_held_still
%              and run_tracking_moving_target
% Inputs     : popTh (one wiring per row), popFreq (frequencies), N (population size), P (tracking parameters)
% Outputs    : M (n x n x N), nd (distinct wirings), top (largest frequency)
% Seed       : none

    cnt = round(popFreq * N); n = P.total;
    assert(sum(cnt) == N, 'clone counts do not sum to N');
    M = zeros(n, n, N); k = 0;
    for i = 1:size(popTh, 1)
        T = reshape(popTh(i, :), n, n);
        M(:, :, k+1:k+cnt(i)) = repmat(T, 1, 1, cnt(i)); k = k + cnt(i);
    end
    nd = size(popTh, 1); top = max(popFreq);
end
