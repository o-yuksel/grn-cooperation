function p = signflip(K, d, rp)
%SIGNFLIP  Two-sided sign-flip p over replicates for a per-network difference.
%
% Implements : 10,000 random sign vectors on the replicate sums; p = (1 + #{|flipped sum| >= |observed|}) /
%              10,001; used by dial_loss_shares and the response-dial analyses
% Inputs     : K (K.seed), d (one value per network), rp (its replicate)
% Outputs    : p
% Seed       : rng(K.seed + 54000, twister) for every contrast; the caller's stream is restored

    [~, ~, loc] = unique(rp); Sr = accumarray(loc, d); T = abs(sum(Sr));
    rs = rng; rng(K.seed + 54000, 'twister'); s = 2 * (rand(numel(Sr), 10000) < 0.5) - 1; rng(rs);
    p = (1 + sum(abs(Sr' * s) >= T - 1e-12)) / (1 + 10000);
end
