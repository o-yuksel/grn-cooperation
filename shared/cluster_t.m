function [e, ci, se] = cluster_t(d, rp)
%CLUSTER_T  Mean of a per-network difference with a 95% t interval, replicates as clusters.
%
% Implements : the cluster interval of the response-dial contrasts (G replicates, G - 1 degrees of freedom);
%              used by dial_loss_shares, response_dial_from_rest, response_dial_recount and fig4_tipping_input
% Inputs     : d (one value per network), rp (its replicate)
% Outputs    : e (mean), ci (95% interval), se
% Seed       : none

    [u, ~, loc] = unique(rp); G = numel(u); N = numel(d);
    Sr = accumarray(loc, d, [G 1]); nr = accumarray(loc, 1, [G 1]); e = sum(Sr) / N;
    se = sqrt(G / (G - 1) * sum((Sr - e * nr) .^ 2)) / N; h = tinv(0.975, G - 1) * se; ci = [e - h, e + h];
end
