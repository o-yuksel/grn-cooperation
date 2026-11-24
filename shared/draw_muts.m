function dr = draw_muts(rs, n)
%DRAW_MUTS  Draws n mutations from the runs' mutation operator.
%
% Implements : toggle with probability 0.2, else a weight change; the slot (g1, g2) uniform; z standard normal
%              (scaled by 0.2 in apply_muts); used by run_safety_search, run_safety_search_without_rules and
%              safety_search_losses
% Inputs     : rs (a RandStream), n
% Outputs    : dr with fields tog, g1, g2, z
% Seed       : the caller's stream rs

    dr.tog = rand(rs, n, 1) < 0.2; dr.g1 = randi(rs, 4, n, 1); dr.g2 = randi(rs, 4, n, 1); dr.z = randn(rs, n, 1);
end
