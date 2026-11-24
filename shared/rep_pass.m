function ok = rep_pass(ra, rx, k, cfg)
%REP_PASS  Whether a split hits its target in the factor-lowering test (alpha or ||x*|| lowered at the split).
%
% Implements : the targeted quantity within targetTol of m and the other within otherTol of 1; used by
%              run_tracking_held_still and run_tracking_moving_target
% Inputs     : ra, rx (alpha and ||x*|| ratios), k (1 for the alpha arm), cfg (m, targetTol, otherTol)
% Outputs    : ok
% Seed       : none

    [tgt, oth] = target_other(ra, rx, k);
    ok = abs(tgt/cfg.m - 1) <= cfg.targetTol && abs(oth - 1) < cfg.otherTol;
end
