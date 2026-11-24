function [tgt, oth] = target_other(ra, rx, k)
%TARGET_OTHER  The targeted and the held quantity of an arm.
%
% Implements : arm 1 targets alpha and holds ||x*||, arm 2 the reverse; used by rep_pass's callers
% Inputs     : ra, rx (alpha and ||x*|| ratios), k (1 or 2)
% Outputs    : tgt, oth
% Seed       : none

    if k == 1, tgt = ra; oth = rx; else, tgt = rx; oth = ra; end
end
