function [tA, tX] = arm_targets(arm, m, geo0)
%ARM_TARGETS  Target alpha and ||x*|| of an arm of the factor-lowering test (alpha or ||x*|| lowered at the split).
%
% Implements : arm A scales alpha by m and holds ||x*||; the other arm holds alpha and scales ||x*||; used by
%              the three tracking runners of the factor-lowering test (alpha or ||x*|| lowered at the split)
% Inputs     : arm ('A' or 'X'), m (the scale), geo0 (the unchanged network's alpha and xn)
% Outputs    : tA, tX
% Seed       : none

    if arm == 'A', tA = m*geo0.alpha; tX = geo0.xn; else, tA = geo0.alpha; tX = m*geo0.xn; end
end
