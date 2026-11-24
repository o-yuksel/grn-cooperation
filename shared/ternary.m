function s = ternary(c, a, b)
%TERNARY  Returns a if c is true, else b.
%
% Implements : a choice in printed text and settings; used by the analysis and experiment scripts
% Inputs     : c (logical), a, b
% Outputs    : a or b
% Seed       : none

    if c, s = a; else, s = b; end
end
