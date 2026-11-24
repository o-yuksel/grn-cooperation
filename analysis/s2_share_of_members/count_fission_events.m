function n = count_fission_events(historyParents, Ng)
% COUNT_FISSION_EVENTS  In how many generations of a group-selection run did exactly one group
% split and replace another? A check of the one split (fission) per generation of the design.
%
% Implements : a check of SI Extended Materials and Methods, "Group selection"; called by
%              group_selection_tests
% Inputs     : historyParents   N x nGens parent record of one population run by evolve_population
%                               under group selection ('fission') or the random-split control
%                               ('fission_sham'); entry (i, t) is the parent of position i
%              Ng               group size
% Outputs    : n   generations with exactly one replaced group; nGens when every generation has one
% Seed       : none (no random draws)
%
% Within a group every parent comes from the group's own positions, so a group whose parents all
% come from other groups' positions is the one a split replaced.

    [N, nGens] = size(historyParents);
    G = N / Ng;
    n = 0;
    for t = 1:nGens
        src = reshape(ceil(double(historyParents(:, t)) / Ng), Ng, G);
        n = n + (sum(all(src ~= (1:G), 1)) == 1);
    end
end
