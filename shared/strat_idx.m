function idx = strat_idx(arch)
%STRAT_IDX  Replicates resampled with replacement within each architecture.
%
% Implements : the stratified bootstrap draw of the clone-weighted fits and shares; used by
%              conditional_cooperation_measures (boot_fit) and cooperate_only_losses
% Inputs     : arch (architecture of each replicate)
% Outputs    : idx (resampled replicate indexes, architecture by architecture)
% Seed       : none set here: the draws come from the caller's stream

    idx = [];
    for m = unique(arch)'
        pool = find(arch == m);
        idx = [idx; pool(randi(numel(pool), numel(pool), 1))]; %#ok<AGROW>
    end
end
