function S = sample_census_members(nPerArch, seed)
% SAMPLE_CENSUS_MEMBERS  Draws, per signaling architecture, a random sample of the
% per-mutation census's members that have a worst-case mutant H, with their 30 recorded draws.
%
% Implements : the member sample of SI S3.3, "The search" (every member with a worst-case mutant
%              when nPerArch exceeds their number, as in the paper's run)
% Inputs     : nPerArch   members drawn per architecture
%              seed       base seed
%              publication/retaliation_fixed_opponent/retaliation_fixed_opponent_census.mat
%                (written by measure_retaliation_fixed_opponent, run by
%                run_fixed_opponent_census)
% Outputs    : S   struct array, one element per sampled member: architecture, census index,
%                  replicate, clone count, wiring W, and per draw the connection (slot), dW,
%                  toggle flag and the census's outcomes (fields named beside each below)
% Seed       : rng(seed + m, 'twister') for architecture m, one randperm of that architecture's
%              members with a worst-case mutant; the sample is its first min(nPerArch, available)
%              entries, so a smaller sample is a prefix of a larger one under the same seed
%
% A member has a worst-case mutant when its retaliation margin (`cret`) is finite, that is, when at
% least one grid mutant defects on it (SI S2.3, "Members with no defecting mutant").

    D = grn_data_file('publication');
    L = load(fullfile(D, 'retaliation_fixed_opponent', 'retaliation_fixed_opponent_census.mat'), 'out');
    S = [];
    for m = 1:4
        res = L.out.arch{m}.res(:); mem = L.out.arch{m}.memInputs(:);
        idx = find(isfinite([res.cret]'));
        rng(seed + m, 'twister'); p = idx(randperm(numel(idx)));
        p = p(1:min(nPerArch, numel(p)));
        for k = 1:numel(p)
            a = p(k);
            e = struct('arch', m, 'idx', a, 'rep', mem(a).rep, 'count', mem(a).count, ...
                'W', mem(a).W, ...                       % memInputs.W, parent wiring M.*B
                'slot', mem(a).slot(:), 'dW', mem(a).dW(:), 'isToggle', logical(mem(a).isToggle(:)), ...
                'ejCensus', logical(mem(a).ejected(:)), ... % `memInputs.ejected`: ends membership
                'lost', logical(res(a).lost(:)), ...     % res.lost: mutant clone readout <= 0
                'rH', res(a).rH(:), ...                  % res.rH: readout of the mutant against H
                'r', res(a).r, 'wr0', res(a).wr0, 'cret', res(a).cret, ...
                'hSlot', res(a).hSlot, 'hDelta', res(a).hDelta);
            if isempty(S), S = e; else, S(end+1) = e; end %#ok<AGROW>
        end
    end
    S = S(:);
end
