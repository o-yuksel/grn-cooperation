function snaps = extract_snapshots(results, reps, cfg)
%EXTRACT_SNAPSHOTS  The census sample of stored populations: snapshots reduced to distinct wirings.
%
% Implements : cfg.nSnaps stored generations at or after cfg.steadyStateFrom per replicate, spaced evenly by
%              index, each reduced to its distinct wirings, their clone counts and one carrier's M and B; used
%              by per_mutation_census, cooperate_only_census and retaliation_margin_vs_rho
% Inputs     : results (a run's results cell), reps (the replicates, a vector), cfg (steadyStateFrom, nSnaps)
% Outputs    : snaps (one cell per replicate in reps: struct rep, S with fields gen, Wu, counts, N, Mfirst,
%              Bfirst)
% Seed       : none

    snaps = cell(numel(reps),1);
    for ii = 1:numel(reps)
        r = reps(ii);
        rep = results{r};
        sg = double(rep.historyGenerations(:))';
        ks = find(sg >= cfg.steadyStateFrom);
        ks = ks(round(linspace(1, numel(ks), cfg.nSnaps)));
        S = struct('gen', {}, 'Wu', {}, 'counts', {}, 'N', {}, 'Mfirst', {}, 'Bfirst', {});
        for s = 1:numel(ks)
            k = ks(s);
            M = double(rep.historyM(:,:,:,k)); B = double(rep.historyB(:,:,:,k));
            N = size(M,3); nG = size(M,1);
            Wall = reshape(M.*B, nG*nG, N)';
            [Wu, ~, ic] = unique(Wall, 'rows');
            counts = accumarray(ic, 1);
            firstIdx = accumarray(ic, (1:N)', [], @min);
            S(s).gen    = sg(k);
            S(s).Wu     = Wu;
            S(s).counts = counts;
            S(s).N      = N;
            S(s).Mfirst = M(:,:,firstIdx);
            S(s).Bfirst = B(:,:,firstIdx);
        end
        snaps{ii} = struct('rep', r, 'S', S);
    end
end
