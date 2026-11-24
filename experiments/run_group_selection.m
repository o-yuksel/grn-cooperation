function out = run_group_selection(arm, phase, nReps, nGens, nWorkers)
% RUN_GROUP_SELECTION  Evolves saturating-contrast populations of N = 1,000 in one arm: in 20
% groups of 50, with one group splitting per generation by payoff (group selection) or at random
% (the random-split control), or undivided.
%
% Produces : the runs behind SI S2.7, Table S2 and Fig. 3 (group selection); Methods, "Game and
%            evolution" (read by group_selection_arms)
% Inputs   : arm       'fission' (group selection), 'random_split' (random-split control) or
%                      'control' (undivided)
%            phase     'main' (the paper's runs; the default and only value)
%            nReps, nGens   defaults 185, 500
%            nWorkers  pool size if no pool is open (default 2)
% Outputs  : publication/group_two_level/group_fission_<arm>_N1000_<nReps>reps_<nGens>gen.mat
%              (an existing file stops the run)
% Run      : run_group_selection(arm, 'main', 185, 500) for each arm (run_group_selection_arms
%            runs all three; about 20 s per population per worker, 2 workers)
% Seed     : replicate r runs under rng(seedBase + r, 'twister'); seedBase 770000 (group
%            selection), 780000 (random-split control), 790000 (undivided)

    if nargin < 2 || isempty(phase), phase = 'main'; end
    if nargin < 3 || isempty(nReps), nReps = 185; end
    if nargin < 4 || isempty(nGens), nGens = 500; end
    if nargin < 5 || isempty(nWorkers), nWorkers = 2; end

    outDir = fullfile(grn_data_file('publication'), 'group_two_level');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    seedBase = seed_base(arm, phase);
    seeds = seedBase + (1:nReps)';
    fname = sprintf('group_fission_%s_N1000_%dreps_%dgen.mat', arm, nReps, nGens);
    fpath = fullfile(outDir, fname);
    assert(~isfile(fpath), 'run_group_selection:Exists', '%s exists; not overwriting', fpath);

    cfg = publication_config('continuous_output', 'populationSize', 1000, 'numGenerations', nGens);
    cfg.numReplicates = nReps;
    cfg.rngSeedBase = seedBase;
    % Group selection (SI "Extended Materials and Methods", "Group selection"): play and
    % fitness-proportional reproduction within fixed groups of 50, no migration; one fission per
    % generation: a group drawn in proportion to its mean selection weight splits into random
    % halves of 25, each doubled, one refilling it and the other replacing a group drawn
    % uniformly among the others (Traulsen and Nowak 2006, generational form). The random-split
    % control draws the splitting group uniformly. Success: last-generation cooperation >= 0.95.
    switch arm
        case {'fission', 'random_split'}
            cfg.groupSize          = 50;
            cfg.groupPairingMode   = 'within';
            cfg.groupSelectionMode = arm;
        case 'control'
            % no group fields: evolve_population runs the population undivided
        otherwise
            error('run_group_selection:BadArm', 'ARM must be fission | random_split | control');
    end
    cfg.conditionName = sprintf('group_fission_%s_%s_N1000', arm, phase);

    open_pool(nWorkers);
    fprintf('\n%s | %s | %d reps x %d gens | seeds %d..%d\n', arm, phase, nReps, nGens, seeds(1), seeds(end));
    t0 = tic;
    res = evolve_population(cfg);
    elapsed = toc(t0);
    nWorkers = gcp('nocreate').NumWorkers;

    finalCoop = cellfun(@(r) double(r.historyCooperation(end)), res(:));
    success = finalCoop >= 0.95;
    wallPerRep = elapsed * nWorkers / nReps;
    fprintf('%s %s: wall %.1f min | %.1f s per replicate per worker\n', ...
        arm, phase, elapsed / 60, wallPerRep);
    fprintf('%s %s: success %d of %d\n', arm, phase, sum(success), nReps);

    save(fpath, 'res', 'cfg', 'seeds', 'finalCoop', 'success', 'elapsed', ...
         'nWorkers', 'wallPerRep', '-v7.3');
    fprintf('saved %s\n', fpath);
    out = struct('file', fpath, 'nSuccess', sum(success), 'n', nReps, ...
                 'elapsed', elapsed, 'wallPerRep', wallPerRep);
end


function b = seed_base(arm, phase)
% Seed base per arm (columns: group selection, the random-split control, undivided).
    T = struct('main', [770000 780000 790000]);
    assert(isfield(T, phase), 'PHASE must be main');
    col = find(strcmp(arm, {'fission', 'random_split', 'control'}));
    assert(~isempty(col), 'ARM must be fission | random_split | control');
    b = T.(phase)(col);
end


function open_pool(nWorkers)
% A pool of nWorkers processes, one thread each; on Windows at BelowNormal priority, and needs 1.5 GB
% of free memory per worker (both skipped on other systems).
    if ~isempty(gcp('nocreate')), return; end
    maxNumCompThreads(2);
    if ispc
        [~, sys] = memory;
        freeGB = sys.PhysicalMemory.Available / 2^30;
        fprintf('available RAM before the pool: %.2f GB (need %.1f)\n', freeGB, 1.5 * nWorkers);
        assert(freeGB >= 1.5 * nWorkers, 'not enough free RAM for %d workers', nWorkers);
    end
    p = parpool('Processes', nWorkers);
    wait(parfevalOnAll(p, @maxNumCompThreads, 0, 1));
    if ispc
        f = parfevalOnAll(p, @() feature('getpid'), 1);
        pids = [feature('getpid'); fetchOutputs(f)];
        st = system(sprintf(['powershell -NoProfile -Command "Get-Process -Id %s | ForEach-Object ', ...
            '{ $_.PriorityClass = ''BelowNormal'' }"'], strjoin(arrayfun(@num2str, pids', 'uni', 0), ',')));
        fprintf('pool: %d workers, 1 thread each; BelowNormal on PIDs %s (status %d)\n', nWorkers, mat2str(pids'), st);
    else
        fprintf('pool: %d workers, 1 thread each; priority and free memory not changed or checked on this system\n', nWorkers);
    end
end
