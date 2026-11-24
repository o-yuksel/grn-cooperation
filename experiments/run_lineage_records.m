function out = run_lineage_records(mode)
% RUN_LINEAGE_RECORDS  Reruns replicates 1 to 20 of the four signaling architectures and records
% which distinct wiring each individual carries at each generation.
%
% Produces : Fig. 2 and Table S4 (the wiring pointer idx and parent record par that
%            selection_inherited, selection_around_collapse and fig2_draw read)
% Inputs   : mode      'full' (replicates 1 to 20) or 'sample' (replicate 1 of each
%                      architecture; the default); both on a pool of 2 workers, opened if none is
%                      open
%            publication/<architecture>_100reps.mat   (results, config; run_signaling_architectures)
% Outputs  : publication/selection_over_time/lineage_records_<arch>_r<NN>.mat   (an existing file is skipped);
%            one printed line per run
% Run      : run_lineage_records('full')   (long, not timed here; 'sample' about 3 min)
% Seed     : replicate r reruns under its stored run's seed, rng(12345 + r, 'twister')
%
% A rerun that differs from the stored run is saved without wiring data (check_gate). cacheW
% lists each distinct wiring W = M .* B once; idx(i, g) points individual i at generation g
% into it.

    if nargin < 1, mode = 'sample'; end
    cfg = struct( ...
        'pub', grn_data_file('publication'), ...
        'archs', {{'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'}}, ...
        'nGen', 500, 'seedBase', 12345, 'nWorkers', 2);
    cfg.matDir = fullfile(cfg.pub, 'selection_over_time');
    switch mode
        case 'sample',  out = run_set(cfg, 1);
        case 'full',    out = run_set(cfg, 1:20);
        otherwise, error('unknown mode %s', mode);
    end
end


% =========================================================================
function out = run_set(cfg, reps)
    jobs = build_jobs(cfg, reps);
    fprintf('%d runs to do (%d already on disk skipped)\n', numel(jobs), ...
        numel(cfg.archs) * numel(reps) - numel(jobs));
    t0 = tic; res = cell(numel(jobs), 1);
    if ~isempty(jobs)
        if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
        parfor j = 1:numel(jobs), res{j} = run_job(jobs{j}); end
    end
    for j = 1:numel(res)
        if isempty(res{j}), continue; end
        t = res{j}.t;
        fprintf('%-20s r%02d  run %5.1f s  unique %4.1f s  total %6.1f s | wirings %5d | gate %s\n', ...
            res{j}.arch, res{j}.r, t.run, t.unique, t.total, res{j}.nCache, res{j}.gate.text);
    end
    fprintf('wall time %.1f min\n', toc(t0) / 60);
    out = res;
end


function jobs = build_jobs(cfg, reps)
% One job per stored replicate without a file yet: its stored config, rerun for 500 generations
% with every generation recorded, and the stored history the identity check compares with.
    jobs = {};
    for a = 1:numel(cfg.archs)
        todo = reps(arrayfun(@(r) ~exist(out_file(cfg, a, r), 'file'), reps));
        if isempty(todo), continue; end
        d = load(fullfile(cfg.pub, [cfg.archs{a} '_100reps.mat']), 'results', 'config');
        for r = todo
            s = d.results{r};
            keep = double(s.historyGenerations(:))' <= cfg.nGen;
            st = struct('hc', s.historyCooperation(1:cfg.nGen), ...
                'par', s.historyParents(:, 1:cfg.nGen), ...
                'gens', double(s.historyGenerations(keep)), ...
                'M', s.historyM(:,:,:,keep), 'B', s.historyB(:,:,:,keep));
            c = d.config;
            c.numGenerations = cfg.nGen; c.numReplicates = 1;
            c.historySampleInterval = 1;
            c.rngSeedBase = cfg.seedBase + r - 1;      % repIdx = 1 -> rng(12345 + r)
            jobs{end+1} = struct('arch', cfg.archs{a}, 'r', r, 'config', c, 'stored', st, ...
                'file', out_file(cfg, a, r)); %#ok<AGROW>
        end
    end
end


function f = out_file(cfg, a, r)
    f = fullfile(cfg.matDir, sprintf('lineage_records_%s_r%02d.mat', cfg.archs{a}, r));
end


function g = check_gate(x, s, nG)
% The identity check: cooperation history, parent record and every stored snapshot of M and B
% equal the stored run's.
    g.hc = isequal(x.historyCooperation(1:nG), s.hc);
    g.par = isequal(x.historyParents(:, 1:nG), s.par);
    [tf, loc] = ismember(s.gens, double(x.historyGenerations(:))');
    g.snap = all(tf) && isequal(x.historyM(:,:,:,loc), s.M) && isequal(x.historyB(:,:,:,loc), s.B);
    g.pass = g.hc && g.par && g.snap;
    dv = find(x.historyCooperation(1:nG) ~= s.hc, 1);
    g.firstDiv = min([dv(:); inf]);
    g.text = sprintf('cooperation %s, parents %s, %d snapshots %s -> %s', pf(g.hc), pf(g.par), ...
        numel(s.gens), pf(g.snap), ternary(g.pass, 'IDENTICAL', 'DIFFERS'));
end


% =========================================================================
function res = run_job(job)
    maxNumCompThreads(1);
    tAll = tic;
    t = tic;
    X = evolve_population(job.config); x = X{1};
    tm.run = toc(t);
    G = job.config.numGenerations;
    gate = check_gate(x, job.stored, G);
    hc = double(x.historyCooperation(:));
    par = x.historyParents(:, 1:G);
    res = struct('arch', job.arch, 'r', job.r, 'gate', gate, 'nCache', 0);
    if ~gate.pass
        tm.unique = 0; tm.total = toc(tAll); res.t = tm;
        save_if_new(job.file, struct('res', res, 'hc', hc, 'hcStored', double(job.stored.hc(:))));
        return;
    end

    % Distinct wirings of the whole run, cached once; idx(i, g) points into the cache.
    t = tic;
    N = size(x.historyM, 3);
    Wall = zeros(N * G, 16);
    for g = 1:G
        Wall((g - 1) * N + (1:N), :) = reshape(double(x.historyM(:,:,:,g)) .* double(x.historyB(:,:,:,g)), 16, N)';
    end
    [cacheW, ~, ic] = unique(Wall, 'rows', 'stable');
    idx = uint32(reshape(ic, N, G));
    clear Wall
    nC = size(cacheW, 1);
    tm.unique = toc(t);

    tm.total = toc(tAll);
    res.nCache = nC; res.t = tm;
    S = struct('res', res, 'hc', hc, 'hcStored', double(job.stored.hc(:)), 'par', par, 'idx', idx, ...
        'cacheW', cacheW);
    save_if_new(job.file, S);
end


function save_if_new(f, S)
    if exist(f, 'file')
        warning('run_lineage_records:exists', '%s appeared during the run; not overwritten', f);
        return;
    end
    save(f, '-struct', 'S', '-v7.3');
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
