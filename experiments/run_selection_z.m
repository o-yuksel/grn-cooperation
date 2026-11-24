function out = run_selection_z(mode)
% RUN_SELECTION_Z  Reruns replicates 1 to 20 of the four signaling architectures and stores each
% individual's cooperation at every generation, with the parent record.
%
% Produces : Fig. 2 and Table S4 (the per-individual cooperation Z that selection_inherited,
%            selection_around_collapse and fig2_draw read)
% Inputs   : mode   'full' (replicates 1 to 20) or 'sample' (replicate 1 of each architecture;
%                   the default); both on a pool of 2 workers, opened if none is open
%            publication/<architecture>_100reps.mat   (results, config; run_signaling_architectures)
%            publication/selection_over_time/selection_sign_<arch>_r<NN>.mat   (S; written by run_selection_sign)
% Outputs  : publication/selection_over_time/selection_z_<arch>_r<NN>.mat   (Z, N-by-500, column g each
%            individual's cooperation in play at generation g; par = historyParents; `gate`;
%            maxDiffS); an existing file is skipped
% Run      : run_selection_z('full') after run_selection_sign('full')   (long, not timed here;
%            'sample' about 2 min)
% Seed     : replicate r reruns under its stored run's seed, rng(12345 + r, 'twister')
%
% A file is written only if the rerun's cooperation history and parent record equal the stored
% ones (else the run prints gate FAIL); if S recomputed from Z differs from the stored S by more
% than 1e-12 (maxDiffS), the run stops before saving.
    if nargin < 1, mode = 'sample'; end
    cfg = struct( ...
        'pub', grn_data_file('publication'), ...
        'archs', {{'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'}}, ...
        'nGen', 500, 'seedBase', 12345, 'nWorkers', 2);
    cfg.matDir = fullfile(cfg.pub, 'selection_over_time');
    if strcmp(mode, 'full'), reps = 1:20; else, reps = 1; end
    jobs = build_jobs(cfg, reps);
    fprintf('%d runs to do\n', numel(jobs));
    t0 = tic; res = cell(numel(jobs), 1);
    if ~isempty(jobs)
        if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
        parfor j = 1:numel(jobs), res{j} = run_job(jobs{j}); end
    end
    for j = 1:numel(res)
        r = res{j};
        fprintf('%-20s r%02d  %5.1f s  gate %s | max|S - S stored| = %.3g (%s)\n', r.arch, r.r, r.sec, ...
            ternary(r.pass, 'PASS', 'FAIL'), r.maxDiffS, ternary(r.maxDiffS <= 1e-12, 'PASS', 'FAIL'));
    end
    fprintf('wall time %.1f min\n', toc(t0) / 60);
    out = res;
end

function jobs = build_jobs(cfg, reps)
    jobs = {};
    for a = 1:numel(cfg.archs)
        todo = reps(arrayfun(@(r) ~exist(zfile(cfg, a, r), 'file'), reps));
        if isempty(todo), continue; end
        d = load(fullfile(cfg.pub, [cfg.archs{a} '_100reps.mat']), 'results', 'config');
        for r = todo
            s = d.results{r};
            st = struct('hc', s.historyCooperation(1:cfg.nGen), 'par', s.historyParents(:, 1:cfg.nGen));
            c = d.config;
            c.numGenerations = cfg.nGen; c.numReplicates = 1; c.historySampleInterval = 1;
            c.snapshotCoopGenerations = 1:cfg.nGen;
            c.rngSeedBase = cfg.seedBase + r - 1;
            jobs{end+1} = struct('arch', cfg.archs{a}, 'r', r, 'config', c, 'stored', st, ...
                'file', zfile(cfg, a, r), 'sfile', fullfile(cfg.matDir, sprintf('selection_sign_%s_r%02d.mat', cfg.archs{a}, r))); %#ok<AGROW>
        end
    end
end

function f = zfile(cfg, a, r)
    f = fullfile(cfg.matDir, sprintf('selection_z_%s_r%02d.mat', cfg.archs{a}, r));
end

function res = run_job(job)
    maxNumCompThreads(1);
    t = tic; G = job.config.numGenerations;
    X = evolve_population(job.config); x = X{1};
    gate.hc = isequal(x.historyCooperation(1:G), job.stored.hc);
    gate.par = isequal(x.historyParents(:, 1:G), job.stored.par);
    pass = gate.hc && gate.par;
    res = struct('arch', job.arch, 'r', job.r, 'pass', pass, 'maxDiffS', Inf, 'sec', 0);
    if pass
        Z = double(x.snapshotCoopVectors); par = x.historyParents(:, 1:G); N = size(Z, 1);
        if ~isequal(size(Z), [N G]), error('Z size %s', mat2str(size(Z))); end
        S = zeros(G, 1);
        for g = 1:G
            n = accumarray(double(par(:, g)), 1, [N 1]);
            S(g) = mean((n - mean(n)) .* (Z(:, g) - mean(Z(:, g))));
        end
        old = load(job.sfile, 'S');
        res.maxDiffS = max(abs(S - old.S));
        if ~(res.maxDiffS <= 1e-12)
            error('run_selection_z:gate', '%s r%02d: S from Z differs from the stored S by %.3g (> 1e-12); not saved', ...
                job.arch, job.r, res.maxDiffS);
        end
        P = struct('Z', Z, 'par', par, 'gate', gate, 'maxDiffS', res.maxDiffS);
        if ~exist(job.file, 'file'), save(job.file, '-struct', 'P', '-v7.3'); end
    end
    res.sec = toc(t);
end
