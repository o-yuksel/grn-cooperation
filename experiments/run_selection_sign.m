function out = run_selection_sign(mode)
% RUN_SELECTION_SIGN  Does selection favor cooperation or defection, generation by generation?
% Reruns replicates 1 to 20 of the four signaling architectures and stores the realized selection
% on cooperation, S (defined at s_of), at every generation.
%
% Produces : Fig. 2 and Table S4 (the stored S that run_selection_z, selection_inherited,
%            selection_around_collapse and fig2_draw read)
% Inputs   : mode      'full' (replicates 1 to 20) or 'sample' (replicate 1 of each
%                      architecture; the default); both on a pool of 2 workers, opened if none is
%                      open; every call first runs the unit test of S on hand-built examples
%            publication/<architecture>_100reps.mat   (results, config; run_signaling_architectures)
% Outputs  : publication/selection_over_time/selection_sign_<arch>_r<NN>.mat   (an existing file is
%            skipped)
% Run      : run_selection_sign('full')   (long, not timed here; 'sample' about 2 min)
% Seed     : replicate r reruns under its stored run's seed, rng(12345 + r, 'twister')
%
% S is stored only if the rerun's cooperation history and parent record equal the stored ones.

    if nargin < 1, mode = 'sample'; end
    cfg = struct( ...
        'pub', grn_data_file('publication'), ...
        'archs', {{'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'}}, ...
        'nGen', 500, 'seedBase', 12345, 'nWorkers', 2);
    cfg.matDir = fullfile(cfg.pub, 'selection_over_time');
    unit_test();
    switch mode
        case 'sample',  out = run_set(cfg, 1);
        case 'full',    out = run_set(cfg, 1:20);
        otherwise, error('unknown mode %s', mode);
    end
end


% =========================================================================
function [S, n] = s_of(par, z)
% Selection term Cov(n, z), population covariance (1/N). par: parent index of each slot of the
% next generation (one column of historyParents); z: cooperation of the N individuals of this
% generation. At generation g, z_i is i's cooperation in play (snapshotCoopVectors) and n_i the
% number of times i appears in historyParents(:, g), the parents of generation g + 1. S includes
% partner luck, which inherited selection S_G (selection_inherited) removes.
    N = numel(z);
    n = accumarray(double(par(:)), 1, [N 1]);
    z = double(z(:));
    S = mean((n - mean(n)) .* (z - mean(z)));
end


function unit_test()
% Hand-built examples, each checked against a hand value. N = 4.
    [S, n] = s_of([1 1 1 4], [0 1 1 0.5]);          % n = [3 0 0 1]; hand value -0.5
    assert(isequal(n, [3; 0; 0; 1]) && abs(S - (-0.5)) < 1e-15, 'unit test 1 failed');
    [S, n] = s_of([2 3 3 4], [0 1 0.5 0.5]);        % n = [0 1 2 1]; dev n [-1 0 1 0], dev z [-.5 .5 0 0]: 0.5/4
    assert(isequal(n, [0; 1; 2; 1]) && abs(S - 0.125) < 1e-15, 'unit test 2 failed');
    [S, n] = s_of([1 2 3 4], [0 1 0.5 0.5]);        % everyone leaves one offspring: S = 0
    assert(sum(n) == 4 && S == 0, 'unit test 3 failed');
    [S, ~] = s_of([1 2 3 4], [0.3 0.3 0.3 0.3]);    % zero variance in z: S = 0
    assert(S == 0, 'unit test 4 failed');
    fprintf('unit test of S: PASS (hand values -0.5, 0.125, 0, 0)\n');
end


% =========================================================================
function out = run_set(cfg, reps)
    jobs = build_jobs(cfg, reps);
    fprintf('%d runs to do (%d already on disk skipped)\n', numel(jobs), numel(cfg.archs) * numel(reps) - numel(jobs));
    t0 = tic; res = cell(numel(jobs), 1);
    if ~isempty(jobs)
        if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
        parfor j = 1:numel(jobs), res{j} = run_job(jobs{j}); end
    end
    for j = 1:numel(res)
        if isempty(res{j}), continue; end
        r = res{j};
        fprintf('%-20s r%02d  %5.1f s  gate %s | sum(n)==N at %d/500 | NaN z %d | S(1) %+.5f S(7) %+.5f\n', ...
            r.arch, r.r, r.sec, r.gate.text, r.nSumOK, r.nNaN, r.S1, r.S7);
    end
    fprintf('wall time %.1f min\n', toc(t0) / 60);
    out = res;
end


function jobs = build_jobs(cfg, reps)
    jobs = {};
    for a = 1:numel(cfg.archs)
        todo = reps(arrayfun(@(r) ~exist(out_file(cfg, a, r), 'file'), reps));
        if isempty(todo), continue; end
        d = load(fullfile(cfg.pub, [cfg.archs{a} '_100reps.mat']), 'results', 'config');
        for r = todo
            s = d.results{r};
            st = struct('hc', s.historyCooperation(1:cfg.nGen), 'par', s.historyParents(:, 1:cfg.nGen));
            c = d.config;
            c.numGenerations = cfg.nGen; c.numReplicates = 1;
            c.historySampleInterval = 1;
            c.snapshotCoopGenerations = 1:cfg.nGen;
            c.rngSeedBase = cfg.seedBase + r - 1;      % repIdx = 1 -> rng(12345 + r)
            jobs{end+1} = struct('arch', cfg.archs{a}, 'r', r, 'config', c, 'stored', st, ...
                'file', out_file(cfg, a, r)); %#ok<AGROW>
        end
    end
end


function f = out_file(cfg, a, r)
    f = fullfile(cfg.matDir, sprintf('selection_sign_%s_r%02d.mat', cfg.archs{a}, r));
end


function res = run_job(job)
    maxNumCompThreads(1);
    t = tic;
    X = evolve_population(job.config); x = X{1};
    G = job.config.numGenerations;
    gate.hc = isequal(x.historyCooperation(1:G), job.stored.hc);
    gate.par = isequal(x.historyParents(:, 1:G), job.stored.par);
    gate.pass = gate.hc && gate.par;
    gate.text = sprintf('cooperation %s, parents %s -> %s', pf(gate.hc), pf(gate.par), ...
        ternary(gate.pass, 'IDENTICAL', 'DIFFERS'));
    hc = double(x.historyCooperation(1:G)); hc = hc(:);
    res = struct('arch', job.arch, 'r', job.r, 'gate', gate, 'nSumOK', 0, 'nNaN', 0, 'S1', NaN, 'S7', NaN);
    P = struct('res', res, 'hc', hc, 'hcStored', double(job.stored.hc(:)));
    if gate.pass
        Z = x.snapshotCoopVectors; N = size(x.historyParents, 1);
        if ~isequal(size(Z), [N G]) || ~isequal(double(x.snapshotCoopGenerations(:)'), 1:G)
            error('snapshotCoopVectors is %s, expected %d x %d with generations 1:%d', mat2str(size(Z)), N, G, G);
        end
        S = zeros(G, 1); meanZ = zeros(G, 1); zeroVar = false(G, 1); nSumOK = 0;
        for g = 1:G
            [S(g), n] = s_of(x.historyParents(:, g), Z(:, g));
            meanZ(g) = mean(Z(:, g));
            zeroVar(g) = all(Z(:, g) == Z(1, g));
            nSumOK = nSumOK + (sum(n) == N);
        end
        res.nSumOK = nSumOK; res.nNaN = sum(isnan(Z(:))); res.S1 = S(1); res.S7 = S(7);
        P.S = S; P.meanZ = meanZ; P.zeroVar = zeroVar; P.N = N;
    end
    res.sec = toc(t);
    P.res = res;
    if exist(job.file, 'file')
        warning('run_selection_sign:exists', '%s appeared during the run; not overwritten', job.file);
    else
        save(job.file, '-struct', 'P', '-v7.3');
    end
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
