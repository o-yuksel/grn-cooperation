function out = run_refusal_start_populations(phase, varargin)
% RUN_REFUSAL_START_POPULATIONS  Rebuilds each saturating-contrast population at the generation its
% cooperation first exceeds 0.90, where the continued runs with mutation on, off or refused begin.
%
% Produces : the start populations of the continued runs behind Fig. 5 and Table S12 (SI S4.6), and
%            of Table S11 (read by run_mutation_off, run_all_refused, run_partial_refusal and
%            tableS11_mutant_encounters; fig5_refusal_curve draws the continued runs)
% Inputs   : publication/saturating_contrast_100reps.mat   (written by run_signaling_architectures)
% Outputs  : publication/refusal_start_populations/branch.mat   (br: per population k the start generation t0,
%              its connections M and weights B at t0, and same, true when the rerun's cooperation
%              history to t0 equals the stored one) and its byte-identical copy
%              publication/mutation_off/branch.mat, the file its readers load
%            publication/refusal_start_populations/main_g2_reps<a>-<b>.mat   ('main': per run the cooperation
%              history hc, the final value, the first generation below 0.5 and final < 0.95)
%            publication/refusal_start_populations/log_<phase>.txt   (the lines this run prints)
% Run      : run_refusal_start_populations('branch', 100)   (all 100 populations; 2 workers);
%            then ('main', 2, 1:10) and ('main', 2, 11:100)   (read by run_mutation_off('report'))
% Seed     : population k reruns as one replicate with rngSeedBase 12345 + k - 1, so under its
%            stored seed rng(12345 + k, 'twister') (set in evolve_population); 'main' continues
%            population k with rngSeedBase 20261023 + k in every arm, run_mutation_off's seed
%
% All 100 populations exceed 0.90; the continued runs ask how populations that built cooperation
% lose it. The rerun stops at t0, so the start is the stored run's own population; same records the
% check. 'main' continues each start population for 150 generations in four arms (mutation on,
% mutation off, and two arms with a gain g on the partner's signal that no paper item reports);
% run_mutation_off('report') checks its own first 150 continued generations against the first two.

    %% Parameters
    pub = grn_data_file('publication');
    odir = fullfile(pub, 'refusal_start_populations');                       % this runner's outputs
    logFile = fullfile(odir, ['log_' phase '.txt']);   % the printed lines, appended
    stamp = char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss');

    %% Inputs
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'results', 'config');

    %% Compute
    if ~exist(odir, 'dir'), mkdir(odir); end
    switch phase
        case 'branch', [out, rep] = branch(d, odir, pub, varargin{:});
        case 'main',   [out, rep] = main_runs(d, odir, varargin{:});
        otherwise, error('run_refusal_start_populations:phase', 'phase must be branch or main');
    end

    %% Report
    rep = [{sprintf('\n===== REFUSAL START POPULATIONS %s %s UTC =====\n', phase, stamp)}, rep];
    fprintf('%s', rep{:});

    %% Save
    fid = fopen(logFile, 'a'); fprintf(fid, '%s', rep{:}); fclose(fid);
end


function [br, rep] = branch(d, odir, pub, nRep)
% Reruns each replicate to the generation t0 it first exceeds 0.90, keeps its population at t0, and
% records whether the cooperation history 1 to t0 equals the stored one. Saves branch.mat and its
% copy in mutation_off/.
    jobs = {}; ks = [];
    for k = 1:nRep
        h = double(d.results{k}.historyCooperation(:))'; t0 = find(h > 0.90, 1);
        if isempty(t0), continue; end
        cfg = d.config; cfg.numReplicates = 1; cfg.rngSeedBase = 12345 + k - 1;
        cfg.numGenerations = t0; cfg.alwaysSaveGenerations = t0; cfg.historySampleInterval = 100000;
        jobs{end + 1} = cfg; ks(end + 1) = k; %#ok<AGROW>
    end
    open_pool(); R = cell(numel(jobs), 1);
    parfor j = 1:numel(jobs)
        t = tic; X = evolve_population(jobs{j}); x = X{1};
        R{j} = struct('M', x.historyM(:, :, :, end), 'B', x.historyB(:, :, :, end), ...
            'hc', double(x.historyCooperation(:))', 'sec', toc(t));
    end
    br = struct('k', {}, 't0', {}, 'M', {}, 'B', {}, 'same', {}); rep = {};
    for j = 1:numel(jobs)
        k = ks(j); t0 = jobs{j}.numGenerations; hs = double(d.results{k}.historyCooperation(1:t0))';
        same = isequal(R{j}.hc, hs);
        br(j) = struct('k', k, 't0', t0, 'M', R{j}.M, 'B', R{j}.B, 'same', same);
        rep{end + 1} = sprintf('rep %2d: branch generation %3d, coop history 1:%d identical to stored %d, coop at branch %.3f, %.0f s (%.2f s/gen)\n', ...
            k, t0, t0, same, R{j}.hc(end), R{j}.sec, R{j}.sec / t0); %#ok<AGROW>
    end
    save(fullfile(odir, 'branch.mat'), 'br');
    cdir = fullfile(pub, 'mutation_off'); if ~exist(cdir, 'dir'), mkdir(cdir); end
    copyfile(fullfile(odir, 'branch.mat'), fullfile(cdir, 'branch.mat'));
    rep{end + 1} = sprintf('All identical: %d\n', all([br.same]));
end


% Arms of 'main': 'control' (mutation on); 'G' (mutation rate 0); 'A' (the partner's signal into the sensors
% scaled by g); 'B' (every weight B and K_M divided by g instead).
function cfg = arm_cfg(base, arm, g)
    cfg = base;
    switch arm
        case 'control'
        case 'G', cfg.mutationRate = 0;
        case 'A', cfg.sensorIntegrationRate = g;
        case 'B', cfg.michaelisConstants = ones(cfg.numGenes, 1) / g;
        otherwise, error('arm');
    end
end


function s = arm_scale(arm, g)
% Factor applied to every stored weight B (and partner wiring) in this arm.
    s = 1; if strcmp(arm, 'B'), s = 1 / g; end
end


function cfg = run_cfg(base, br1, arm, g, ngen)
    cfg = arm_cfg(base, arm, g); cfg.numReplicates = 1; cfg.numGenerations = ngen;
    cfg.rngSeedBase = 20261023 + br1.k;   % same stream for every arm of this population
    cfg.alwaysSaveGenerations = ngen; cfg.historySampleInterval = 100000;
    cfg.initPopulationM = br1.M; cfg.initPopulationB = single(double(br1.B) * arm_scale(arm, g));
end


function r = one_run(cfg)
    t = tic; X = evolve_population(cfg); x = X{1};
    r.hc = double(x.historyCooperation(:))'; r.final = r.hc(end);
    r.collapse = first_or_inf(r.hc < 0.5); r.fail = r.final < 0.95;
    r.M = x.finalNetworkM; r.B = x.finalNetworkB; r.sec = toc(t);
end


function [O, rep] = main_runs(d, odir, g, reps, ngen)
% Start populations reps, ngen generations (default 150) in the four arms; saves
% main_g<g>_reps<a>-<b>.mat.
    if nargin < 5, ngen = 150; end
    L = load(fullfile(odir, 'branch.mat')); br = L.br(reps);
    arms = {'control', 'G', 'A', 'B'}; jobs = {}; lab = {};
    for j = 1:numel(br)
        for a = 1:4, jobs{end + 1} = run_cfg(d.config, br(j), arms{a}, g, ngen); lab(end + 1, :) = {br(j).k, arms{a}}; end %#ok<AGROW>
    end
    open_pool(); R = cell(numel(jobs), 1); t = tic;
    parfor j = 1:numel(jobs), R{j} = one_run(jobs{j}); end
    rep = {};
    for j = 1:numel(jobs)
        rep{end + 1} = sprintf('rep %2d %-7s final %.3f collapse gen %s fail %d (%.0f s)\n', lab{j, 1}, lab{j, 2}, R{j}.final, ...
            num2str(R{j}.collapse), R{j}.fail, R{j}.sec); %#ok<AGROW>
    end
    rep{end + 1} = sprintf('Wall time %.1f min\n', toc(t) / 60);
    fn = fullfile(odir, sprintf('main_g%g_reps%d-%d.mat', g, reps(1), reps(end)));
    for j = 1:numel(R), R{j} = rmfield(R{j}, {'M', 'B'}); end   % keep the .mat small
    save(fn, 'R', 'lab', 'g', 'ngen'); O = R;
end


function t = first_or_inf(v)
    t = find(v, 1); if isempty(t), t = inf; end
end


function open_pool()
    maxNumCompThreads(1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', 2); end
    wait(parfevalOnAll(pool, @maxNumCompThreads, 0, 1));
end
