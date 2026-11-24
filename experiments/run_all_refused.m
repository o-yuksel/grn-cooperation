function out = run_all_refused(phase, varargin)
% RUN_ALL_REFUSED  How often do built populations lose cooperation when every mutation that ends
% cooperation with the clone is refused? Continues the 100 saturating-contrast populations to
% generation 500 with all such mutations refused (all refused, q = 1), the end point of Fig. 5.
%
% Produces : the arm with all refused in the first 100 runs, for Fig. 5, Table S12, SI S4.6 and main
%            text "How cooperating populations lose cooperation" (read by fig5_refusal_curve and
%            run_partial_refusal)
% Inputs   : publication/saturating_contrast_100reps.mat (its config; run_signaling_architectures);
%            publication/mutation_off/branch.mat, the start populations (run_refusal_start_populations)
% Outputs  : publication/refused_by_class/rruns_<arm>_<from>-<to>.mat (per run the cooperation
%              history and the refusal's log) and split_log_<phase>_<arm>.txt (the printed lines)
% Run      : run_all_refused('timing_run', 'coopend') is a timing run: it runs populations 1-4 and projects the time, then
%              run_all_refused('main', 'coopend') runs 5-100 (about 22 min on 4 workers); the same
%              two calls with 'response'
% Seed     : rngSeedBase = 20261023 + k for population k, Mersenne twister (in evolve_population),
%            as in run_mutation_off; a refused mutation is redrawn from the same stream
%
% The arm 'coopend' refuses, through mutate_with_refusal 'REFUSE_ALL', every mutation that ends
% cooperation with the clone and leaves free those that end only retaliation. The arm 'response'
% refuses only those that end it through the response to the partner (mutate_with_refusal's class
% 1); no paper item reports it, but run_partial_refusal('report') prints its collapses beside the
% partial-refusal arms.

    %% Parameters
    pub = grn_data_file('publication');
    sdir = fullfile(pub, 'refused_by_class');                 % this runner's outputs
    arm = varargin{1};
    if ~any(strcmp(arm, {'coopend', 'response'})), error('run_all_refused:arm', 'arm must be coopend or response'); end
    logFile = fullfile(sdir, ['split_log_' phase '_' arm '.txt']);   % the printed lines, appended
    stamp = char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss');

    %% Inputs
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'config');
    assert(d.config.geneCooperator == 1 && d.config.geneDefector == 2, 'gene order');
    assert(~isfield(d.config, 'geneInputBias') || isempty(d.config.geneInputBias), 'config carries geneInputBias');
    L = load(fullfile(pub, 'mutation_off', 'branch.mat')); br = L.br;

    %% Compute
    if ~exist(sdir, 'dir'), mkdir(sdir); end
    switch phase
        case 'timing_run', [out, rep] = timing_run(d.config, br, sdir, arm);
        case 'main',  [out, rep] = main_runs(d.config, br, sdir, arm, varargin{2:end});
        otherwise, error('run_all_refused:phase', 'unknown phase %s', phase);
    end

    %% Report
    rep = [{sprintf('\n===== RESPONSE SPLIT %s %s %s UTC =====\n', phase, arm, stamp), ...
        sprintf('operator: %s\n', which('mutate_with_refusal'))}, rep];
    fprintf('%s', rep{:});

    %% Save
    fid = fopen(logFile, 'a'); fprintf(fid, '%s', rep{:}); fclose(fid);
end


function cfg = run_cfg(base, b, arm, ngen)
% One continued run of population b.k, as run_mutation_off's run_cfg, with mutate_with_refusal in
% mode 'REFUSE_ALL' (arm 'coopend') or 'REFUSE_RESPONSE' (arm 'response'); populations saved at
% continued generations 1, 41, 151 and the last, which only record.
    cfg = base; cfg.numReplicates = 1; cfg.numGenerations = ngen;
    % the continuation seed of population k, the same in every arm
    cfg.rngSeedBase = 20261023 + b.k;
    cfg.alwaysSaveGenerations = ngen; cfg.historySampleInterval = 100000;
    cfg.initPopulationM = b.M; cfg.initPopulationB = b.B;
    modes = struct('coopend', 'REFUSE_ALL', 'response', 'REFUSE_RESPONSE');
    cfg.mutationOperator = @mutate_with_refusal; cfg.refusalMode = modes.(arm);
    cfg.alwaysSaveGenerations = unique(min([1 41 151 ngen], ngen));
end


function r = r3_one_run(cfg)
% One refusing run: the cooperation history, the share defecting against its own clone at each
% saved generation, and the refusal's log (the same function as in run_mutation_off).
    t = tic; mutate_with_refusal('reset');
    X = evolve_population(cfg); x = X{1}; L = mutate_with_refusal('log');
    r.hc = double(x.historyCooperation(:))'; r.sec = toc(t);
    r.hg = double(x.historyGenerations(:))'; r.defShare = nan(1, numel(r.hg));
    for q = 1:numel(r.hg), r.defShare(q) = def_share(x.historyM(:, :, :, q), x.historyB(:, :, :, q), cfg); end
    r.sampleP = L.sampleP; r.sampleM = L.sampleM; r.sampleS = L.sampleS;
    r.log = rmfield(L, {'sampleP', 'sampleM', 'sampleS'});
end


function sh = def_share(M, B, cfg)
% Share of individuals, among the commonest genotypes covering 95% of the population, whose
% cooperation against their clone is below 0.5.
    N = size(M, 3); [Wu, ~, ic] = unique(reshape(double(M) .* double(B), 16, N)', 'rows'); cnt = accumarray(ic, 1);
    [cs, ord] = sort(cnt, 'descend'); k = find(cumsum(cs) >= 0.95 * N, 1); sh = 0;
    for u = ord(1:k)'
        W = reshape(Wu(u, :), 4, 4);
        if mutate_with_refusal('coop', W, W, cfg) < 0.5, sh = sh + cnt(u); end
    end
    sh = sh / N;
end


function open_pool(nw)
% A pool of nw workers (default 2), one thread each, opened only when none is open.
    if nargin < 1, nw = 2; end
    maxNumCompThreads(1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', nw); end
    wait(parfevalOnAll(pool, @maxNumCompThreads, 0, 1));
end


function t = first_or_inf(v)
% The first index where v is true, Inf if none.
    t = find(v, 1); if isempty(t), t = inf; end
end


function [O, rep] = timing_run(base, br, sdir, arm)
% Populations 1-4 of one arm, and the time it projects for populations 5-100.
    open_pool(4); jobs = arrayfun(@(j) run_cfg(base, br(j), arm, 501 - br(j).t0), 1:4, 'UniformOutput', false);
    R = cell(4, 1); t = tic;
    parfor j = 1:4, R{j} = r3_one_run(jobs{j}); end
    O.sec = toc(t); gens = arrayfun(@(j) 501 - br(j).t0, 1:4);
    lab = [num2cell([br(1:4).k]'), repmat({arm}, 4, 1), num2cell([br(1:4).t0]')];
    save(fullfile(sdir, sprintf('rruns_%s_001-004.mat', arm)), 'R', 'lab');
    rest = sum(501 - [br(5:end).t0]); perGen = max(cellfun(@(r) r.sec, R)) / max(gens);
    rep = {sprintf('[TIMING RUN %s] reps 1-4 (%s generations) on 4 workers: wall %.1f min; per run %s s; finals %s; refused first draws %s\n', ...
        upper(arm), mat2str(gens), O.sec / 60, mat2str(round(cellfun(@(r) r.sec, R))'), mat2str(cellfun(@(r) r.hc(end), R)', 3), ...
        mat2str(cellfun(@(r) r.log.refFirst / r.log.mut, R)', 3))};
    rep{end + 1} = sprintf('[TIMING RUN %s] ~%.3f s per generation per worker; reps 5-100 (%d generations) on 4 workers: ~%.0f min\n', ...
        upper(arm), perGen, rest, perGen * rest / 4 / 60);
end


function [O, rep] = main_runs(base, br, sdir, arm, chunk)
% Populations 5-100 of one arm in chunks (default 16); a chunk already saved is skipped.
    if nargin < 5, chunk = 16; end
    open_pool(4); tAll = tic; rep = {};
    for c0 = 5:chunk:numel(br)
        idx = c0:min(c0 + chunk - 1, numel(br));
        fn = fullfile(sdir, sprintf('rruns_%s_%03d-%03d.mat', arm, br(idx(1)).k, br(idx(end)).k));
        if exist(fn, 'file'), rep{end + 1} = sprintf('skip %s (exists)\n', fn); continue; end %#ok<AGROW>
        jobs = arrayfun(@(j) run_cfg(base, br(j), arm, 501 - br(j).t0), idx, 'UniformOutput', false);
        lab = [num2cell([br(idx).k]'), repmat({arm}, numel(idx), 1), num2cell([br(idx).t0]')];
        R = cell(numel(jobs), 1); t = tic;
        parfor j = 1:numel(jobs), R{j} = r3_one_run(jobs{j}); end
        for j = 1:numel(jobs)
            h = R{j}.hc; cg = first_or_inf(h < 0.5);
            rep{end + 1} = sprintf('rep %3d %s t0 %3d gens %3d final %.3f collapse at %s refused %.3f (%.0f s)\n', lab{j, 1}, ...
                arm, lab{j, 3}, numel(h), h(end), num2str(cg), R{j}.log.refFirst / R{j}.log.mut, R{j}.sec); %#ok<AGROW>
        end
        save(fn, 'R', 'lab');
        rep{end + 1} = sprintf('chunk %s: %.1f min (total %.1f min)\n', fn, toc(t) / 60, toc(tAll) / 60); %#ok<AGROW>
    end
    O = [];
end
