function out = run_mutation_off(phase, varargin)
% RUN_MUTATION_OFF  Do populations that built cooperation lose it less often when mutation stops?
% Continues the 100 saturating-contrast populations to generation 500 with mutation on and with
% mutation off, and counts the collapses of the cooperation they built.
%
% Produces : the arms mutation on and mutation off of the first 100 runs, for Fig. 5, Table S12,
%            SI S4.6, Abstract, Significance, main text "How cooperating populations lose
%            cooperation" and Methods "Failure route" (fig5_refusal_curve reads them)
% Inputs   : publication/saturating_contrast_100reps.mat (run_signaling_architectures);
%            publication/mutation_off/branch.mat, the start populations (a copy of refusal_start_populations/
%              branch.mat, run_refusal_start_populations('branch', 100)); refusal_start_populations/main_g2_reps*.mat
%              ('report' only; run_refusal_start_populations('main', 2, ...))
% Outputs  : publication/mutation_off/runs_<from>-<to>.mat (each run's cooperation history hc) and
%              report.mat (per population cg, the first continued generation below 0.5, Inf if none,
%              and fin, the final cooperation); publication/refuse_less_conditional/rruns_*.mat (third
%              arm) and report.mat (also T.control and T.off, the file fig5_refusal_curve reads);
%              each phase's printed lines in mutation_off_log_<phase>.txt or refuse_less_conditional_log_<phase>.txt
% Run      : run_mutation_off('main', 10) (both arms in chunks of 10), then ('report')   (about 50 min
%              on 2 workers); then the third arm: 'timing_run' (populations 1-4),
%              ('refuse_less_conditional', 16) (the rest), 'report3'   (about 40 min on 4 workers)
% Seed     : rngSeedBase = 20261023 + k for population k, Mersenne twister (in evolve_population),
%            the same in every arm
%
% Each run starts at the first generation t0 its cooperation exceeds 0.90 (all do) and runs 501 - t0
% generations, to generation 500 of its original run. Mutation off sets the mutation rate to zero
% and changes nothing else. A collapse is the first drop below 0.5. The third arm keeps mutation on
% and refuses every mutation that lowers the network's gap (its conditionality; a lower gap is a less
% conditional network); no paper item reports it, but 'report3'
% writes refuse_less_conditional/report.mat, which carries the counts of mutation on and off that
% fig5_refusal_curve draws.

    %% Parameters
    pub = grn_data_file('publication');
    odir = fullfile(pub, 'mutation_off');                      % arms mutation on and off
    rdir = fullfile(pub, 'refuse_less_conditional');                % third arm
    third = any(strcmp(phase, {'timing_run', 'refuse_less_conditional', 'report3'}));
    if third, logFile = fullfile(rdir, ['refuse_less_conditional_log_' phase '.txt']);   % the printed lines, appended
    else, logFile = fullfile(odir, ['mutation_off_log_' phase '.txt']); end
    stamp = char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss');

    %% Inputs
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'results', 'config');
    L = load(fullfile(odir, 'branch.mat')); br = L.br;

    %% Compute
    if third && ~exist(rdir, 'dir'), mkdir(rdir); end
    switch phase
        case 'main',    [out, rep] = main_runs(d.config, br, odir, varargin{:});
        case 'report',  [out, rep] = report(d, br, odir, pub);
        case 'timing_run', [out, rep] = r3_timing_run(d.config, br, rdir);
        case 'refuse_less_conditional', [out, rep] = r3_main(d.config, br, rdir, varargin{:});
        case 'report3', [out, rep] = r3_report(d, br, odir, rdir);
        otherwise, error('run_mutation_off:phase', 'phase must be main, report, timing_run, refuse_less_conditional or report3');
    end

    %% Report
    rep = [{sprintf('\n===== MUTATION OFF %s %s UTC =====\n', phase, stamp)}, rep];
    fprintf('%s', rep{:});

    %% Save
    fid = fopen(logFile, 'a'); fprintf(fid, '%s', rep{:}); fclose(fid);
end


function cfg = run_cfg(base, b, arm, ngen)
% One continued run of population b.k for arm 'control' (mutation on), 'off' or 'refuse_less_conditional'.
    cfg = base; cfg.numReplicates = 1; cfg.numGenerations = ngen;
    if strcmp(arm, 'off'), cfg.mutationRate = 0; end
    % the continuation seed of population k, the same in every arm
    cfg.rngSeedBase = 20261023 + b.k;
    cfg.alwaysSaveGenerations = ngen; cfg.historySampleInterval = 100000;
    cfg.initPopulationM = b.M; cfg.initPopulationB = b.B;
    % third arm: the refusing mutation step, and populations saved at continued generations 1, 41, 151
    % and the last, which only record
    if strcmp(arm, 'refuse_less_conditional')
        cfg.mutationOperator = @mutate_with_refusal; cfg.refusalMode = 'REFUSE_LESS_CONDITIONAL';
        cfg.alwaysSaveGenerations = unique(min([1 41 151 ngen], ngen));
    end
end


function r = one_run(cfg)
    t = tic; X = evolve_population(cfg); x = X{1};
    r.hc = double(x.historyCooperation(:))'; r.sec = toc(t);
end


function [O, rep] = main_runs(base, br, odir, chunk)
% Chunks of `chunk` populations (both arms); a chunk already saved is skipped.
    arms = {'control', 'off'}; open_pool(); tAll = tic; rep = {};
    for c0 = 1:chunk:numel(br)
        idx = c0:min(c0 + chunk - 1, numel(br));
        fn = fullfile(odir, sprintf('runs_%03d-%03d.mat', br(idx(1)).k, br(idx(end)).k));
        if exist(fn, 'file'), rep{end + 1} = sprintf('skip %s (exists)\n', fn); continue; end %#ok<AGROW>
        jobs = {}; lab = {};
        for j = idx
            for a = 1:2
                jobs{end + 1} = run_cfg(base, br(j), arms{a}, 501 - br(j).t0); %#ok<AGROW>
                lab(end + 1, :) = {br(j).k, arms{a}, br(j).t0}; %#ok<AGROW>
            end
        end
        R = cell(numel(jobs), 1); t = tic;
        parfor j = 1:numel(jobs), R{j} = one_run(jobs{j}); end
        for j = 1:numel(jobs)
            h = R{j}.hc; cg = find(h < 0.5, 1); if isempty(cg), cg = inf; end
            rep{end + 1} = sprintf('rep %3d %-7s t0 %3d gens %3d final %.3f collapse at %s (%.0f s)\n', lab{j, 1}, lab{j, 2}, ...
                lab{j, 3}, numel(h), h(end), num2str(cg), R{j}.sec); %#ok<AGROW>
        end
        save(fn, 'R', 'lab');
        rep{end + 1} = sprintf('chunk %s: %.1f min (total %.1f min)\n', fn, toc(t) / 60, toc(tAll) / 60); %#ok<AGROW>
    end
    O = [];
end


function [O, rep] = report(d, br, odir, pub)
% Counts per arm with Clopper-Pearson 95% intervals and two-sided Fisher tests between the arms.
    [R, lab] = load_runs(odir);
    arms = {'control', 'off'}; T = struct();
    for a = 1:2
        q = find(strcmp(lab(:, 2), arms{a})); [~, o] = sort(cell2mat(lab(q, 1))); q = q(o);
        H = cellfun(@(r) r.hc, R(q), 'UniformOutput', false);
        cg = cellfun(@(h) first_or_inf(h < 0.5), H); fin = cellfun(@(h) h(end), H);
        % a: collapsed; b: collapsed and below 0.95 at generation 500; c: below 0.95 at 500
        T.(arms{a}) = struct('k', cell2mat(lab(q, 1)), 'H', {H}, 'cg', cg, 'fin', fin, ...
            'a', isfinite(cg), 'b', isfinite(cg) & fin < 0.95, 'c', fin < 0.95);
    end
    n = numel(T.control.k); assert(isequal(T.control.k, T.off.k) && n == numel(br));
    rep = {sprintf('\n[COUNTS] %d populations per arm, branch at first > 0.90, run to generation 500\n', n)};
    rep{end + 1} = sprintf('  %-8s %-34s %-34s %-34s\n', 'arm', '(a) collapse', '(b) lost, not recovered', '(c) final < 0.95');
    for a = 1:2
        S = T.(arms{a}); line = sprintf('  %-8s', arms{a});
        for f = {'a', 'b', 'c'}
            x = sum(S.(f{1})); [~, ci] = binofit(x, n);
            line = [line sprintf(' %3d / %d [%.3f, %.3f]        ', x, n, ci(1), ci(2))]; %#ok<AGROW>
        end
        rep{end + 1} = sprintf('%s\n', line); %#ok<AGROW>
    end
    for f = {'a', 'b', 'c'}
        x1 = sum(T.control.(f{1})); x2 = sum(T.off.(f{1}));
        [~, p] = fishertest([x1, n - x1; x2, n - x2]);
        rep{end + 1} = sprintf('  (%s) control %d vs off %d: delta %+d, two-sided Fisher p = %.3g\n', f{1}, x1, x2, x2 - x1, p); %#ok<AGROW>
    end
    rep{end + 1} = sprintf('\n[COLLAPSE WINDOWS] continuation generation of first drop below 0.5 (also reported)\n');
    for a = 1:2
        cg = T.(arms{a}).cg;
        rep{end + 1} = sprintf('  %-8s 1-30: %d | 31-150: %d | 151 on: %d | median final %.3f\n', arms{a}, sum(cg <= 30), ...
            sum(cg > 30 & cg <= 150), sum(cg > 150 & isfinite(cg)), median(T.(arms{a}).fin)); %#ok<AGROW>
    end
    rep = [rep, stored_window(d, br), start_populations_check(T, pub)];
    save(fullfile(odir, 'report.mat'), 'T');
    O = T;
end


function rep = stored_window(d, br)
% The same counts on the stored original runs over generations t0 to 500, a check on control.
    a = 0; b = 0; c = 0;
    for j = 1:numel(br)
        h = double(d.results{br(j).k}.historyCooperation(:))'; w = h(br(j).t0:500);
        col = any(w < 0.5); a = a + col; b = b + (col && w(end) < 0.95); c = c + (w(end) < 0.95);
    end
    rep = {sprintf('\n[STORED RUNS, same window t0..500, cross-check on control] (a) %d | (b) %d | (c) %d of %d\n', a, b, c, numel(br))};
end


function rep = start_populations_check(T, pub)
% The first 150 continued generations of control and off against run_refusal_start_populations's 'main' runs
% with g = 2 (its arms 'control', mutation on, and 'G', mutation rate 0), which use the same seeds.
    f = dir(fullfile(pub, 'refusal_start_populations', 'main_g2_reps*.mat')); P = {}; pl = {};
    if isempty(f)
        error('run_mutation_off:input', 'no %s: run run_refusal_start_populations(''main'', 2, reps) first', ...
            fullfile(pub, 'refusal_start_populations', 'main_g2_reps*.mat'));
    end
    for i = 1:numel(f), S = load(fullfile(f(i).folder, f(i).name)); P = [P; S.R]; pl = [pl; S.lab]; end %#ok<AGROW>
    map = {'control', 'control'; 'off', 'G'}; rep = {};
    for a = 1:2
        S = T.(map{a, 1}); same = 0; m = 0;
        for j = 1:numel(S.k)
            q = find(cell2mat(pl(:, 1)) == S.k(j) & strcmp(pl(:, 2), map{a, 2}));
            if isempty(q), continue; end
            hp = P{q}.hc; L = min(numel(hp), numel(S.H{j})); m = m + 1;
            same = same + isequal(hp(1:L), S.H{j}(1:L));
        end
        rep{end + 1} = sprintf('[REPRODUCTION] %s vs run_refusal_start_populations(''main'', 2) arm %s: first 150 generations identical in %d of %d\n', map{a, 1}, map{a, 2}, same, m); %#ok<AGROW>
    end
end


function [R, lab] = load_runs(odir)
    f = dir(fullfile(odir, 'runs_*.mat')); R = {}; lab = {};
    for i = 1:numel(f), S = load(fullfile(odir, f(i).name)); R = [R; S.R]; lab = [lab; S.lab]; end %#ok<AGROW>
end


function t = first_or_inf(v)
    t = find(v, 1); if isempty(t), t = inf; end
end


function open_pool(nw)
    if nargin < 1, nw = 2; end
    maxNumCompThreads(1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', nw); end
    wait(parfevalOnAll(pool, @maxNumCompThreads, 0, 1));
end


% Third arm 'refuse_less_conditional': mutation on, but a mutation that lowers the network's gap
% (mutate_with_refusal with 'REFUSE_LESS_CONDITIONAL') is refused and redrawn. Same start points,
% seeds and horizon as the other arms.
% The refused mutant would be less conditional than its parent. No paper item reports this arm, but 'report3'
% needs its runs
% to write refuse_less_conditional/report.mat, the file fig5_refusal_curve reads.
function r = r3_one_run(cfg)
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


function [O, rep] = r3_timing_run(base, br, rdir)
    open_pool(4); jobs = arrayfun(@(j) run_cfg(base, br(j), 'refuse_less_conditional', 501 - br(j).t0), 1:4, 'UniformOutput', false);
    R = cell(4, 1); t = tic;
    parfor j = 1:4, R{j} = r3_one_run(jobs{j}); end
    O.sec = toc(t); gens = arrayfun(@(j) 501 - br(j).t0, 1:4);
    lab = [num2cell([br(1:4).k]'), repmat({'refuse_less_conditional'}, 4, 1), num2cell([br(1:4).t0]')];
    save(fullfile(rdir, 'rruns_001-004.mat'), 'R', 'lab');
    rest = sum(501 - [br(5:end).t0]); perGen = max(cellfun(@(r) r.sec, R)) / max(gens);
    rep = {sprintf('[TIMING RUN] reps 1-4 (%s generations) on 4 workers: wall %.1f min; per run %s s; finals %s; refused first draws %s\n', ...
        mat2str(gens), O.sec / 60, mat2str(round(cellfun(@(r) r.sec, R))'), mat2str(cellfun(@(r) r.hc(end), R)', 3), ...
        mat2str(cellfun(@(r) r.log.refFirst / r.log.mut, R)', 3))};
    rep{end + 1} = sprintf('[TIMING RUN] ~%.3f s per generation per worker; reps 5-100 (%d generations) on 4 workers: ~%.0f min\n', ...
        perGen, rest, perGen * rest / 4 / 60);
end


function [O, rep] = r3_main(base, br, rdir, chunk)
    open_pool(4); tAll = tic; rep = {};
    for c0 = 5:chunk:numel(br)
        idx = c0:min(c0 + chunk - 1, numel(br));
        fn = fullfile(rdir, sprintf('rruns_%03d-%03d.mat', br(idx(1)).k, br(idx(end)).k));
        if exist(fn, 'file'), rep{end + 1} = sprintf('skip %s (exists)\n', fn); continue; end %#ok<AGROW>
        jobs = arrayfun(@(j) run_cfg(base, br(j), 'refuse_less_conditional', 501 - br(j).t0), idx, 'UniformOutput', false);
        lab = [num2cell([br(idx).k]'), repmat({'refuse_less_conditional'}, numel(idx), 1), num2cell([br(idx).t0]')];
        R = cell(numel(jobs), 1); t = tic;
        parfor j = 1:numel(jobs), R{j} = r3_one_run(jobs{j}); end
        for j = 1:numel(jobs)
            h = R{j}.hc; cg = first_or_inf(h < 0.5);
            rep{end + 1} = sprintf('rep %3d refuse_less_conditional t0 %3d gens %3d final %.3f collapse at %s refused %.3f (%.0f s)\n', lab{j, 1}, ...
                lab{j, 3}, numel(h), h(end), num2str(cg), R{j}.log.refFirst / R{j}.log.mut, R{j}.sec); %#ok<AGROW>
        end
        save(fn, 'R', 'lab');
        rep{end + 1} = sprintf('chunk %s: %.1f min (total %.1f min)\n', fn, toc(t) / 60, toc(tAll) / 60); %#ok<AGROW>
    end
    O = [];
end


function s = pf(ok), if ok, s = 'PASS'; else, s = 'FAIL'; end, end


function [O, rep] = r3_report(d, br, odir, rdir)
% Three-arm counts: stored control and off (mutation_off/report.mat) against the refuse_less_conditional runs.
    S = load(fullfile(odir, 'report.mat')); T = S.T;
    [R, lab] = load_runs3(rdir); k = cell2mat(lab(:, 1)); [k, o] = sort(k); R = R(o);
    assert(isequal(k(:), T.control.k(:)) && numel(k) == numel(br), 'refuse_less_conditional runs incomplete');
    H = cellfun(@(r) r.hc, R, 'UniformOutput', false); cg = cellfun(@(h) first_or_inf(h < 0.5), H);
    fin = cellfun(@(h) h(end), H);
    T.refuse_less_conditional = struct('k', k, 'H', {H}, 'cg', cg, 'fin', fin, 'a', isfinite(cg), 'b', isfinite(cg) & fin < 0.95, 'c', fin < 0.95);
    arms = {'control', 'off', 'refuse_less_conditional'}; n = numel(k);
    for a = 1:3, T.(arms{a}).d = T.(arms{a}).cg <= 40; end
    f = {'a', 'b', 'c', 'd'}; nm = {'(a) collapse', '(b) collapsed, < 0.95 at 500', '(c) < 0.95 at 500, any', '(d) collapse within 40'};
    O.cnt = zeros(4, 3); O.p = zeros(4, 2);
    rep = {sprintf('\n[COUNTS] %d populations per arm (control, off: stored runs)\n', n)};
    for q = 1:4
        x = cellfun(@(a) sum(T.(a).(f{q})), arms); O.cnt(q, :) = x;
        O.p(q, :) = [fp(x(3), x(1), n), fp(x(3), x(2), n)]; line = sprintf('  %-30s', nm{q});
        for a = 1:3, [~, ci] = binofit(x(a), n); line = [line sprintf(' %-8s %3d [%.3f, %.3f]', arms{a}, x(a), ci(1), ci(2))]; end %#ok<AGROW>
        rep{end + 1} = sprintf('%s | p refuse_less_conditional vs on %.3g, vs off %.3g\n', line, O.p(q, 1), O.p(q, 2)); %#ok<AGROW>
    end
    rep{end + 1} = sprintf('\n[COLLAPSE WINDOWS] continuation generation of the first drop below 0.5\n');
    for a = 1:3
        c = T.(arms{a}).cg; fa = T.(arms{a}).fin; nc = ~isfinite(c);
        rep{end + 1} = sprintf('  %-8s 1-30: %d | 31-40: %d | 41-150: %d | 151 on: %d | recovered to >= 0.95: %d | not collapsed but < 0.95 (plateau): %d, their finals %s | median final %.3f\n', ...
            arms{a}, sum(c <= 30), sum(c > 30 & c <= 40), sum(c > 40 & c <= 150), sum(c > 150 & isfinite(c)), ...
            sum(~nc & fa >= 0.95), sum(nc & fa < 0.95), mat2str(sort(fa(nc & fa < 0.95))', 2), median(fa)); %#ok<AGROW>
    end
    late = O.cnt(1, :) - O.cnt(4, :);
    rep{end + 1} = sprintf('  collapses after 40: on %d, off %d, refuse_less_conditional %d; p refuse_less_conditional vs on %.3g, vs off %.3g\n', late, ...
        fp(late(3), late(1), n), fp(late(3), late(2), n));
    [O.ops, r1] = manip(R, br); [O.rescore, r2] = rescore(R, d.config); rep = [rep, r1, r2];
    save(fullfile(rdir, 'report.mat'), 'T', 'O');
end


function [o, rep] = manip(R, br)
    lg = cellfun(@(r) r.log, R); rf = [lg.refFirst] ./ [lg.mut];
    rep = {sprintf('\n[OPERATOR] mutations %d; refused first draws %d (%.3f; per run %.3f-%.3f, median %.3f); redraws %d; redraw-limit failures %d (%.4f of mutations)\n', ...
        sum([lg.mut]), sum([lg.refFirst]), sum([lg.refFirst]) / sum([lg.mut]), min(rf), max(rf), median(rf), ...
        sum([lg.redraws]), sum([lg.fail]), sum([lg.fail]) / sum([lg.mut]))};
    rep{end + 1} = sprintf('  operator calls with no population defector in the partner set: %d of %d\n', sum([lg.noPop]), sum([lg.calls]));
    g = [1 41 151]; sh = nan(numel(R), 4);
    for j = 1:numel(R)
        r = R{j};
        for q = 1:3, i = find(r.hg == g(q), 1); if ~isempty(i), sh(j, q) = r.defShare(i); end, end
        sh(j, 4) = r.defShare(end);
    end
    rep{end + 1} = sprintf('  defect-against-clone share (95%% genotypes) at continuation generations 1 / 41 / 151 / last: mean %s; populations holding any %s of %s\n', ...
        mat2str(mean(sh, 1, 'omitnan'), 3), mat2str(sum(sh > 0, 1)), mat2str(sum(~isnan(sh), 1)));
    o = struct('refFirst', rf, 'sh', sh, 'nBr', numel(br));
end


function [a, rep] = rescore(R, cfg)
% Accepted mutants from each run's sample, scored again with fresh caches: none may have a smaller
% gap than its parent.
    P = {}; Mu = {}; S = {};
    for j = 1:numel(R), P = [P, R{j}.sampleP]; Mu = [Mu, R{j}.sampleM]; S = [S, R{j}.sampleS]; end %#ok<AGROW>
    n = numel(Mu); v = zeros(n, 2); cfg.refusalMode = 'REFUSE_LESS_CONDITIONAL'; open_pool(4);
    parfor q = 1:n
        mutate_with_refusal('reset');
        v(q, :) = [mutate_with_refusal('gap', P{q}, S{q}, cfg), mutate_with_refusal('gap', Mu{q}, S{q}, cfg)];
    end
    a.n = n; a.down = sum(v(:, 2) < v(:, 1) - 1e-9); a.pass = a.down == 0;
    rep = {sprintf('[RESCORE] %d accepted reservoir mutants, fresh caches: smaller gap %d -> %s\n', n, a.down, pf(a.pass))};
end


function [R, lab] = load_runs3(rdir)
    f = dir(fullfile(rdir, 'rruns_*.mat')); R = {}; lab = {};
    for i = 1:numel(f), S = load(fullfile(rdir, f(i).name)); R = [R; S.R]; lab = [lab; S.lab]; end %#ok<AGROW>
end


function p = fp(a, b, n)
    [~, p] = fishertest([a, n - a; b, n - b]);
end
