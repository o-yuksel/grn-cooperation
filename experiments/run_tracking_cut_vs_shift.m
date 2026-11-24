function out = run_tracking_cut_vs_shift(mode, runName, opt)
% RUN_TRACKING_CUT_VS_SHIFT  With the target held still, does the cut in the mutation load follow
% the shift one mutation causes rather than the design? Measures the mean shift per mutation of
% each design at each mutation size, and runs design 2 (and design 1) at further sizes.
%
% Produces : Fig. 6; SI S5.1; Table S15 (shift_by_size.txt, shift_by_size_design1.txt and cut_vs_shift_summary.txt are also
%            read by tracking_load_rank_correlation)
% Inputs   : `mode`, the phase: 'x' (the shifts at every point; finds by bisection the size at
%              which design 2's mean shift matches design 1's), 'run' (one run, runName 'A', 'C',
%              'B' or 'Amean'), 'final' (table and comparisons), 'all' (default: 'x', the four
%              runs, 'final'), 'x3' (design 1 at the seven sizes of Fig. 6B; runName, if given, is
%              the publication folder to read); `opt`, optional, for 'run': repList (1:20, the
%              replicates) and nWorkers (3)
%            'x', 'x3', 'run': publication/tracking_runs/tracking_runs_N1000.mat
%              (run_tracking_task); 'run' also lowering_settings.mat (run_tracking_single_genotypes)
%            'final': tracking_load_regime/shift_by_size.mat and <run folder>/main.mat (this file; folders design2_sd0p001821, design1_sd0p060, design2_sd0p0063, design2_matched_shift for runs A, C, B, Amean), and
%              under publication/, tracking_cut_after_150/design1 and design2/main.mat
%              (tracking_cut_after_150), tracking_mutation_regime/size_10x and rare/main.mat
%              (run_tracking_mutation_regime)
% Outputs  : in publication/tracking_load_regime/: 'x' shift_by_size.mat, shift_by_size.txt; 'x3' shift_by_size_design1.txt;
%              'run' <run folder>/main.mat, main.txt, lines appended to cut_vs_shift_summary.txt; 'final'
%              the table and comparisons appended to cut_vs_shift_summary.txt (replacing an earlier 'final'
%              block)
% Run      : run_tracking_cut_vs_shift('all')   ('x' about 2 min; each run 22 to 47 min on 3
%              workers; 'final' about 30 s); run_tracking_cut_vs_shift('x3')   (about 1 min)
% Seed     : the mean shift per mutation of best genotype r: rng(20260917 + 1000 f + r, 'twister'),
%            f = 1 in design 1 and 3 in design 2 (run_tracking_task's seeds); the runs:
%            run_tracking_held_still's
%
% The mean shift per mutation is the mean distance outputs 3 and 4 move under one application of
% the mutation distribution (4000 draws, W33 and W44 held) at the probe input u = [0, 0.25, 0, 0],
% on the 30 generation-350 best genotypes with the target held still, at each point's own size;
% the axis is its median over genotypes (the mean printed beside it).

    if nargin < 1 || isempty(mode), mode = 'all'; end
    if nargin < 3 || isempty(opt), opt = struct(); end
    if ~isfield(opt, 'repList'), opt.repList = 1:20; end
    if ~isfield(opt, 'nWorkers'), opt.nWorkers = 3; end
    pub = fullfile(grn_data_file('publication'));
    od = fullfile(pub, 'tracking_load_regime');
    switch mode
        case 'x', out = measure_x(pub, od);
        case 'x3'
            if nargin >= 2 && ~isempty(runName), pub = runName; end
            out = measure_x3(pub, od);
        case 'run', out = one_run(runName, od, opt);
        case 'final', out = final(pub, od);
        case 'all'
            measure_x(pub, od);
            for r = {'A', 'C', 'B', 'Amean'}, one_run(r{1}, od, opt); end
            out = final(pub, od);
    end
end


% ===================== runs ========================================
function R = run_defs(od)
% Each run is run_tracking_held_still with the target held still: 150 generations with mutation
% after the split and 150 without, per-entry probability 0.08, replicates 1 to 20, tolerance 0.0106 (design 1,
% option designIndex 1) or 0.0110 (design 2, designIndex 3); alpha or ||x*|| is lowered to
% 0.70. A: design 2, sd 0.001821. C: design 1, sd 0.060. B: design 2, sd 0.0063. Amean: design 2,
% the sd phase 'x' found by bisection.
    R = struct('name', {'A', 'C', 'B', 'Amean'}, 'designIndex', {3, 1, 3, 3}, 'design', {2, 1, 2, 2}, ...
        'dir', {'design2_sd0p001821', 'design1_sd0p060', 'design2_sd0p0063', 'design2_matched_shift'}, ...
        'tol', {0.0110, 0.0106, 0.0110, 0.0110}, 'sd', {0.001821, 0.060, 0.0063, NaN});
    f = fullfile(od, 'shift_by_size.mat');
    if exist(f, 'file'), X = getfield(load(f), 'out'); R(4).sd = X.ameanSd; end
end


function S = one_run(name, od, opt)
    R = run_defs(od); k = find(strcmp({R.name}, name)); r = R(k);
    outDir = fullfile(od, r.dir);
    if exist(fullfile(outDir, 'main.mat'), 'file')
        fprintf('%s/main.mat exists: run %s skipped (never overwritten)\n', r.dir, name);
        S = []; return
    end
    assert(isfinite(r.sd), 'run %s has no mutSd (run mode ''x'' first)', name);
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    c = struct('designIndex', r.designIndex, 'tolerance', r.tol, 'nOff', 150, 'mutProb', 0.08, 'mutSd', r.sd, ...
        'repList', opt.repList, 'nWorkers', opt.nWorkers, 'outDir', outDir, 'pilot', false, 'm', 0.7);
    t0 = tic;
    run_tracking_held_still(c);
    S = append_run(name, r, outDir, od, toc(t0) / 60);
end


function S = append_run(name, r, outDir, od, minutes)
    o = getfield(load(fullfile(outDir, 'main.mat')), 'out'); S = run_stats(o);
    fid = fopen(fullfile(od, 'cut_vs_shift_summary.txt'), 'a'); pr = @(varargin) both(fid, varargin{:});
    pr('\n==== RUN %s (design %d, a = %.2f) | mutation distribution p %.3f sd %.6f | m %.2f | %d replicates | %.1f min (%s) ====\n', ...
        name, r.design, o.E.a, o.P.mutation_prob, o.P.mutation_sd, o.cfg.m, numel(o.cfg.repList), minutes, ...
        datestr(now, 'yyyy-mm-dd HH:MM'));
    pr('runner checks: check against the stored run %s [best genotype bestCheck %.3g, knobs solved A %d X %d] | window score %s | knob verification kept A %d/%d, X %d/%d | control check %s | AT BALANCE [C C2 A X] %s\n', ...
        pf(o.c1.pass), o.c1.bestCheck, o.c1.nOK(1), o.c1.nOK(2), pf(o.c0b.pass), sum(o.verif.rep(:, 1)), ...
        size(o.verif.rep, 1), sum(o.verif.rep(:, 2)), size(o.verif.rep, 1), pf(o.analysis.c1c.pass), ...
        mat2str(o.analysis.plateau.atBalance));
    pr('CONTROL mutation load (median) %.5f\n', S.c0);
    pr('ALPHA-DOWN cut %.1f%% [%.1f, %.1f] (n %d) | XNORM-DOWN cut %.1f%% [%.1f, %.1f] (n %d)\n', ...
        S.cutA, S.ciA, S.nA, S.cutX, S.ciX, S.nX);
    pr('ALPHA/XNORM effect ratio %.3f [%.3f, %.3f]\n', S.ratio, S.ratioCi);
    pr('CONTROL2 - CONTROL d %+.5f [%+.5f, %+.5f] = %+.1f%% of CONTROL\n', S.d2, S.ci2, S.pct2);
    fclose(fid);
end


% ===================== x-axis =======================================
function X = measure_x(pub, od)
    f = fullfile(od, 'shift_by_size.mat');
    if exist(f, 'file'), X = getfield(load(f), 'out'); fprintf('shift_by_size.mat exists: x not remeasured\n'); return; end
    txt = fullfile(od, 'shift_by_size.txt'); if exist(txt, 'file'), delete(txt); end
    fid = fopen(txt, 'w'); pr = @(varargin) both(fid, varargin{:});
    nThreads = maxNumCompThreads(1); restoreThreads = onCleanup(@() maxNumCompThreads(nThreads));
    t0 = tic;
    D = getfield(load(fullfile(pub, 'tracking_runs', 'tracking_runs_N1000.mat')), 'out');
    P = D.P; uProbe = [0; P.trajectory_radius; 0; 0]; nDraws = 4000; seed0 = 20260917;
    % The shifts of design f's 30 best genotypes with the target held still, at mutation size sd
    % and per-entry probability p
    meas = @(f, sd, p) measure(D.fixed.ThBest(f, :), P, sd, p, uProbe, nDraws, seed0, f);
    pr('SHIFT BY SIZE, %s | estimator fragility_fast (this file) | %d draws | %d best genotypes per design\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), nDraws, D.nReps);
    ref = meas(1, 0.020, 0.08); refMean = mean(ref);
    pr('design 1 reference at sd 0.020: mean %.5f, median %.5f\n', refMean, median(ref));

    % A-mean: the sd at which design 2's mean shift is within 1.5 times design 1's at sd 0.020.
    % The sd is divided by the ratio until the target is bracketed, then bisected on a log
    % scale; at most 7 tries
    pr('\nA-MEAN BISECTION (design 2 mean shift / design 1 mean %.5f, target within 1.5x):\n', refMean);
    sd = 0.000562; lo = NaN; hi = NaN; steps = zeros(0, 3);
    for it = 0:6
        v = meas(3, sd, 0.08); rM = mean(v) / refMean; ok = rM >= 1/1.5 && rM <= 1.5;
        steps(end+1, :) = [sd, mean(v), rM]; %#ok<AGROW>
        pr('step %d: sd %.6f -> mean %.5f (median %.5f), ratio %.3f -> %s\n', it, sd, mean(v), median(v), rM, ...
            ternary(ok, 'WITHIN 1.5', 'OUTSIDE'));
        if ok, break; end
        if rM > 1, hi = sd; else, lo = sd; end
        if isnan(lo) || isnan(hi), sd = sd / rM; else, sd = sqrt(lo * hi); end
    end
    ameanSd = ternary(ok, sd, NaN);
    pr('A-mean mutSd: %s\n', ternary(ok, sprintf('%.6f', sd), 'NONE (bisection failed)'));

    % Every point: label, design (f = 1 or 3), sd and per-entry probability p (p 0.008: mutations
    % ten times rarer)
    C = struct('label', {'d1 sd0.020', 'd1 sd0.060', 'd1 sd0.200', 'd1 rare p0.008 sd0.020', ...
        'd2 sd0.020', 'd2 sd0.0063', 'd2 sd0.001821', 'd2 A-mean'}, 'f', {1, 1, 1, 1, 3, 3, 3, 3}, ...
        'sd', {0.020, 0.060, 0.200, 0.020, 0.020, 0.0063, 0.001821, ameanSd}, ...
        'p', {0.08, 0.08, 0.08, 0.008, 0.08, 0.08, 0.08, 0.08});
    pr('\nCONTEXTS (shift per mutation over %d best genotypes):\n', D.nReps);
    pr('%-24s %5s %9s %6s %10s %10s %10s %10s\n', 'context', 'index', 'sd', 'p', 'median', 'mean', 'IQR lo', 'IQR hi');
    for i = 1:numel(C)
        if ~isfinite(C(i).sd), C(i).med = NaN; C(i).mean = NaN; continue; end
        v = meas(C(i).f, C(i).sd, C(i).p); C(i).med = median(v); C(i).mean = mean(v); C(i).v = v;
        pr('%-24s %5d %9.6f %6.3f %10.5f %10.5f %10.5f %10.5f\n', C(i).label, C(i).f, C(i).sd, C(i).p, ...
            C(i).med, C(i).mean, prctile(v, 25), prctile(v, 75));
    end
    pr('elapsed %.1f min\n', toc(t0) / 60);
    fclose(fid);
    X = struct('ref', ref, 'steps', steps, 'ameanSd', ameanSd, 'ctx', C, 'finishedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    out = X; save(f, 'out');
end


function X = measure_x3(pub, od)
% Stops with an error unless the three sizes already in shift_by_size.txt (sd 0.020, 0.060, 0.200)
% give the same median and mean to 5 decimals
    old = fileread(fullfile(od, 'shift_by_size.txt'));
    fid = fopen(fullfile(od, 'shift_by_size_design1.txt'), 'w'); pr = @(varargin) both(fid, varargin{:});
    nThreads = maxNumCompThreads(1); restoreThreads = onCleanup(@() maxNumCompThreads(nThreads));
    t0 = tic;
    D = getfield(load(fullfile(pub, 'tracking_runs', 'tracking_runs_N1000.mat')), 'out');
    P = D.P; uProbe = [0; P.trajectory_radius; 0; 0]; nDraws = 4000; seed0 = 20260917;
    pr('SHIFT BY SIZE, DESIGN 1, %s | estimator and seeds of mode ''x'' (shift_by_size.txt) | %d draws | %d best genotypes | design 1, p 0.080\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), nDraws, D.nReps);
    sds = [0.006 0.010 0.020 0.040 0.060 0.100 0.200];
    pr('%-12s %9s %10s %10s %10s %10s\n', 'size', 'sd', 'median', 'mean', 'IQR lo', 'IQR hi');
    X = struct('sd', num2cell(sds), 'med', NaN, 'mean', NaN); ok = true; nOld = 0;
    for i = 1:numel(sds)
        v = measure(D.fixed.ThBest(1, :), P, sds(i), 0.08, uProbe, nDraws, seed0, 1);
        X(i).med = median(v); X(i).mean = mean(v);
        pr('%-12s %9.6f %10.5f %10.5f %10.5f %10.5f\n', sprintf('%.4gx', sds(i) / 0.020), sds(i), X(i).med, X(i).mean, ...
            prctile(v, 25), prctile(v, 75));
        % A size already in shift_by_size.txt (its design-1 row at p 0.080) must give the same median
        % and mean within 5e-6
        tok = regexp(old, ['d1 sd' sprintf('%.3f', sds(i)) '\s+1\s+[\d.]+\s+0\.080\s+([\d.]+)\s+([\d.]+)'], 'tokens', 'once');
        if ~isempty(tok)
            nOld = nOld + 1;
            ok = ok && all(abs([X(i).med X(i).mean] - str2double(tok)) < 5e-6);
        end
    end
    ok = ok && nOld == 3;
    pr('CONTROL sd 0.020 / 0.060 / 0.200 reproduce shift_by_size.txt median and mean (5 dp, %d of 3 found): %s\n', nOld, pf(ok));
    pr('elapsed %.1f min\n', toc(t0) / 60);
    fclose(fid);
    if ~ok, error('run_tracking_cut_vs_shift:control', 'mode x3 does not reproduce shift_by_size.txt'); end
end


function b = measure(ThB, P, sd, p, uProbe, nDraws, seed0, f)
% The shift per mutation of each genotype in ThB, its draws seeded from design f and replicate r
    P.mutation_sd = sd; P.mutation_prob = p; b = nan(1, numel(ThB));
    for r = 1:numel(ThB)
        rng(seed0 + 1000*f + r, 'twister');
        b(r) = fragility_fast(ThB{r}, uProbe, P, nDraws);
    end
end


% ===================== final ========================================
function F = final(pub, od)
    X = getfield(load(fullfile(od, 'shift_by_size.mat')), 'out'); C = X.ctx; lab = {C.label};
    src = {'d1 sd0.020', fullfile(pub, 'tracking_cut_after_150', 'design1', 'main.mat'), 'reused tail-150 design 1';
           'd1 sd0.060', fullfile(od, 'design1_sd0p060', 'main.mat'), 'run C';
           'd1 sd0.200', fullfile(pub, 'tracking_mutation_regime', 'size_10x', 'main.mat'), 'reused tracking_mutation_regime/size_10x';
           'd2 sd0.020', fullfile(pub, 'tracking_cut_after_150', 'design2', 'main.mat'), 'reused tail-150 design 2';
           'd2 sd0.0063', fullfile(od, 'design2_sd0p0063', 'main.mat'), 'run B';
           'd2 sd0.001821', fullfile(od, 'design2_sd0p001821', 'main.mat'), 'run A';
           'd2 A-mean', fullfile(od, 'design2_matched_shift', 'main.mat'), 'run A-mean';
           'd1 rare p0.008 sd0.020', fullfile(pub, 'tracking_mutation_regime', 'rate_0p1x', 'main.mat'), 'reused tracking_mutation_regime/rate_0p1x (separate marker)'};
    n = size(src, 1); Pt = struct('label', src(:, 1), 'src', src(:, 3));
    for i = 1:n
        j = find(strcmp(lab, src{i, 1}));
        Pt(i).sd = C(j).sd; Pt(i).p = C(j).p; Pt(i).f = C(j).f; Pt(i).med = C(j).med; Pt(i).mean = C(j).mean;
        Pt(i).ok = exist(src{i, 2}, 'file') == 2;
        if Pt(i).ok, Pt(i).S = run_stats(getfield(load(src{i, 2}), 'out')); else, Pt(i).S = []; end
    end
    drop_final_block(fullfile(od, 'cut_vs_shift_summary.txt'), '================ CUT VS SHIFT FINAL');
    fid = fopen(fullfile(od, 'cut_vs_shift_summary.txt'), 'a'); pr = @(varargin) both(fid, varargin{:});
    pr('\n\n================ CUT VS SHIFT FINAL (%s) ================\n', datestr(now, 'yyyy-mm-dd HH:MM'));
    pr('x = shift per mutation over the 30 gen-350 best genotypes (shift_by_size.txt); cut = ALPHA-DOWN, %% of CONTROL [paired-bootstrap 95%%]\n');
    pr('%-24s %5s %9s %6s %10s %10s %22s %22s  %s\n', 'context', 'index', 'mutSd', 'p', 'median x', 'mean x', ...
        'ALPHA cut [CI]', 'XNORM cut [CI]', 'source');
    for i = 1:n
        if Pt(i).ok
            s = Pt(i).S;
            ca = sprintf('%.1f [%.1f, %.1f]', s.cutA, s.ciA); cx = sprintf('%.1f [%.1f, %.1f]', s.cutX, s.ciX);
        else, ca = 'missing'; cx = 'missing';
        end
        pr('%-24s %5d %9.6f %6.3f %10.5f %10.5f %22s %22s  %s\n', Pt(i).label, Pt(i).f, Pt(i).sd, Pt(i).p, ...
            Pt(i).med, Pt(i).mean, ca, cx, Pt(i).src);
    end
    F = rules(Pt, pr);
    fclose(fid);
end


function F = rules(Pt, pr)
    get = @(l) Pt(strcmp({Pt.label}, l));   % local name; shadows the graphics get() here only
    pr('\n---- COMPARISONS ----\n');
    d1 = get('d1 sd0.020'); A = get('d2 sd0.001821'); Am = get('d2 A-mean');
    % Run A (design 2 at a median shift close to design 1's) against design 1
    if A.ok
        pr('run A cut %.1f%% [%.1f, %.1f]; design 1 [%.1f, %.1f]\n', A.S.cutA, A.S.ciA, d1.S.ciA);
    end
    % Spearman correlation of the ALPHA cut with the median shift over the seven size points; the
    % same with the mean shift is printed beside it
    sz = {'d1 sd0.020', 'd1 sd0.060', 'd1 sd0.200', 'd2 sd0.020', 'd2 sd0.0063', 'd2 sd0.001821', 'd2 A-mean'};
    Q = cellfun(@(l) get(l), sz, 'UniformOutput', false); Q = [Q{:}]; okq = [Q.ok];
    x = [Q(okq).med]'; y = arrayfun(@(q) q.S.cutA, Q(okq))';
    rho = corr(x, y, 'type', 'Spearman');
    F.rho = rho;
    pr('Spearman(median x, cut) over %d of 7 size points = %+.3f%s\n', sum(okq), rho, ...
        ternary(all(okq), '', ' (INCOMPLETE: missing points)'));
    xm = [Q(okq).mean]'; pr('        Spearman(mean x, cut) = %+.3f\n', corr(xm, y, 'type', 'Spearman'));
    % The choice of average: the A-mean run against run A
    if A.ok && Am.ok
        a = A.S.ciA; b = Am.S.ciA;
        pr('A-mean cut %.1f%% [%.1f, %.1f] (sd %.6f) vs A %.1f%% [%.1f, %.1f] (sd %.6f); median-matched %.1f%%, mean-matched %.1f%%\n', ...
            Am.S.cutA, b, Am.sd, A.S.cutA, a, A.sd, A.S.cutA, Am.S.cutA);
    end
    % The ALPHA/XNORM ratio of the cuts in each run
    for l = {'d2 sd0.001821', 'A'; 'd2 sd0.0063', 'B'; 'd1 sd0.060', 'C'}'
        q = get(l{1});
        if ~q.ok, continue; end
        pr('run %s: ALPHA/XNORM ratio %.3f [%.3f, %.3f]\n', l{2}, q.S.ratio, q.S.ratioCi);
    end
    % The second control against the control, in each run of this file
    for l = {'d2 sd0.001821', 'A'; 'd1 sd0.060', 'C'; 'd2 sd0.0063', 'B'; 'd2 A-mean', 'A-mean'}'
        q = get(l{1});
        if ~q.ok, continue; end
        pr('run %s: CONTROL2 - CONTROL d %+.5f [%+.5f, %+.5f] = %+.1f%% of CONTROL %.5f\n', ...
            l{2}, q.S.d2, q.S.ci2, q.S.pct2, q.S.c0);
    end
end


function both(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
