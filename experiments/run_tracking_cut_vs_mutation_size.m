function out = run_tracking_cut_vs_mutation_size(mode, runName, opt)
% RUN_TRACKING_CUT_VS_MUTATION_SIZE  With the target held still in design 1, how does the cut in
% the mutation load from lowering alpha or ||x*|| to 0.70 change with mutation size? Runs four
% further sizes and prints the cuts at all seven (sd 0.006 to 0.200).
%
% Produces : Methods, "Non-social tracking task"; main text, "Responsiveness carries a mutation
%            load in a non-social trait"; Fig. 6B; SI S5.1; Table S15 (cut_vs_size_summary.txt is also
%            read by fig6_tracking_load and tracking_load_rank_correlation)
% Inputs   : `mode`, the phase: 'run' (one run, runName 'E', 'F', 'G' or 'H'), 'final' (the
%              table), 'all' (default: runs E to H, then 'final'); `opt`, optional, for 'run':
%              repList (1:20, the replicates) and nWorkers (3)
%            'run': those of run_tracking_held_still (tracking_runs_N1000.mat, lowering_settings.mat)
%            'final', under publication/: tracking_load_regime/design1_sd0p006 to design1_sd0p100/main.mat
%              ('run'), tracking_load_regime/design1_sd0p060/main.mat (run_tracking_cut_vs_shift),
%              tracking_cut_after_150/design1/main.mat (tracking_cut_after_150) and
%              tracking_mutation_regime/size_10x/main.mat (run_tracking_mutation_regime)
% Outputs  : 'run': publication/tracking_load_regime/<run folder>/main.mat and main.txt (design1_sd0p006, design1_sd0p010, design1_sd0p040, design1_sd0p100 for runs E, F, G, H), and the
%              run's lines appended to tracking_load_regime/cut_vs_size_summary.txt
%            'final': the seven-size table and comparisons appended to cut_vs_size_summary.txt (replacing
%              an earlier 'final' block)
% Run      : run_tracking_cut_vs_mutation_size('all')   (runs E to H, about 30 min each on 3
%              workers, then 'final'); run_tracking_cut_vs_mutation_size('final')   (seconds)
% Seed     : run_tracking_held_still's (replicate r continues with rng(20261001 + r, 'twister'));
%            'final' draws nothing
%
% Each run is run_tracking_held_still in design 1 (designIndex 1, a = 0.70): 150 generations with
% mutation after the split and 150 without, per-entry probability 0.08, replicates 1 to 20,
% tolerance 0.0106. Only the mutation size after the split changes: E sd 0.006 (0.3 times the
% standard 0.020), F 0.010 (0.5), G 0.040 (2), H 0.100 (5); 'final' adds 1, 3 and 10 times from
% the runs above. A run whose main.mat exists is skipped, never overwritten.

    if nargin < 1 || isempty(mode), mode = 'all'; end
    if nargin < 3 || isempty(opt), opt = struct(); end
    if ~isfield(opt, 'repList'), opt.repList = 1:20; end
    if ~isfield(opt, 'nWorkers'), opt.nWorkers = 3; end
    pub = fullfile(grn_data_file('publication'));
    od = fullfile(pub, 'tracking_load_regime');
    switch mode
        case 'run', out = one_run(runName, od, opt);
        case 'final', out = final(pub, od);
        case 'all'
            for r = {'E', 'F', 'G', 'H'}, one_run(r{1}, od, opt); end
            out = final(pub, od);
    end
end


function R = run_defs()
    R = struct('name', {'E', 'F', 'G', 'H'}, 'sd', {0.006, 0.010, 0.040, 0.100}, ...
        'dir', {'design1_sd0p006', 'design1_sd0p010', 'design1_sd0p040', 'design1_sd0p100'});
end


function S = one_run(name, od, opt)
    R = run_defs(); r = R(strcmp({R.name}, name));
    outDir = fullfile(od, r.dir);
    if exist(fullfile(outDir, 'main.mat'), 'file')
        fprintf('%s/main.mat exists: run %s skipped (never overwritten)\n', r.dir, name);
        S = []; return
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    c = struct('designIndex', 1, 'tolerance', 0.0106, 'nOff', 150, 'mutProb', 0.08, 'mutSd', r.sd, ...
        'repList', opt.repList, 'nWorkers', opt.nWorkers, 'outDir', outDir, 'pilot', false, 'm', 0.7);
    t0 = tic;
    run_tracking_held_still(c);
    S = append_run(name, outDir, od, toc(t0) / 60);
end


function S = append_run(name, outDir, od, minutes)
    o = getfield(load(fullfile(outDir, 'main.mat')), 'out'); S = run_stats(o);
    fid = fopen(fullfile(od, 'cut_vs_size_summary.txt'), 'a'); pr = @(varargin) both(fid, varargin{:});
    pr('\n==== RUN %s (design 1, a = %.2f) | mutation distribution p %.3f sd %.6f | m %.2f | %d replicates | %.1f min (%s) ====\n', ...
        name, o.E.a, o.P.mutation_prob, o.P.mutation_sd, o.cfg.m, numel(o.cfg.repList), minutes, ...
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


function F = final(pub, od)
% The seven sizes of Fig. 6B: label, size as a multiple of the standard sd 0.020, the run's file
    L = {'0.3x (E)', 0.3, fullfile(od, 'design1_sd0p006', 'main.mat');
         '0.5x (F)', 0.5, fullfile(od, 'design1_sd0p010', 'main.mat');
         '1x (d1 sd0.020)', 1, fullfile(pub, 'tracking_cut_after_150', 'design1', 'main.mat');
         '2x (G)', 2, fullfile(od, 'design1_sd0p040', 'main.mat');
         '3x (d1 sd0.060)', 3, fullfile(od, 'design1_sd0p060', 'main.mat');
         '5x (H)', 5, fullfile(od, 'design1_sd0p100', 'main.mat');
         '10x (d1 sd0.200)', 10, fullfile(pub, 'tracking_mutation_regime', 'size_10x', 'main.mat')};
    n = size(L, 1); S = cell(n, 1);
    for i = 1:n, S{i} = run_stats(getfield(load(L{i, 3}), 'out')); end
    drop_final_block(fullfile(od, 'cut_vs_size_summary.txt'), '================ CUT VS SIZE FINAL');
    fid = fopen(fullfile(od, 'cut_vs_size_summary.txt'), 'a'); pr = @(varargin) both(fid, varargin{:});
    pr('\n\n================ CUT VS SIZE FINAL (%s) ================\n', datestr(now, 'yyyy-mm-dd HH:MM'));
    pr('design 1, target held still, responsiveness x0.70; cut = %% of CONTROL load [paired-bootstrap 95%%]\n');
    pr('%-18s %7s %22s %22s %4s %4s %11s %9s\n', 'size', 'x std', 'ALPHA cut [CI]', 'XNORM cut [CI]', 'nA', 'nX', 'ALPHA/XNORM', 'ctrl load');
    for i = 1:n
        s = S{i};
        pr('%-18s %7.1f %22s %22s %4d %4d %11.3f %9.5f\n', L{i, 1}, L{i, 2}, ...
            sprintf('%.1f [%.1f, %.1f]', s.cutA, s.ciA), sprintf('%.1f [%.1f, %.1f]', s.cutX, s.ciX), ...
            s.nA, s.nX, s.ratio, s.c0);
    end
    x = cell2mat(L(:, 2)); y = cellfun(@(s) s.cutA, S); ci = cell2mat(cellfun(@(s) s.ciA, S, 'UniformOutput', false));
    pr('\n---- COMPARISONS ----\n');
    % Spearman correlation of the ALPHA cut with size; a reversal is an adjacent pair where the
    % larger size has the larger cut and the two intervals are disjoint
    rho = corr(x, y, 'type', 'Spearman');
    F.rho = rho;
    pr('Spearman(size, ALPHA-DOWN cut) over 7 sizes = %+.3f\n', rho);
    nrev = 0;
    for i = 1:n-1
        if y(i+1) > y(i) && ci(i+1, 1) > ci(i, 2)
            nrev = nrev + 1;
            pr('        REVERSAL: %s cut %.1f [%.1f, %.1f] -> %s cut %.1f [%.1f, %.1f] (larger size, larger cut, intervals disjoint)\n', ...
                L{i, 1}, y(i), ci(i, :), L{i+1, 1}, y(i+1), ci(i+1, :));
        end
    end
    if nrev == 0, pr('        reversals (adjacent pair, larger size, larger cut, disjoint intervals): none\n'); end
    % The second control against the control in each new run (E to H)
    for i = 1:n
        nm = regexp(L{i, 1}, '\(([EFGH])\)', 'tokens', 'once');
        if ~isempty(nm)
            pr('run %s: CONTROL2 - CONTROL d %+.5f [%+.5f, %+.5f] = %+.1f%% of CONTROL %.5f\n', ...
                nm{1}, S{i}.d2, S{i}.ci2, S{i}.pct2, S{i}.c0);
        end
    end
    % The cut's interval at 0.3 times the standard size against that at the standard
    i03 = 1; i1 = 3;
    pr('        cut at 0.3x %.1f [%.1f, %.1f] against 1x %.1f [%.1f, %.1f]\n', y(i03), ci(i03, :), y(i1), ci(i1, :));
    fclose(fid);
end


function both(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
