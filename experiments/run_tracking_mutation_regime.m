function S = run_tracking_mutation_regime(mode, opt)
% RUN_TRACKING_MUTATION_REGIME  With the target held still in design 1, does lowering
% responsiveness still cut the mutation load when each mutation is ten times larger, or when
% mutations are ten times rarer? Reruns the population test of the load both ways.
%
% Produces : Table S14 (the rows with mutations ten times larger and ten times rarer); Table S15
%            (the sd 0.200 row); Fig. 6B (the point at ten times the standard size); SI S5.1;
%            Methods, "Non-social tracking task" (main.mat is read by fig6_tracking_load,
%            run_tracking_cut_vs_shift and run_tracking_cut_vs_mutation_size)
% Inputs   : `mode`, the phase: 'all' (default: both runs, then the summary), 'main' (both runs),
%              'summary' (reads the two main.mat and prints); `opt`, optional: repList (1:20, the
%              replicates of the runs) and nWorkers (2)
%            the inputs of run_tracking_held_still (publication/tracking_runs/
%              tracking_runs_N1000.mat and publication/tracking_single_genotypes/lowering_settings.mat)
% Outputs  : publication/tracking_mutation_regime/size_10x/main.mat and main.txt   (sd 0.200)
%            publication/tracking_mutation_regime/rate_0p1x/main.mat and main.txt  (probability 0.008)
%            publication/tracking_mutation_regime/summary.txt   (what 'summary' prints)
% Run      : run_tracking_mutation_regime('all')   (about 40 min for big and 56 min for rare, on
%              2 workers, then the summary)
% Seed     : those of run_tracking_held_still (replicate r continues with rng(20261001 + r,
%            'twister')); the summary draws nothing
%
% Each run is run_tracking_held_still in design 1 (a = 0.70), replicates 1 to 20, 150 generations
% with mutation and 150 without, with one change after the split: the per-entry probability or
% the size of a mutation. Nothing before the split reads them, so the split and its (b, delta)
% are those of the standard run. The runner's check that the control's mean fitness reproduces
% the stored one is expected to fail here, since mutation changed, and is no check of these runs.

    if nargin < 1 || isempty(mode), mode = 'all'; end
    if nargin < 2 || isempty(opt), opt = struct(); end
    if ~isfield(opt, 'repList'), opt.repList = 1:20; end
    if ~isfield(opt, 'nWorkers'), opt.nWorkers = 2; end
    outRoot = fullfile(grn_data_file('publication'), 'tracking_mutation_regime');
    A = struct('name', {'size_10x', 'rate_0p1x'}, 'mutProb', {0.08, 0.008}, 'mutSd', {0.200, 0.020});
    S = struct();
    if any(strcmp(mode, {'main', 'all'}))
        for k = 1:numel(A)
            c = struct('designIndex', 1, 'tolerance', 0.0106, 'outDir', fullfile(outRoot, A(k).name), ...
                'repList', opt.repList, 'nWorkers', opt.nWorkers, 'mutProb', A(k).mutProb, ...
                'mutSd', A(k).mutSd, 'nOff', 150);
            run_tracking_held_still(c);
        end
    end
    if any(strcmp(mode, {'summary', 'all'}))
        S = summary(A, outRoot);
    end
end


function S = summary(A, outRoot)
% SUMMARY  Per run: the mutation load of the control and of the arm with alpha lowered, their
% paired difference and the ratio r of the two loads.
    f =fullfile(outRoot, 'summary.txt'); if exist(f, 'file'), delete(f); end
    fid = fopen(f, 'a'); pr = @(varargin) both(fid, varargin{:});
    pr('LOAD-REGIME CHECK, %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
    for k = 1:numel(A)
        o = getfield(load(fullfile(outRoot, A(k).name, 'main.mat')), 'out');
        Y = o.analysis.Y; T = o.analysis.table; keep = o.verif.rep;
        ok = isfinite(Y(:, 3, 1)) & isfinite(Y(:, 1, 1)) & keep(:, 1);
        c0 = median(Y(:, 1, 1), 'omitnan'); a0 = median(Y(ok, 3, 1));
        rMed = a0 / c0; rPer = median(Y(ok, 3, 1) ./ Y(ok, 1, 1));
        ciA = T(1, 3).ci; ciN = T(1, 2).ci;
        fitC = median(Y(:, 1, 6), 'omitnan'); bal = o.analysis.plateau.atBalance;
        pr('\n==== %s (p %.3f, sd %.3f; mutation distribution as run: p %.3f sd %.3f) ====\n', upper(A(k).name), ...
            A(k).mutProb, A(k).mutSd, o.P.mutation_prob, o.P.mutation_sd);
        pr('mutation load: CONTROL %.4f | ALPHA-DOWN %.4f (n %d) | d %+.4f [%+.4f, %+.4f] = %+.1f%% of CONTROL [%+.1f, %+.1f]\n', ...
            c0, a0, sum(ok), T(1, 3).d, ciA, 100*T(1, 3).d/c0, 100*ciA/c0);
        pr('r = ALPHA-DOWN/CONTROL: ratio of medians %.3f (1 - r %.3f) | median of per-replicate ratios %.3f (1 - r %.3f)\n', ...
            rMed, 1 - rMed, rPer, 1 - rPer);
        pr('CONTROL2 - CONTROL d %+.4f [%+.4f, %+.4f] | CONTROL mean fitness (load window, median) %.4f | mutation cost/gen CONTROL %.4f | AT BALANCE [C C2 A X] %s\n', ...
            T(1, 2).d, ciN, fitC, median(Y(:, 1, 3), 'omitnan'), mat2str(bal));
        S.(A(k).name) = struct('rMed', rMed, 'rPer', rPer);
    end
    fclose(fid);
end


function both(fid, varargin)
    fprintf(varargin{:}); fprintf(fid, varargin{:});
end
