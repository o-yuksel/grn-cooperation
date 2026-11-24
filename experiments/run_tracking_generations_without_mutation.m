function out = run_tracking_generations_without_mutation(cfg)
% RUN_TRACKING_GENERATIONS_WITHOUT_MUTATION  With the target held still, do the cuts in the
% mutation load depend on how long selection runs without mutation before the load left is read?
% Runs the population tests with 150 such generations and reads the cuts after 30 to 150.
%
% Produces : Fig. 6; SI S5.1; Table S14 ('run' writes the runs tracking_cut_after_150 reads;
%            'summary' prints the cuts after 30 generations without mutation and the fall of the
%            load left from 120 to 150 that SI S5.1 states)
% Inputs   : cfg (optional): `mode`, the phase ('run', or 'summary' by default); outDir
%              (publication/tracking_held_still_150gen), nOff (150), Ts ([30 60 90 120 150]), designs
%              ([1 2]), repList (1:20) and nWorkers (2) of the runs
%            'run': those of run_tracking_held_still (tracking_runs_N1000.mat, lowering_settings.mat)
%            'summary': <outDir>/design1 and design2/main.mat ('run'); and the runs with 30
%              generations without mutation, publication/tracking_held_still_30gen/design1/main.mat
%              (run_tracking_held_still) and tracking_held_still_30gen/design2/main.mat
%              (run_tracking_held_still_design2)
% Outputs  : 'run': <outDir>/design1 and design2, main.mat and main.txt
%            'summary': <outDir>/summary.txt (what it prints) and reanalysis_log.txt (the
%              runner's analysis at each T)
% Run      : run_tracking_generations_without_mutation(struct('mode', 'run'))   (design 1 about
%              37 min, design 2 about 68 min, on 2 workers), then 'summary' (about 25 s)
% Seed     : 'run': run_tracking_held_still's (replicate r continues with rng(20261001 + r,
%            'twister')); 'summary' draws nothing here; the runner's bootstrap uses 20261002
%
% 'run' calls run_tracking_held_still with nOff = 150 for design 1 (designIndex 1, a = 0.70) and design
% 2 (designIndex 3, a = 0.20), all else the runner's default. The generations without mutation continue
% the same random stream, so their first 30 equal those of the runs with 30; 'summary' checks
% this (outcome matrix and NaN pattern). For each T it cuts the saved run to its first T such
% generations and passes it to the runner's analysis (`cfg.reanalyse`), which reads the load left
% over the last 10 of them.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    pub = fullfile(grn_data_file('publication'));
    def = struct('mode', 'summary', 'outDir', fullfile(pub, 'tracking_held_still_150gen'), 'nOff', 150, ...
        'Ts', [30 60 90 120 150], 'designs', [1 2], 'repList', 1:20, 'nWorkers', 2);
    f = fieldnames(def);
    for i = 1:numel(f), if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end, end
    D = struct('designIndex', {1, 3}, 'tol', {0.0106, 0.0110}, 'dir', {'design1', 'design2'}, ...
        'pubDir', {'tracking_held_still_30gen/design1', 'tracking_held_still_30gen/design2'});
    switch cfg.mode
        case 'run'
            for k = cfg.designs
                c = struct('designIndex', D(k).designIndex, 'tolerance', D(k).tol, 'nOff', cfg.nOff, ...
                    'outDir', fullfile(cfg.outDir, D(k).dir), 'repList', cfg.repList, 'nWorkers', cfg.nWorkers);
                out = run_tracking_held_still(c);
            end
        case 'summary'
            out = summary(cfg, D, pub);
    end
end


function out = summary(cfg, D, pub)
    out = struct('ctrl', false(1, 2), 'settled', false(1, 2));
    logF = fullfile(cfg.outDir, 'reanalysis_log.txt'); if exist(logF, 'file'), delete(logF); end
    sumF = fullfile(cfg.outDir, 'summary.txt'); if exist(sumF, 'file'), delete(sumF); end
    arms = {'CONTROL', 'CONTROL2', 'ALPHA-DOWN', 'XNORM-DOWN'};
    for k = cfg.designs
        S = load(fullfile(cfg.outDir, D(k).dir, 'main.mat')); run = S.out;
        keep = run.verif.rep;
        nT = numel(cfg.Ts); Yt = cell(nT, 1); Tab = cell(nT, 1);
        for t = 1:nT
            T = cfg.Ts(t); o = run; o.cfg.nOff = T;
            for j = 1:numel(o.jobs)
                z = o.jobs{j};
                if isempty(z.off), continue; end
                for fn = {'meanFit', 'wBar', 'wBest', 'varFit'}
                    z.off.(fn{1}) = z.off.(fn{1})(1:T);
                end
                o.jobs{j} = z;
            end
            S2 = struct('out', o); tmp = [tempname '.mat']; save(tmp, '-struct', 'S2');
            c = struct('reanalyse', tmp, 'tolerance', D(k).tol, 'designIndex', D(k).designIndex, 'outDir', cfg.outDir);
            txt = evalc('R = run_tracking_held_still(c);');
            delete(tmp);
            append_text(logF, sprintf('\n\n######## design %d, T = %d ########\n%s', k, T, txt));
            Yt{t} = R.analysis.Y; Tab{t} = R.analysis.table;
        end
        fid = fopen(sumF, 'a');
        pr = @(varargin) fprintf_both(fid, varargin{:});
        pr('\n================ design %d (cfg.designIndex %d) ================\n', k, D(k).designIndex);
        % Check: at T = 30 the outcome matrix and its NaN pattern must equal those of the run
        % with 30 generations without mutation
        P = load(fullfile(pub, D(k).pubDir, 'main.mat')); Yp = P.out.analysis.Y;
        dmax = max(abs(Yt{1}(:) - Yp(:)), [], 'omitnan');
        sameNaN = isequal(isnan(Yt{1}), isnan(Yp));
        ctrlPass = dmax == 0 && sameNaN;
        out.ctrl(k) = ctrlPass;
        pr('CONTROL T = 30 vs the stored run %s/main.mat: max |diff| over the outcome matrix %.3g, NaN pattern equal %d -> %s\n', ...
            D(k).pubDir, dmax, sameNaN, pf(ctrlPass));
        % Cuts at every T: the mean paired difference from CONTROL over CONTROL's median
        % mutation load (outcome 1)
        pr('mutation load and cuts by T, the generations without mutation (cut = paired mean difference / CONTROL median, percentile bootstrap):\n');
        for t = 1:nT
            Tb = Tab{t}; c0 = median(Yt{t}(:, 1, 1), 'omitnan');
            pr('  T = %3d: CONTROL mutation load %.4f | ALPHA-DOWN d %+.4f [%+.4f, %+.4f] = %+.1f%% [%+.1f, %+.1f] (n %d) | XNORM-DOWN d %+.4f [%+.4f, %+.4f] = %+.1f%% [%+.1f, %+.1f] (n %d) | CONTROL2 d %+.4f [%+.4f, %+.4f]\n', ...
                cfg.Ts(t), c0, Tb(1,3).d, Tb(1,3).ci, 100*Tb(1,3).d/c0, 100*Tb(1,3).ci/c0, Tb(1,3).n, ...
                Tb(1,4).d, Tb(1,4).ci, 100*Tb(1,4).d/c0, 100*Tb(1,4).ci/c0, Tb(1,4).n, Tb(1,2).d, Tb(1,2).ci);
        end
        % The load left at every T (outcome 5, offLoad; the runner's outcomes are mutLoad, load,
        % mutCost, benefit, offLoad, maint)
        pr('load left (offLoad, last 10 generations without mutation), median over replicates [mean +- SE]:\n');
        for a = 1:4
            ok = arm_ok(Yt{1}, a, keep);
            row = '';
            for t = 1:nT
                v = Yt{t}(ok, a, 5);
                row = [row sprintf(' T%d %.4f [%.4f+-%.4f]', cfg.Ts(t), median(v), mean(v), std(v)/sqrt(numel(v)))]; %#ok<AGROW>
            end
            pr('  %-10s n %2d:%s\n', arms{a}, sum(ok), row);
        end
        % Settled rule, per arm: the mean change of the load left from T = 120 to 150 is smaller
        % than 5% of its mean at 120 or than the standard error of that change
        pr('SETTLED RULE (per arm): |mean(offLoad150 - offLoad120)| < max(0.05 * mean(offLoad120), SE of the paired difference)\n');
        allS = true;
        for a = 1:4
            ok = arm_ok(Yt{1}, a, keep);
            y120 = Yt{cfg.Ts == 120}(ok, a, 5); y150 = Yt{cfg.Ts == 150}(ok, a, 5);
            dd = y150 - y120; md = mean(dd); se = std(dd) / sqrt(numel(dd)); thr = max(0.05 * mean(y120), se);
            s = abs(md) < thr; allS = allS && s;
            pr('  %-10s mean diff %+.5f (%+.1f%% of offLoad120 %.5f), SE %.5f, threshold %.5f -> %s\n', ...
                arms{a}, md, 100*md/mean(y120), mean(y120), se, thr, ternary(s, 'settled', 'NOT settled'));
        end
        out.settled(k) = allS;
        pr('DESIGN %d RESULT: %s\n', k, ternary(allS, 'SPLIT SETTLED (cuts at T = 150 above are the settled-tail cuts)', 'NOT SETTLED (curve reported; no settled cut)'));
        % Population fitness variance along the generations without mutation
        pr('population fitness variance in the generations without mutation (median over replicates) at generation [1 30 60 90 120 150] of them:\n');
        for a = 1:4
            ok = find(arm_ok(Yt{1}, a, keep)); V = nan(numel(ok), 6);
            J = [run.jobs{:}];
            for q = 1:numel(ok)
                z = J([J.r] == run.cfg.repList(ok(q)) & [J.arm] == a);
                V(q, :) = z.off.varFit([1 30 60 90 120 150]);
            end
            pr('  %-10s %s\n', arms{a}, mat2str(median(V, 1, 'omitnan'), 3));
        end
        fclose(fid);
    end
end


function ok = arm_ok(Y, a, keep)
% Replicates with a finite load left in arm a and in CONTROL; in the two lowered arms (a = 3, 4)
% only those the runner kept (keep = out.verif.rep, the solve for the two settings (b, delta) that lower alpha
% or ||x*||, verified at the split)
    ok = isfinite(Y(:, a, 5)) & isfinite(Y(:, 1, 5));
    if a >= 3, ok = ok & keep(:, a-2); end
end


function fprintf_both(fid, varargin)
    fprintf(varargin{:}); fprintf(fid, varargin{:});
end


function append_text(f, s)
    fid = fopen(f, 'a'); fprintf(fid, '%s', s); fclose(fid);
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
