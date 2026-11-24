function S = repair_test_summary(mode)
% REPAIR_TEST_SUMMARY  Does restoring retaliation cut collapse by more than the sham? Compares
% collapse within 100 generations in the observed arm, the repaired arm and the sham, with and
% without the defector that later took over, gives the part of the reduction specific to restoring
% retaliation, and prints the checks of the design.
%
% Produces : SI S4.3 and Table S10; Fig. S1 (figS1_draw calls this script for its per-population
%            shares); main text, "How cooperating populations lose cooperation"
% Inputs   : mode, 'main' (the only one)
%            publication/repair_matched/repair_matched_build_targets.mat and repair_matched_main.mat
%              (run_repair_test, modes 'build' and 'main'; W99s: targets whose cooperation is at
%              least 0.99)
%            publication/repair_test_observed/repair_test_observed.mat   (run_repair_test_observed)
% Outputs  : publication/repair_matched/repair_matched.txt   (the printed report: [0] checks of
%              the design, [1] sham minus repair, [2] the other differences and the specific
%              share, [G1] starting cooperation, [AR] by architecture); S, the per-population
%              shares by arm (returned)
% Run      : repair_test_summary('main')   (about 1 min)
% Seed     : rng(20261013, 'twister'), then randi(n, n, 10000) over the n populations, the
%            resamples of every interval
%
% Arm labels, as in the report: A, A0 observed (with and without the defector); B', B0' repaired
% (`Bp`, `B0p`); S', S0' the sham (`Sp`, `S0p`). Every arm runs on the same seeds, ten runs per
% population.

    if nargin < 1, mode = 'main'; end
    assert(strcmp(mode, 'main'), 'repair_test_summary:mode', 'mode must be ''main''');
    pub = grn_data_file('publication');
    dir0 = fullfile(pub, 'repair_matched');
    Bd = load(fullfile(dir0, 'repair_matched_build_targets.mat'), 'out'); Bd = Bd.out;   % the build: targets with c >= 0.99
    R = load(fullfile(dir0, 'repair_matched_main.mat'), 'out'); R = R.out;
    I = load(fullfile(pub, 'repair_test_observed', 'repair_test_observed.mat'), 'out');
    T = Bd.targets;
    T = T([T.inRun] & ~[T.excluded]);
    assert(isequal(sort({T.id}), sort(R.events)), 'run events differ from the build''s included events');
    txtFile = fullfile(dir0, 'repair_matched.txt');
    if exist(txtFile, 'file'), delete(txtFile); end
    diary(txtFile);
    cl = onCleanup(@() diary('off'));
    header(Bd, R, T);
    X = collect(T, R.runs, I.out.runs);
    rng(20261013, 'twister');
    bi = randi(numel(T), numel(T), 10000);
    S = struct('X', X);
    S.primary = primary(X, bi);
    S.secondary = secondary(X, bi);
    gen1(X, T);
    by_arch(X, T);
    fprintf('\n');
end


% =========================================================================
function header(Bd, R, T)
    fprintf('REPAIR TEST, TARGETS FILTERED UNDER THE RUNS'' ENCOUNTER RULE (MAIN)  runs: experiments/run_repair_test.m,\n');
    fprintf('  report: analysis/s4_losing_cooperation/repair_test_summary.m\n');
    fprintf('  design: tau = %.2f, pool %s, source-matched sham %d\n', ...
        Bd.cfg.tau, Bd.cfg.pool, Bd.cfg.matchSource);
    fprintf('  events %d; arms %s x k = 10 seeds (1..10); %d jobs in %.1f min; %.2f worker-s per run\n', ...
        numel(T), strjoin(R.arms, ' '), R.nJobs, R.runMin, R.secPerRun);
    fprintf('  stored arms A, A0 (repair_test_observed.mat), same seeds\n\n');
    x = Bd.exclusion; d = Bd.cdist; s = Bd.cstart;
    fprintf('[0] CONTROLS (from repair_matched_build_targets.mat)\n');
    fprintf('  exclusion: %d of 68 main events excluded (%s); run-time events excluded %d; stop rule fired: %s\n', ...
        x.nMain, strjoin(x.idsMain, ', '), x.nPilot, yn(x.stop));
    fprintf('  C-code (unfiltered target code reproduces the stored repair and sham populations): %s\n', ...
        pf(all([Bd.targets.codeOK])));
    fprintf('  C-dist: sham/repair distance ratio median %.4f [IQR %.3f, %.3f], %d individuals; same individuals: YES: %s\n', ...
        d.median, d.iqr, d.n, pf(d.pass));
    fprintf(['  C-start (%d main events): median(repair - observed) %+.4f, median(sham - observed) %+.4f, ', ...
        'repair >= observed - 0.05 in %d (%.1f%%): %s\n\n'], s.n, s.medBA, s.medSA, s.nB, 100 * s.shareB, pf(s.pass));
end


% =========================================================================
function X = collect(T, runs, inj)
% Per population: collapse share by arm (NaN where the arm was not run) and recorded generation-1
% cooperation.
    n = numel(T);
    arms = {'A', 'A0', 'Bp', 'Sp', 'B0p', 'S0p'};
    for a = 1:numel(arms), X.(['p' arms{a}]) = nan(n, 1); X.(['h' arms{a}]) = nan(n, 1); end
    X.arch = [T.arch]'; X.id = {T.id}';
    X.coop1A = [T.coop1A]'; X.coop1Bp = [T.coop1Bp]'; X.coop1Sp = [T.coop1Sp]';
    for i = 1:n
        m = T(i).num;
        for a = 1:numel(arms)
            if any(strcmp(arms{a}, {'A', 'A0'})), pool = inj; else, pool = runs; end
            r = pool([pool.event] == m & strcmp({pool.arm}, arms{a}));
            if isempty(r), continue; end
            kExp = 10;
            assert(numel(r) == 1 && numel(r.r.collapse) == kExp, 'event %s arm %s: not one job of %d runs', T(i).id, arms{a}, kExp);
            X.(['p' arms{a}])(i) = mean(r.r.collapse);
            X.(['h' arms{a}])(i) = mean(r.r.hc(:, 1));
        end
    end
end


% =========================================================================
% D = mean(p_S' - p_B'), the part specific to restoring retaliation.
function V = primary(X, bi)
    [D, ci] = boot(X.pSp - X.pBp, bi);
    fprintf('[1] D = mean(p_S'' - p_B''), the retaliation-specific effect (n = %d events)\n', size(bi, 1));
    fprintf('  mean p_B'' %.3f, mean p_S'' %.3f; D = %+.3f [%+.3f, %+.3f]\n\n', mean(X.pBp), mean(X.pSp), D, ci);
    V = struct('D', D, 'ci', ci);
end


function V = secondary(X, bi)
    fprintf('[2] also reported\n');
    V.AS = stat_line('A - S''   (the sham''s own effect)', X.pA - X.pSp, bi);
    V.AB = stat_line('A - B''   (the total repair effect)', X.pA - X.pBp, bi);
    V.ratio = ratio_line('mean p_A / mean p_B''', X.pA, X.pBp, bi);
    V.A0S0 = stat_line('A0 - S0'' (no defector added)', X.pA0 - X.pS0p, bi);
    V.A0B0 = stat_line('A0 - B0'' (no defector added)', X.pA0 - X.pB0p, bi);
    V.ratio0 = ratio_line('mean p_A0 / mean p_B0''', X.pA0, X.pB0p, bi);
    V.share = ratio_line('D / (A - B'') (retaliation-specific share)', X.pSp - X.pBp, X.pA - X.pBp, bi);
    fprintf('  means: p_A %.3f, p_S'' %.3f, p_B'' %.3f | p_A0 %.3f, p_S0'' %.3f, p_B0'' %.3f\n\n', ...
        mean(X.pA), mean(X.pSp), mean(X.pBp), mean(X.pA0), mean(X.pS0p), mean(X.pB0p));
end


function v = stat_line(name, x, bi)
    [m, ci] = boot(x, bi);
    fprintf('  %-46s %+.3f [%+.3f, %+.3f]\n', name, m, ci);
    v = [m ci];
end


function v = ratio_line(name, num, den, bi)
    m = mean(num) / mean(den);
    b = mean(num(bi), 1) ./ mean(den(bi), 1);
    ci = quantile(b, [0.025 0.975]);
    fprintf('  %-46s %.3f [%.3f, %.3f]\n', name, m, ci);
    v = [m ci];
end


% =========================================================================
function gen1(X, T)
    fprintf('[G1] GENERATION-1 COOPERATION (reported)\n');
    fprintf('  the runs'' encounter rule, no defector added, mean of 5 evaluations (median over events): A %.3f | B'' %.3f | S'' %.3f\n', ...
        median(X.coop1A), median(X.coop1Bp), median(X.coop1Sp));
    fprintf('  recorded generation-1 cooperation of the runs, seeds 1..10 (mean over events): A %.3f | B'' %.3f | S'' %.3f\n', ...
        mean(X.hA), mean(X.hBp), mean(X.hSp));
    fprintf('  individuals replaced per event: median %d [%d, %d]\n\n', median([T.nReplaced]), ...
        min([T.nReplaced]), max([T.nReplaced]));
end


function by_arch(X, T)
    fprintf('[AR] BY ARCHITECTURE (descriptive): n | p_A | p_S'' | p_B'' | S'' - B''\n');
    for a = 1:4
        k = X.arch == a;
        if ~any(k), continue; end
        fprintf('    arch %d (%s) n %2d | %.3f | %.3f | %.3f | %+.3f\n', a, T(find(k, 1)).id(1:2), sum(k), ...
            mean(X.pA(k)), mean(X.pSp(k)), mean(X.pBp(k)), mean(X.pSp(k) - X.pBp(k)));
    end
end


function [m, ci] = boot(x, bi)
    m = mean(x); b = mean(x(bi), 1); ci = quantile(b, [0.025 0.975]);
end


function s = yn(t), if t, s = 'YES'; else, s = 'NO'; end, end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
