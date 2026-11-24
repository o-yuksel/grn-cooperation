function R = selection_around_collapse(dataRoot)
% SELECTION_AROUND_COLLAPSE  Does inherited selection turn against cooperation before a population
% loses it? Reads inherited selection on cooperation, S_G, in windows around each collapse, and
% over the late generations of populations that held cooperation.
%
% Produces : the lost and held counts per population (Table S4; SI S2.4)
% Inputs   : dataRoot (default: the data root's publication folder), under which, for the four
%            architectures and replicates 1 to 20:
%            selection_over_time/selection_z_<arch>_r<NN>.mat   (run_selection_z; cooperation Z and parents)
%            selection_over_time/lineage_records_<arch>_r<NN>.mat   (run_lineage_records; each individual's wiring)
%            selection_over_time/selection_sign_<arch>_r<NN>.mat   (run_selection_sign; S, for the check)
%            <arch>_100reps.mat   (run_signaling_architectures; cooperation curves, for the check)
% Outputs  : publication/selection_around_collapse/selection_around_collapse.txt   (per architecture and pooled, the
%            window counts and onsets; one row per collapsed population; nothing on the console)
% Run      : selection_around_collapse   (not timed)
% Seed     : rng(1) (MATLAB's default generator), once; the block bootstrap draws in a fixed order
%
% S = Cov(n, z) and S_G = Cov(n, g), with n_i the offspring of individual i, z_i its cooperation and
% g_i the mean z of the individuals sharing its wiring; S_G < 0 is selection against cooperation.
% g90 is the first generation with mean cooperation above 0.90; the collapse gc, the first from g90
% on below 0.5. Collapsed: g90 and gc exist; held: g90 only. Table S4 reads the window means (lost:
% mean S_G < 0 over gc-30..gc-1; held: mean S_G > 0 over g90..500), which do not depend on the
% random stream; the bootstrap intervals do.
    if nargin < 1 || isempty(dataRoot)
        dataRoot = fullfile(grn_data_file('publication'));
    end
    if isempty(dir(fullfile(dataRoot, 'selection_over_time', 'selection_z_*.mat')))
        error('selection_around_collapse:inputs', 'Missing inputs: no selection_z_*.mat under %s\\selection_over_time', dataRoot);
    end
    outDir = fullfile(grn_data_file('publication'), 'selection_around_collapse');
    if ~isfolder(outDir), mkdir(outDir); end
    rng(1);
    R = collect(dataRoot);
    fid = fopen(fullfile(outDir, 'selection_around_collapse.txt'), 'w'); cl = onCleanup(@() fclose(fid));
    report(fid, R, dataRoot);
end

function R = collect(dataRoot)
% Per population: S and S_G at every generation, the controls, and the windows. Control 1 (stops
% the run above 1e-10): S recomputed from Z equals the stored S. Control 2 (reported): the mean
% of Z against the stored cooperation curve of the 100-replicate file.
    modes = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'};
    R = []; dS = 0; dC = 0; ed = fullfile(dataRoot, 'selection_over_time');
    for m = 1:4
        cf = load_coop(dataRoot, modes{m});
        for r = 1:20
            tag = sprintf('%s_r%02d', modes{m}, r);
            Zf = load(fullfile(ed, ['selection_z_' tag '.mat']), 'Z', 'par');
            Rf = load(fullfile(ed, ['lineage_records_' tag '.mat']), 'idx');
            Sf = load(fullfile(ed, ['selection_sign_' tag '.mat']), 'S');
            [S, SG] = series(Zf, Rf);
            mc = mean(double(Zf.Z), 1)';
            dS = max(dS, max(abs(S - Sf.S)));
            dC = max(dC, max(abs(mc - cf(r, 1:500)')));
            e = test_membership(modes{m}, r, mc, SG);
            if isempty(R), R = e; else, R(end+1) = e; end %#ok<AGROW>
        end
    end
    if dS > 1e-10, error('selection_around_collapse:control', 'recomputed S differs from stored S: %g', dS); end
    for k = 1:numel(R), R(k).dS = dS; R(k).dC = dC; end
end

function cf = load_coop(dataRoot, mode)
    f = fullfile(dataRoot, [mode '_100reps.mat']);
    if ~isfile(f), error('selection_around_collapse:inputs', 'Missing input %s', f); end
    d = load(f, 'results');
    cf = zeros(numel(d.results), 500);
    for q = 1:numel(d.results), cf(q, :) = d.results{q}.historyCooperation(1:500); end
end

function [S, SG] = series(Zf, Rf)
    S = zeros(500, 1); SG = S;
    for g = 1:500, [S(g), SG(g)] = parts_of(Zf.par(:, g), Zf.Z(:, g), Rf.idx(:, g)); end
end

function e = test_membership(arch, r, mc, SG)
% Windows of a collapsed population: g90..gc-50, gc-30..gc-1, gc..gc+30, gc+31..500; of a held
% one: g90..500 and 300..500 (each cut to generations 1 to 500). A population that never exceeds
% 0.90 is 'other' and gets no window.
    e = struct('arch', arch, 'r', r, 'cls', 'other', 'g90', NaN, 'gc', NaN, ...
        'win', [], 'pn', NaN, 'dS', 0, 'dC', 0);
    g90 = find(mc > 0.90, 1);
    if isempty(g90), return; end
    e.g90 = g90; gc = g90 - 1 + find(mc(g90:end) < 0.5, 1);
    if isempty(gc), e.cls = 'held'; ws = {g90:500; 300:500};
    else
        e.cls = 'collapsed'; e.gc = gc;
        ws = {g90:gc-50; max(gc-30, 1):gc-1; gc:min(gc+30, 500); gc+31:500};
        e.pn = onset(SG, g90, gc) - gc;
    end
    for k = 1:numel(ws)
        ii = ws{k}; ii = ii(ii >= 1 & ii <= 500);
        s = summ(SG(ii), mean(mc(ii)));
        if isempty(e.win), e.win = s; else, e.win(end+1) = s; end %#ok<AGROW>
    end
end

function pn = onset(SG, g90, gc)
% Persistent-negative onset: the first g >= g90 + 9 from which the 10-generation trailing mean of
% S_G stays below 0 through gc (NaN if none). It finds negativity that lasts to the collapse, so
% it leans toward late onsets.
    rm = movmean(SG, [9 0]); pn = NaN;
    for g = g90+9:gc
        if all(rm(g:gc) < 0), pn = g; return; end
    end
end

function s = summ(x, mcoop)
% Window mean of S_G with a 95% block-bootstrap interval over generations (blocks of 10,
% 2,000 resamples), and the shares of its generations below and above zero.
    x = x(:); n = numel(x);
    s = struct('n', n, 'mean', NaN, 'lo', NaN, 'hi', NaN, 'neg', NaN, 'pos', NaN, 'mcoop', mcoop);
    if n == 0, return; end
    B = 2000; L = min(10, n); nb = ceil(n/L); bm = zeros(B, 1);
    for b = 1:B
        st = randi(max(n-L+1, 1), nb, 1);
        ix = bsxfun(@plus, st, 0:L-1)'; ix = ix(:); ix = ix(1:n);
        bm(b) = mean(x(ix));
    end
    ci = prctile(bm, [2.5 97.5]);
    s.mean = mean(x); s.lo = ci(1); s.hi = ci(2); s.neg = mean(x < 0); s.pos = mean(x > 0);
end

function report(fid, R, dataRoot)
    P = @(varargin) fprintf(fid, varargin{:});
    P('selection_around_collapse (generated by analysis/s2_share_of_members/selection_around_collapse.m)\n');
    P('inputs: %s\n', dataRoot);
    P('Block bootstrap (block 10 gens, 2000 resamples, rng(1)), 95%% CI on window mean of S_G.\n\n');
    P('CONTROL 1: recomputed S equals stored S within 1e-10 for all %d populations (max abs diff %.3g): PASS\n', numel(R), R(1).dS);
    P('CONTROL 2: mean of Z vs stored per-replicate cooperation curve (100reps file, replicates 1-20, gens 1-500): max abs diff %.3g\n\n', R(1).dC);
    arch = {R.arch}; cls = {R.cls};
    for a = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast', 'ALL'}
        sel = strcmp(cls, 'collapsed') & (strcmp(a{1}, 'ALL') | strcmp(arch, a{1}));
        P('%s collapsed n=%d\n', a{1}, sum(sel)); block(P, R(sel));
    end
    for a = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast', 'ALL'}
        sel = strcmp(cls, 'held') & (strcmp(a{1}, 'ALL') | strcmp(arch, a{1}));
        P('%s held n=%d\n', a{1}, sum(sel)); hblock(P, R(sel));
    end
    P('\nTABLE collapsed populations: arch r g90 gc | A=g90..gc-50, B=gc-30..gc-1, C=gc..gc+30, D=gc+31..500: mean S_G [CI] share<0 | onset rel. collapse\n');
    letters = 'ABCD';
    for e = R(strcmp(cls, 'collapsed'))
        P('%s r%02d %d %d', e.arch, e.r, e.g90, e.gc);
        for k = 1:numel(e.win)
            w = e.win(k); P(' | %s:', letters(k));
            if w.n == 0, P('n/a'); else, P('%+.1e [%+.1e,%+.1e] <0:%.2f', w.mean, w.lo, w.hi, w.neg); end
        end
        P(' | onset %g\n', e.pn);
    end
end

function block(P, C)
    names = {'g90..gc-50', 'gc-30..gc-1', 'gc..gc+30', 'gc+31..500'};
    for k = 1:4
        X = arrayfun(@(e) numel(e.win) >= k && e.win(k).n > 0, C);
        if ~any(X), continue; end
        W = arrayfun(@(e) e.win(k), C(X));
        P('  %-12s runs=%d mean<0 in %d; CI wholly<0 in %d; CI wholly>0 in %d; median window mean=%+.2e; median share<0=%.2f\n', ...
            names{k}, numel(W), sum([W.mean] < 0), sum([W.hi] < 0), sum([W.lo] > 0), median([W.mean]), median([W.neg]));
    end
    pn = [C.pn]; pn = pn(~isnan(pn));
    P('  persistent-neg onset rel collapse: found in %d of %d; median %g; list %s\n', numel(pn), numel(C), median(pn), mat2str(sort(pn)));
end

function hblock(P, H)
    names = {'g90..500', '300..500'};
    if isempty(H), return; end
    for k = 1:2
        W = arrayfun(@(e) e.win(k), H);
        P('  %-9s runs=%d CI wholly>0 in %d; CI wholly<0 in %d; mean>0 in %d; median window mean=%+.2e; median share>0=%.2f\n', ...
            names{k}, numel(W), sum([W.lo] > 0), sum([W.hi] < 0), sum([W.mean] > 0), median([W.mean]), median([W.pos]));
    end
end
