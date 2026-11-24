function out = selection_inherited(outFile)
% SELECTION_INHERITED  When does the part of selection on cooperation that offspring inherit turn
% in cooperation's favor? Splits selection on cooperation, generation by generation, in 20
% populations of each signaling architecture.
%
% Produces : the first generation favoring cooperation, the shares of generations 16 to 500 against
%            cooperation, not detectable and for it, and the generation of lowest mean cooperation
%            (Table S4; SI S2.4); the summary the strip of Fig. 2A is checked against
% Inputs   : outFile   (default publication/selection_over_time/summary_selection_inherited.txt)
%            publication/selection_over_time/, for the four architectures and replicates 1 to 20:
%              selection_z_<arch>_r<NN>.mat (run_selection_z; cooperation Z and parents),
%              selection_sign_<arch>_r<NN>.mat (run_selection_sign; stored S, mean cooperation),
%              lineage_records_<arch>_r<NN>.mat (run_lineage_records; each individual's wiring)
% Outputs  : outFile and the same text on the console: the checks per run, the saturating-contrast
%            states at generations 10 to 20, the table, the correlations, the means with
%            intervals by generation
% Run      : selection_inherited   (not timed)
% Seed     : rng(20261022) (MATLAB's default generator) at the start of each bootstrap, run twice
%            per measure to check that it reproduces
%
% With n_i the offspring of individual i and z_i its cooperation: realized S = Cov(n, z); inherited
% S_G = Cov(n, g), with g_i the mean z of the individuals of that generation sharing i's wiring; the
% rest S - S_G = Cov(n, z - g) (all with 1/N). A generation favors cooperation when the 95%
% bootstrap interval of the 20-run mean (10,000 resamples of runs) lies above zero.
    cfg = struct( ...
        'pub', grn_data_file('publication'), ...
        'archs', {{'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'}}, ...
        'nGen', 500, 'bootN', 10000, 'bootSeed', 20261022);
    cfg.matDir = fullfile(cfg.pub, 'selection_over_time');
    if nargin < 1 || isempty(outFile)
        outFile = fullfile(grn_data_file('publication'), 'selection_over_time', 'summary_selection_inherited.txt');
    end
    unit_test();
    out = summary_of(cfg, outFile);
end

function unit_test()
% Individuals 1 and 2 share wiring 1 and differ only in partner luck (z 1 and 0); 1 leaves two
% offspring and 2 none, so S > 0 and S_G = 0.
    [S, SG, SL, nd] = parts_of([1 1 3 4], [1 0 0.5 0.5], [1 1 2 3]);   % n = [2 0 1 1]
    assert(abs(S - 0.25) < 1e-15 && abs(SG) < 1e-15 && abs(SL - 0.25) < 1e-15 && nd == 3, 'unit test 1 failed');
    % distinct wirings, no luck: g = z, S_G = S, S_L = 0
    [S, SG, SL] = parts_of([1 1 1 4], [0 1 1 0.5], [1 2 3 4]);          % hand value -0.5
    assert(abs(S + 0.5) < 1e-15 && abs(SG + 0.5) < 1e-15 && abs(SL) < 1e-15, 'unit test 2 failed');
    fprintf('unit test: PASS (identical wirings differing only in partner luck: S = 0.25, S_G = 0, S_L = 0.25; distinct wirings: S_G = S)\n');
end

function out = summary_of(cfg, outFile)
    fid = fopen(outFile, 'w'); cl = onCleanup(@() fclose(fid));
    w = @(varargin) wboth(fid, varargin{:});
    nA = numel(cfg.archs); G = cfg.nGen;
    gl = [1 3 5 7 10 15 20 30 50 100 200 500];
    w('selection_inherited summary (%s)\n', datestr(now, 'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>
    w('g_i = mean z over individuals sharing i''s wiring (idx pointer of run_lineage_records file); S_G = Cov(n, g); S_L = S - S_G\n\n');
    w('SELF-CHECKS\n  1. unit test (identical wirings differing only in partner luck: S > 0, S_G = 0): PASS (asserted at start)\n');
    w('  2. per run: z-file vs run_lineage_records parent record identical; stored S reproduced (max abs diff); S_L vs Cov(n, z - g) (max abs diff against S - S_G)\n');
    A = cell(nA, 1); bad = {}; maxS = 0; maxL = 0;
    for a = 1:nA
        A{a} = struct('S', [], 'SG', [], 'SL', [], 'meanZ', [], 'nd', [], 'hc', []);
        for r = 1:20
            tag = sprintf('%s_r%02d', cfg.archs{a}, r);
            fz = fullfile(cfg.matDir, ['selection_z_' tag '.mat']);
            fs = fullfile(cfg.matDir, ['selection_sign_' tag '.mat']);
            fr = fullfile(cfg.matDir, ['lineage_records_' tag '.mat']);
            if ~all(cellfun(@(f) exist(f, 'file'), {fz, fs, fr})), bad{end+1} = [tag ' missing']; continue; end %#ok<AGROW>
            Zf = load(fz, 'Z', 'par'); Sf = load(fs, 'S', 'meanZ', 'hc'); Rf = load(fr, 'par', 'idx');
            joinOK = isequal(Zf.par, Rf.par);
            if ~joinOK, w('     %s  PARENT RECORDS DIFFER: excluded\n', tag); bad{end+1} = [tag ' join']; continue; end %#ok<AGROW>
            S = zeros(G, 1); SG = S; SL = S; nd = S; dL = 0;
            for g = 1:G
                [S(g), SG(g), SL(g), nd(g)] = parts_of(Zf.par(:, g), Zf.Z(:, g), Rf.idx(:, g));
                dL = max(dL, abs(SL(g) - (S(g) - SG(g))));
            end
            dS = max(abs(S - Sf.S));
            ok = joinOK && dS <= 1e-12 && dL <= 1e-12;
            w('     %-20s r%02d  parent join %s | max|S - stored S| %.2g | max|S_L - (S - S_G)| %.2g | %s\n', ...
                cfg.archs{a}, r, pf(joinOK), dS, dL, ternary(ok, 'used', 'EXCLUDED'));
            maxS = max(maxS, dS); maxL = max(maxL, dL);
            if ~ok, bad{end+1} = tag; continue; end %#ok<AGROW>
            A{a}.S(:, end+1) = S; A{a}.SG(:, end+1) = SG; A{a}.SL(:, end+1) = SL;
            A{a}.meanZ(:, end+1) = Sf.meanZ; A{a}.nd(:, end+1) = nd; A{a}.hc(:, end+1) = Sf.hc;
        end
    end
    w('     runs used:'); for a = 1:nA, w('  %s %d', cfg.archs{a}, size(A{a}.S, 2)); end
    w('\n     excluded or missing: %d%s; max over runs |S - stored S| %.2g, |S_L - Cov(n, z - g)| %.2g (tolerance 1e-12)\n', ...
        numel(bad), ternary(isempty(bad), '', [' -> ' strjoin(bad, '; ')]), maxS, maxL);
    if any(cellfun(@(x) size(x.S, 2) == 0, A)), w('  an architecture has no usable run: stop\n'); out = A; return; end
    flds = {'S', 'SG', 'SL'}; ci = struct(); st = struct(); repro = true;
    for q = 1:3
        f = flds{q}; ci.(f) = cell(nA, 1); st.(f) = cell(nA, 1);
        c1 = boot_all(cellfun(@(x) x.(f), A, 'UniformOutput', false), cfg);
        c2 = boot_all(cellfun(@(x) x.(f), A, 'UniformOutput', false), cfg);
        repro = repro && isequal(c1, c2);
        for a = 1:nA
            ci.(f){a} = c1{a}; s = zeros(G, 1); s(c1{a}(:, 1) > 0) = 1; s(c1{a}(:, 2) < 0) = -1; st.(f){a} = s;
        end
    end
    w('  3. bootstrap reproducible (rng(%d), %d resamples, run twice, S, S_G, S_L): %s\n\n', cfg.bootSeed, cfg.bootN, pf(repro));
    sn = @(s) ternary(s > 0, 'coop pays', ternary(s < 0, 'defect pays', 'no clear adv'));

    % Saturating contrast: the S_G state at each generation 10 to 20
    k = sum(st.SG{4}(10:20) == 1);
    w('    saturating_contrast: states at g10..20 under S_G:'); w(' %d', st.SG{4}(10:20)); w('   (1 coop pays, 0 no clear, -1 defect pays)\n');
    w('    cooperating pays in %d of 11\n\n', k);

    % Table S4: the first generation favoring cooperation under S_G and S, the generation of the
    % lowest mean cooperation, the shares of generations 16 to 500 in each state; then the
    % correlations with the change in mean z and the means by generation
    w('REPORT\n  First "cooperating pays" generation (CI above 0), S_G against S, and the 20-run mean curve minimum:\n');
    for a = 1:nA
        fS = find(st.S{a} == 1, 1); fG = find(st.SG{a} == 1, 1);
        [~, gm] = min(mean(A{a}.hc, 2));
        w('    %-20s S_G %s | S %s | curve minimum g %d\n', cfg.archs{a}, num_(fG), num_(fS), gm);
    end
    w('  Share of generations 16-500 in each state:\n');
    for a = 1:nA
        for q = [2 1]
            s2 = st.(flds{q}){a}(16:G);
            w('    %-20s %-3s defecting pays %.3f, no clear advantage %.3f, cooperating pays %.3f\n', cfg.archs{a}, strrep(flds{q}, 'SG', 'S_G'), mean(s2 == -1), mean(s2 == 0), mean(s2 == 1));
        end
    end
    w('  Correlation over generations 1-499 with the change in mean z:\n');
    for a = 1:nA
        dz = diff(A{a}.meanZ, 1, 1); mdz = mean(dz, 2);
        cG = corr(mean(A{a}.SG(1:G-1, :), 2), mdz); cS = corr(mean(A{a}.S(1:G-1, :), 2), mdz);
        rG = arrayfun(@(k) corr(A{a}.SG(1:G-1, k), dz(:, k)), 1:size(dz, 2));
        w('    %-20s ensemble-mean series: S_G %.3f, S %.3f | per-run S_G: median %.3f [min %.3f, max %.3f]\n', cfg.archs{a}, cG, cS, median(rG), min(rG), max(rG));
    end
    w('  Ensemble mean [95%% bootstrap CI] by generation; S_G state; distinct wirings (median over runs)\n');
    for a = 1:nA
        w('  %s (%d runs)\n    gen   S                          S_G                        S_L                        S_G state     distinct wirings\n', cfg.archs{a}, size(A{a}.S, 2));
        for g = gl
            w('    %3d', g);
            for q = 1:3
                f = flds{q}; w('  %+.4f [%+.4f,%+.4f]', mean(A{a}.(f)(g, :)), ci.(f){a}(g, :));
            end
            w('  %-12s  %d\n', sn(st.SG{a}(g)), round(median(A{a}.nd(g, :))));
        end
    end
    out = A;
end

function s = num_(v), if isempty(v), s = 'never'; else, s = sprintf('%d', v); end, end

function ci = boot_all(Sc, cfg)
% Bootstrap over runs, as in run_selection_sign: 2.5 and 97.5 percentiles of the mean over a
% resample of runs, at every generation; the seed is reset at each call.
    rng(cfg.bootSeed);
    ci = cell(numel(Sc), 1); chunk = 1000;
    for a = 1:numel(Sc)
        S = Sc{a}; nr = size(S, 2); G = size(S, 1);
        bm = zeros(G, cfg.bootN);
        for b0 = 1:chunk:cfg.bootN
            nb = min(chunk, cfg.bootN - b0 + 1);
            idx = randi(nr, nr, nb);
            X = reshape(S(:, idx(:)), G, nr, nb);
            bm(:, b0:b0 + nb - 1) = reshape(mean(X, 2), G, nb);
        end
        ci{a} = prctile(bm, [2.5 97.5], 2);
    end
end

function wboth(fid, varargin)
    fprintf(varargin{:}); fprintf(fid, varargin{:});
end
function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
