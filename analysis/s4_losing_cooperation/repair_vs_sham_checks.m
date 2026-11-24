function out = repair_vs_sham_checks(cfg)
% REPAIR_VS_SHAM_CHECKS  Do the repair and the sham differ in retaliation and little else? Compares
% the wirings each moves individuals to (the targets): their cooperation and retaliation margins
% and mutation-sensitivity (check 1), the genetic diversity each arm leaves (check 2), and whether
% the targets retaliate under the runs' encounter rule (check 3). No run evolves; no collapse is read.
%
% Produces : SI S4.3, "Repair and sham targets differ in retaliation and little else"
% Inputs   : cfg (optional struct; defaults: nWorkers 3 (at most 3), nSeed 5, nBoot 10000,
%              bootSeed 20261013, m 10, rngSeed 20260928, nMemCtl 4, tag '')
%            publication/repair_matched/repair_matched_build_targets.mat   (run_repair_test, mode
%              'build': per population the replaced individuals and their repair and sham targets;
%              W99s: targets whose cooperation is at least 0.99)
%            publication/repair_restart_points/repair_restart_points_events.mat   (run_repair_restart_points)
% Outputs  : publication/repair_matched/repair_fairness.mat   (out: cfg, check1 to check3, minutes)
%            publication/repair_matched/repair_fairness.txt   (the report; the console prints
%              progress only)
% Run      : repair_vs_sham_checks()   (about 2.3 min on 3 workers). With no
%            pool open it opens one; the free-memory check and the lowered process priority
%            (memory, PowerShell) run on Windows only and are skipped elsewhere.
% Seed     : bootstrap resamples from rng(bootSeed, 'twister'); the slots of the ten added
%            defector copies from rng(rngSeed + e) for population e, the session's generator, as in
%            the runs; check 3's seeds beside check3; membership_test seeds each wiring itself
%
% Per replaced individual: itself at the restart point (a cooperator that does not retaliate; A
% in the report), its repair target (B') and its sham target (S'); control: the members at the
% restart point. A target taken by many individuals counts once per individual.

    if nargin < 1, cfg = struct(); end
    def = struct('nWorkers', 3, 'nSeed', 5, 'nBoot', 10000, 'bootSeed', 20261013, 'm', 10, ...
        'rngSeed', 20260928, 'nMemCtl', 4, 'tag', '');
    f = fieldnames(def);
    for i = 1:numel(f), if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end, end
    assert(cfg.nWorkers <= 3, 'use at most 3 workers (cfg.nWorkers)');
    pub = fullfile(grn_data_file('publication'));
    outDir = fullfile(pub, 'repair_matched');
    open_pool(cfg);
    t0 = tic;
    B = load(fullfile(outDir, 'repair_matched_build_targets.mat'), 'out'); T = B.out.targets;
    E = load(fullfile(pub, 'repair_restart_points', 'repair_restart_points_events.mat'), 'out'); ev = E.out.events; ev = ev([ev.usable]);
    T = T([T.inRun] & ~[T.excluded]);
    [~, ie] = ismember({T.id}, {ev.id}); ev = ev(ie);
    assert(numel(T) == 65 && isequal({T.id}, {ev.id}), 'expected the 65 included main events');
    for i = 1:numel(T)
        assert(isequal(T(i).ind(:), find(any(reshape(ev(i).Mpre ~= T(i).Mb, 16, []), 1) | ...
            any(reshape(ev(i).Bpre ~= T(i).Bb, 16, []), 1))'), '%s: replaced set', T(i).id);
    end

    % ---- gather populations and the distinct (arch, wiring) list
    P = struct('id', {}, 'arch', {}, 'ind', {}, 'GA', {}, 'GB', {}, 'GS', {}, 'GAi', {}, 'GBi', {}, 'GSi', {}, 'GMem', {}, 'nMem', {});
    for i = 1:numel(T)
        e = ev(i); t = T(i);
        P(i).id = t.id; P(i).arch = e.arch; P(i).ind = t.ind(:);
        P(i).GA = geno(e.Mpre, e.Bpre); P(i).GB = geno(t.Mb, t.Bb); P(i).GS = geno(t.Ms, t.Bs);
        [Mi, Bi] = with_later_defector(e, e.Mpre, e.Bpre, cfg); P(i).GAi = geno(Mi, Bi);
        [Mi, Bi] = with_later_defector(e, t.Mb, t.Bb, cfg);     P(i).GBi = geno(Mi, Bi);
        [Mi, Bi] = with_later_defector(e, t.Ms, t.Bs, cfg);     P(i).GSi = geno(Mi, Bi);
    end
    spA = arrayfun(@(e) e.sp, ev);
    spA = arrayfun(@(s) setfield(s, 'gradedC2', true), spA); %#ok<SFLD>

    % The members at each restart point (the control group), by the membership test
    fprintf('testing membership of the restart-point genotypes for the member control ...\n');
    [memRows, memCnt] = restart_members(P, spA);
    for i = 1:numel(P), P(i).GMem = memRows{i}; P(i).nMem = memCnt{i}; end

    % ---- CHECK 1: per target wiring, the cooperation margin (readout x_C - x_D against the clone at
    % x*), beta (membership_test: Frobenius norm of dx*/dW), the retaliation margin (minus the highest
    % readout against a defecting grid mutant; NaN when none defects) and beta over the cooperation
    % margin; per population the median over its replaced individuals, then B' - S' paired
    allArch = []; allW = {};
    kmap = containers.Map();
    for i = 1:numel(P)
        for G = {P(i).GA(P(i).ind, :), P(i).GB(P(i).ind, :), P(i).GS(P(i).ind, :), P(i).GMem}
            Gu = unique(G{1}, 'rows');
            for r = 1:size(Gu, 1)
                k = wkey(P(i).arch, Gu(r, :));
                if ~isKey(kmap, k), kmap(k) = kmap.Count + 1; allArch(end+1) = P(i).arch; allW{end+1} = reshape(Gu(r, 1:16) .* Gu(r, 17:32), 4, 4); end %#ok<AGROW>
            end
        end
    end
    nW = numel(allW);
    fprintf('CHECK 1: %d distinct (architecture, wiring) targets/controls to describe\n', nW);
    spArch = cell(1, 4);
    for i = 1:numel(T), spArch{ev(i).arch} = spA(i); end
    D = nan(nW, 6); parts = cell(nW, 1);
    parfor q = 1:nW
        parts{q} = describe(allW{q}, spArch{allArch(q)});
    end
    for q = 1:nW, D(q, :) = parts{q}; end   % margin, beta, retaliation margin, C1, member, worst_retal
    fprintf('  described in %.1f min\n', toc(t0) / 60);
    c1 = check1(P, kmap, D, cfg);

    % ---- CHECK 2: distinct genotypes (M, B) and wirings W entering generation 1, with and without
    % the added defector
    c2 = check2(P, cfg);

    % ---- CHECK 3: the class of each distinct target under the runs' encounter rule
    c3 = check3(P, spArch, cfg);

    out = struct('cfg', cfg, 'check1', c1, 'check2', c2, 'check3', c3, 'minutes', toc(t0) / 60);
    save(fullfile(outDir, ['repair_fairness' cfg.tag '.mat']), 'out', '-v7.3');
    write_report(out, fullfile(outDir, ['repair_fairness' cfg.tag '.txt']), numel(T));
    fprintf('done in %.1f min\n', out.minutes);
end


% =========================================================================
function G = geno(M, B)
    n = size(M, 3);
    G = [reshape(double(M), 16, n); reshape(double(B), 16, n)]';
end


function [M, B] = with_later_defector(e, M, B, cfg)
% The added defector placed as the runs place it (run_repair_test): cfg.m copies at random slots.
    rng(cfg.rngSeed + e.num);
    idx = randperm(size(M, 3), cfg.m);
    M(:,:,idx) = repmat(reshape(logical(e.invG(1:16)), 4, 4), 1, 1, cfg.m);
    B(:,:,idx) = repmat(reshape(single(e.invG(17:32)), 4, 4), 1, 1, cfg.m);
end


function k = wkey(a, g)
    w = g(1:16) .* g(17:32);
    k = sprintf('%d|%s', a, sprintf('%.9g,', w));
end


function open_pool(cfg)
% Opens the worker pool after a free-memory check: one thread per worker, below-normal priority
% (memory and the PowerShell call work on Windows only; on other systems the memory check and the
% priority change are skipped).
    if ~isempty(gcp('nocreate')), return; end
    if ispc
        [~, sys] = memory;
        freeGB = sys.PhysicalMemory.Available / 2^30; need = 1.2 * cfg.nWorkers;
        fprintf('available RAM before the pool: %.2f GB (need %.1f)\n', freeGB, need);
        assert(freeGB >= need, 'not enough free RAM for %d workers', cfg.nWorkers);
    end
    p = parpool('Processes', cfg.nWorkers);
    wait(parfevalOnAll(p, @maxNumCompThreads, 0, 1));
    if ispc
        f = parfevalOnAll(p, @() feature('getpid'), 1);
        pids = [feature('getpid'); fetchOutputs(f)];
        cmd = sprintf(['powershell -NoProfile -Command "Get-Process -Id %s | ForEach-Object ', ...
            '{ $_.PriorityClass = ''BelowNormal'' }"'], strjoin(arrayfun(@num2str, pids', 'uni', 0), ','));
        st = system(cmd);
        fprintf('pool: %d workers, 1 thread each; BelowNormal on PIDs %s (status %d)\n', cfg.nWorkers, mat2str(pids'), st);
    else
        fprintf('pool: %d workers, 1 thread each; priority and free memory not changed or checked on this system\n', cfg.nWorkers);
    end
end


function d = describe(W, sp)
% [cooperation margin, beta, retaliation margin, C1, member, worst_retal] of one wiring, by the
% membership test with every grid mutant scored (gradedC2).
    r = membership_test(W, sp);
    d = [r.x_star(1) - r.x_star(2), r.beta, -r.worst_retal, r.is_weak, r.is_css, r.worst_retal];
end


function [rows, cnt] = restart_members(P, spA)
% The members (C1 and C2 by the membership test) at each restart point: distinct genotypes, counts.
    n = numel(P); rows = cell(n, 1); cnt = cell(n, 1);
    parfor i = 1:n
        [Gu, ~, ig] = unique(P(i).GA, 'rows'); c = accumarray(ig, 1);
        keep = false(size(Gu, 1), 1);
        for u = 1:size(Gu, 1)
            r = membership_test(reshape(Gu(u, 1:16) .* Gu(u, 17:32), 4, 4), spA(i));
            keep(u) = r.is_weak && r.is_css;
        end
        rows{i} = Gu(keep, :); cnt{i} = c(keep);
    end
end


% =========================================================================
function c = check1(P, kmap, D, cfg)
    names = {'margin', 'beta', 'cret', 'beta_over_margin'};
    arms = {'A', 'Bp', 'Sp', 'Mem'};
    nE = numel(P);
    pooled = struct(); perEv = struct(); pdiff = nan(nE, 4); plog = nan(nE, 4); PD = nan(0, 4); PL = nan(0, 4);
    for a = arms, pooled.(a{1}) = nan(0, 4); perEv.(a{1}) = nan(nE, 4); end
    for i = 1:nE
        for a = arms
            switch a{1}
                case 'A',   G = P(i).GA(P(i).ind, :);
                case 'Bp',  G = P(i).GB(P(i).ind, :);
                case 'Sp',  G = P(i).GS(P(i).ind, :);
                case 'Mem', G = repelem(P(i).GMem, P(i).nMem, 1);
            end
            V = nan(size(G, 1), 4);
            for r = 1:size(G, 1)
                d = D(kmap(wkey(P(i).arch, G(r, :))), :);
                V(r, :) = [d(1) d(2) d(3) d(2) / d(1)];
            end
            pooled.(a{1}) = [pooled.(a{1}); V];
            perEv.(a{1})(i, :) = median(V, 1, 'omitnan');
            if strcmp(a{1}, 'Bp'), VB = V; elseif strcmp(a{1}, 'Sp'), VS = V; end
        end
        % paired per replaced individual (same slot: B' target minus S' target), then the event's median
        dif = VB - VS; dif(~isfinite(dif)) = NaN;
        rt = VB ./ VS; rt(~isfinite(rt) | rt <= 0) = NaN; lr = log(rt);   % ratios used for beta and beta/margin only
        pdiff(i, :) = median(dif, 1, 'omitnan'); plog(i, :) = median(lr, 1, 'omitnan');
        PD = [PD; dif]; PL = [PL; lr]; %#ok<AGROW>
    end
    rng(cfg.bootSeed, 'twister'); bi = randi(nE, nE, cfg.nBoot);
    c = struct('names', {names}, 'pooled', pooled, 'perEv', perEv);
    c.pooledDiff = PD; c.pooledLog = PL;
    c.summary = struct();
    for k = 1:4
        for a = arms
            x = perEv.(a{1})(:, k);
            c.summary.(names{k}).(a{1}) = [median(x, 'omitnan'), quantile(x, [0.25 0.75]), sum(isnan(x))];
        end
        c.summary.(names{k}).diff = boot(pdiff(:, k), bi);
        if k == 2 || k == 4, c.summary.(names{k}).logratio = boot(plog(:, k), bi); end
    end
    for a = arms
        V = pooled.(a{1});
        c.shareCretPos.(a{1}) = mean(V(~isnan(V(:, 3)), 3) > 0);
        c.nCretNaN.(a{1}) = sum(isnan(V(:, 3)));
        c.nInd.(a{1}) = size(V, 1);
    end
    c.arms = arms;
end


function r = boot(dd, bi)
% Paired differences over populations: median and mean, with 95% percentile bootstrap intervals.
% With missing values, the resamples are drawn again, from the current stream, over those present.
    ok = ~isnan(dd);
    if ~all(ok), dd = dd(ok); bi = randi(numel(dd), numel(dd), size(bi, 2)); end
    bm = median(dd(bi), 1); bu = mean(dd(bi), 1);
    r = struct('median', median(dd), 'medianCI', quantile(bm, [0.025 0.975]), ...
        'mean', mean(dd), 'meanCI', quantile(bu, [0.025 0.975]), 'n', numel(dd), 'nPos', sum(dd > 0), 'nNeg', sum(dd < 0));
end


% =========================================================================
function c = check2(P, cfg)
    nE = numel(P);
    pops = {'A', 'Bp', 'Sp'}; fld = struct('A', {{'GA', 'GAi'}}, 'Bp', {{'GB', 'GBi'}}, 'Sp', {{'GS', 'GSi'}});
    c = struct();
    % positive control for the counter: a clonal population counts 1, a population of N distinct rows N
    assert(ndist(repmat(P(1).GA(1, :), 1000, 1)) == 1 && ndist([P(1).GA(1:5, :); P(1).GA(1:5, :)]) == ndist(P(1).GA(1:5, :)), 'counter');
    for inv = 0:1
        for a = pops
            f = fld.(a{1}){inv + 1};
            g = zeros(nE, 1); w = g; top = g; eff = g;
            for i = 1:nE
                G = P(i).(f); g(i) = ndist(G);
                W = G(:, 1:16) .* G(:, 17:32); w(i) = ndist(W);
                [~, ~, ig] = unique(G, 'rows'); p = accumarray(ig, 1) / size(G, 1);
                top(i) = max(p); eff(i) = 1 / sum(p .^ 2);
            end
            c.(sprintf('inv%d', inv)).(a{1}) = struct('geno', g, 'wiring', w, 'top', top, 'eff', eff);
        end
    end
    rng(cfg.bootSeed, 'twister'); bi = randi(nE, nE, cfg.nBoot);
    for inv = 0:1
        s = c.(sprintf('inv%d', inv));
        for m = {'geno', 'wiring', 'top', 'eff'}
            for a = pops, c.(sprintf('inv%d', inv)).summary.(m{1}).(a{1}) = [median(s.(a{1}).(m{1})), quantile(s.(a{1}).(m{1}), [0.25 0.75])]; end
            c.(sprintf('inv%d', inv)).summary.(m{1}).diff = boot(s.Bp.(m{1}) - s.Sp.(m{1}), bi);
            c.(sprintf('inv%d', inv)).summary.(m{1}).BminusA = boot(s.Bp.(m{1}) - s.A.(m{1}), bi);
            c.(sprintf('inv%d', inv)).summary.(m{1}).SminusA = boot(s.Sp.(m{1}) - s.A.(m{1}), bi);
        end
    end
end


function n = ndist(G)
    n = size(unique(G, 'rows'), 1);
end


% =========================================================================
function c = check3(P, spArch, cfg)
% Class under the runs' encounter rule of each distinct target wiring and control, cfg.nSeed seeds.
% Seeds: the members sampled as controls from rng(bootSeed + 3, 'twister'); wiring q, seed s, from
% 90000000 + 100 q + s.
    groups = {'Sp', 'Bp', 'A', 'Mem'};
    L = struct('arch', {}, 'W', {}, 'key', {}); kmap = containers.Map();
    use = struct();
    for g = groups, use.(g{1}) = cell(numel(P), 1); end
    rng(cfg.bootSeed + 3, 'twister');
    for i = 1:numel(P)
        for g = groups
            switch g{1}
                case 'A',   G = P(i).GA(P(i).ind, :);
                case 'Bp',  G = P(i).GB(P(i).ind, :);
                case 'Sp',  G = P(i).GS(P(i).ind, :);
                case 'Mem'
                    G = P(i).GMem; if size(G, 1) > cfg.nMemCtl, G = G(randperm(size(G, 1), cfg.nMemCtl), :); end
                    G = [G; zeros(0, 32)];
            end
            [Gu, ~, ig] = unique(G, 'rows'); cnt = accumarray(ig, 1);
            id = zeros(size(Gu, 1), 1);
            for r = 1:size(Gu, 1)
                k = wkey(P(i).arch, Gu(r, :));
                if ~isKey(kmap, k)
                    L(end+1) = struct('arch', P(i).arch, 'W', reshape(Gu(r, 1:16) .* Gu(r, 17:32), 4, 4), 'key', k); %#ok<AGROW>
                    kmap(k) = numel(L);
                end
                id(r) = kmap(k);
            end
            use.(g{1}){i} = [id, cnt];
        end
    end
    nL = numel(L);
    pick = 1:nL;
    fprintf('CHECK 3: %d distinct wirings to test for membership under the runs'' encounter rule, %d seeds each (%d wirings run)\n', nL, cfg.nSeed, numel(pick));
    cls = nan(nL, cfg.nSeed); pc1 = nan(nL, cfg.nSeed);
    t1 = tic; Lp = L(pick); clsP = nan(numel(pick), cfg.nSeed); pc1P = clsP;
    nSeed = cfg.nSeed;
    parfor q = 1:numel(pick)
        cr = nan(1, nSeed); pr = nan(1, nSeed);
        for s = 1:nSeed
            r = class_under_encounter_rule(Lp(q).W, spArch{Lp(q).arch}, 90000000 + 100 * pick(q) + s, true);
            cr(s) = r.cls; pr(s) = r.pC1;
        end
        clsP(q, :) = cr; pc1P(q, :) = pr;
    end
    cls(pick, :) = clsP; pc1(pick, :) = pc1P;
    fprintf('  membership tested in %.1f min\n', toc(t1) / 60);
    c = struct('cls', cls, 'pC1', pc1, 'pick', pick, 'use', use, 'nL', nL);
    % summaries: share retaliating (cls == 1) per group, over distinct wirings and per replaced individual
    for g = groups
        c.summary.(g{1}) = summarise3(cls, pc1, use.(g{1}), pick);
    end
end


function s = summarise3(cls, pc1, use, pick)
% Per distinct wiring: member by majority of seeds (ties = not a member); per individual weights = counts.
    picked = false(size(cls, 1), 1); picked(pick) = true;
    nS = size(cls, 2);
    ids = []; w = []; ev = [];
    for i = 1:numel(use)
        u = use{i}; if isempty(u), continue; end
        keep = picked(u(:, 1));
        ids = [ids; u(keep, 1)]; w = [w; u(keep, 2)]; ev = [ev; repmat(i, sum(keep), 1)]; %#ok<AGROW>
    end
    nMem = sum(cls(ids, :) == 1, 2); nC1 = sum(cls(ids, :) == 2, 2); nF = sum(cls(ids, :) == 3, 2);
    maj = nMem > nS / 2;
    s = struct('nDistinctPerEvent', numel(ids), 'shareWiringMajority', mean(maj), ...
        'shareWiringAnySeed', mean(nMem > 0), 'shareWiringAllSeeds', mean(nMem == nS), ...
        'shareSeedLevel', mean(cls(ids, :) == 1, 'all'), 'shareFailC1', mean(cls(ids, :) == 3, 'all'), ...
        'shareIndMajority', sum(w .* maj) / sum(w), 'shareIndSeedLevel', sum(w .* mean(cls(ids, :) == 1, 2)) / sum(w), ...
        'nInd', sum(w), 'nWiringRows', numel(ids), 'meanPC1', sum(w .* mean(pc1(ids, :), 2)) / sum(w));
    evs = unique(ev); pe = nan(numel(evs), 1);
    for j = 1:numel(evs)
        k = ev == evs(j); pe(j) = sum(w(k) .* mean(cls(ids(k), :) == 1, 2)) / sum(w(k));
    end
    s.perEventShare = [median(pe), quantile(pe, [0.25 0.75])];
    s.perEventAll = pe; s.events = evs;
    s.nC1 = sum(nC1(:)); s.nF = sum(nF(:));
end


% =========================================================================
function write_report(o, file, nE)
    fid = fopen(file, 'w'); cl = onCleanup(@() fclose(fid));
    p = @(varargin) fprintf(fid, varargin{:});
    p('REPAIR FAIRNESS CHECKS (no collapse outcome read)\n');
    p('Input: publication/repair_matched/repair_matched_build_targets.mat (run_repair_test, mode ''build'').\n');
    p('%d included events (populations). No evolution run, no collapse outcome read. Arms: A = the population at the restart point,\n', nE);
    p('B'' = repaired, S'' = sham; the same individuals are replaced in B'' and S''. Bootstrap: 95%% percentile over events, %d resamples, rng(%d).\n\n', o.cfg.nBoot, o.cfg.bootSeed);
    p('WHAT EACH CHECK COULD HAVE SHOWN AND THAT THE CODE CAN SHOW IT\n');
    p('  Check 1: B'' targets with a larger cooperation margin (''margin'' below) or a smaller beta than S'' targets (a robustness confound), or the reverse. The same\n');
    p('    code gives c_ret > 0 for members at the restart point and c_ret < 0 for the replaced cooperators that do not retaliate (A), so it separates classes where they differ; the\n');
    p('    margin and beta columns come from the same membership_test call.\n');
    p('  Check 2: fewer distinct genotypes under one arm (less raw material for selection). Positive control in the code: a clonal population counts 1,\n');
    p('    a duplicated set counts as itself (asserted). The replacement maps many genotypes of cooperators that do not retaliate to few targets, so A > B'', S'' is expected.\n');
    p('  Check 3: sham targets that retaliate under the runs'' encounter rule (the sham would then repair too). The rule returns ''does not retaliate''\n');
    p('    for most replaced wirings that cooperate without retaliating (group A) and ''retaliates'' for most members at the restart point, so both answers occur.\n');
    p('\n');
    c = o.check1;
    p('=== CHECK 1. Robustness of the targets (membership test: membership_test, simulate_pair, fixed start 0.005, 600 steps) ===\n');
    p('Per event: median over its replaced individuals (a target taken by many individuals counts many times). Then median [IQR] over events.\n');
    p('Individuals: A %d, B'' %d, S'' %d (replaced); control: members at the restart point %d (clone-weighted).\n', c.nInd.A, c.nInd.Bp, c.nInd.Sp, c.nInd.Mem);
    unit = {'margin (x_C - x_D vs clone)', 'beta (mutation-sensitivity)', 'c_ret (minus readout vs the worst-case mutant)', 'beta / margin'};
    for k = 1:4
        s = c.summary.(c.names{k});
        p('\n%s\n', unit{k});
        for a = {'A', 'Bp', 'Sp', 'Mem'}
            v = s.(a{1}); p('  %-4s median %8.4f  IQR [%8.4f, %8.4f]  (events with no value: %d)\n', a{1}, v(1), v(2), v(3), v(4));
        end
        d = s.diff;
        p('  paired B'' - S'' (per replaced individual, then the event median): median %+.4f [%+.4f, %+.4f]; mean %+.4f [%+.4f, %+.4f]; events B''>S'' %d, B''<S'' %d of %d\n', ...
            d.median, d.medianCI, d.mean, d.meanCI, d.nPos, d.nNeg, d.n);
        if k == 2 || k == 4
            d = s.logratio;
            p('  paired ln(B''/S''): median %+.4f [%+.4f, %+.4f] (ratio %.3f [%.3f, %.3f]); events B''>S'' %d, B''<S'' %d of %d\n', ...
                d.median, d.medianCI, exp(d.median), exp(d.medianCI), d.nPos, d.nNeg, d.n);
        end
    end
    p('\nPooled over all replaced individuals (not per event): median [IQR]\n');
    for k = 1:4
        p('  %-18s', c.names{k});
        for a = {'A', 'Bp', 'Sp', 'Mem'}
            v = c.pooled.(a{1})(:, k); p(' %s %.4f [%.4f, %.4f];', a{1}, median(v, 'omitnan'), quantile(v, 0.25), quantile(v, 0.75));
        end
        p('\n');
    end
    p('\nShare of individuals with c_ret > 0 (where defined): A %.3f, B'' %.3f, S'' %.3f, members %.3f; individuals with no defecting grid mutant (c_ret undefined): A %d, B'' %d, S'' %d, members %d\n', ...
        c.shareCretPos.A, c.shareCretPos.Bp, c.shareCretPos.Sp, c.shareCretPos.Mem, c.nCretNaN.A, c.nCretNaN.Bp, c.nCretNaN.Sp, c.nCretNaN.Mem);
    p('\nPooled individual-level pairs (B'' target minus S'' target, same replaced individual): quantiles 10/25/50/75/90%%\n');
    for k = 1:4
        v = c.pooledDiff(:, k); v = v(~isnan(v));
        p('  %-18s diff %s;  share |diff| < 1e-3: %.3f', c.names{k}, mat2str(quantile(v, [.1 .25 .5 .75 .9]), 4), mean(abs(v) < 1e-3));
        if k == 2 || k == 4
            l = c.pooledLog(:, k); l = l(~isnan(l));
            p(';  ln ratio %s;  share |ln ratio| < 0.01: %.3f, ratio > 1: %.3f', mat2str(quantile(l, [.1 .25 .5 .75 .9]), 4), mean(abs(l) < 0.01), mean(l > 0));
        end
        p('\n');
    end
    p('(Means of paired beta and beta/margin differences are driven by outliers with a near-zero margin or a very large beta; read medians.)\n');
    d2 = o.check2;
    p('\n=== CHECK 2. Genetic diversity of the population entering generation 1 ===\n');
    for inv = 0:1
        if inv == 0, p('\nWithout the added defector (arms A0, B0'', S0'')\n'); else, p('\nWith the 10 added copies of the defector that later took over, at the runs'' slots (arms A, B'', S'')\n'); end
        s = d2.(sprintf('inv%d', inv)).summary;
        for m = {'geno', 'wiring', 'top', 'eff'}
            lab = struct('geno', 'distinct genotypes (M,B)', 'wiring', 'distinct wirings W = M.*B', 'top', 'share of the commonest genotype', 'eff', 'effective number of genotypes (1/sum p^2)');
            p('  %s\n', lab.(m{1}));
            for a = {'A', 'Bp', 'Sp'}
                v = s.(m{1}).(a{1}); p('    %-3s median %8.3f  IQR [%8.3f, %8.3f]\n', a{1}, v);
            end
            for q = {'diff', 'BminusA', 'SminusA'}
                d = s.(m{1}).(q{1}); nm = struct('diff', 'B'' - S''', 'BminusA', 'B'' - A', 'SminusA', 'S'' - A');
                p('    paired %-7s median %+8.3f [%+8.3f, %+8.3f]; mean %+8.3f [%+8.3f, %+8.3f]; events >0: %d, <0: %d of %d\n', ...
                    nm.(q{1}), d.median, d.medianCI, d.mean, d.meanCI, d.nPos, d.nNeg, d.n);
            end
        end
    end
    if ~isempty(o.check3)
        c3 = o.check3;
        p('\n=== CHECK 3. Retaliation under the runs'' encounter rule (class_under_encounter_rule: U(0,0.01) starts, 300 steps, 20 start pairs) ===\n');
        p('%d distinct wirings in all groups; %d run (%d seeds each).\n', c3.nL, numel(c3.pick), o.cfg.nSeed);
        p('cls 1 = member (cooperates with its clone AND defects back against sampled defecting grid mutants, no exploitation) = retaliates; 2 = cooperates, does not retaliate.\n');
        lab = struct('Sp', 'S'' sham targets', 'Bp', 'B'' repair targets', 'A', 'A: the replaced cooperators that do not retaliate', 'Mem', 'members at the restart point (control, up to 4 distinct per event)');
        for g = {'Sp', 'Bp', 'A', 'Mem'}
            s = c3.summary.(g{1});
            p('\n%s: %d (event, wiring) rows, %d individuals\n', lab.(g{1}), s.nWiringRows, s.nInd);
            p('  share of distinct wirings that retaliate: majority of seeds %.3f; in all seeds %.3f; in at least one seed %.3f; per seed %.3f\n', ...
                s.shareWiringMajority, s.shareWiringAllSeeds, s.shareWiringAnySeed, s.shareSeedLevel);
            p('  per replaced individual: majority %.3f, seed-averaged %.3f; per event median [IQR] %.3f [%.3f, %.3f]\n', ...
                s.shareIndMajority, s.shareIndSeedLevel, s.perEventShare);
            p('  share of draws failing C1 under the runs'' encounter rule (cls 3): %.3f; mean p_C1 %.3f\n', s.shareFailC1, s.meanPC1);
        end
        if isfield(c3.summary.Sp, 'perEventAll')
            ps = c3.summary.Sp; pb = c3.summary.Bp;
            [~, ia, ib] = intersect(ps.events, pb.events);
            dd = pb.perEventAll(ib) - ps.perEventAll(ia);
            rng(o.cfg.bootSeed, 'twister'); bi = randi(numel(dd), numel(dd), o.cfg.nBoot); r = boot(dd, bi);
            p('\nPaired per event, share of replaced individuals whose target retaliates, B'' - S'': median %+.3f [%+.3f, %+.3f]; mean %+.3f [%+.3f, %+.3f]\n', ...
                r.median, r.medianCI, r.mean, r.meanCI);
        end
    end
    p('\nRuntime %.1f min.\n', o.minutes);
end
