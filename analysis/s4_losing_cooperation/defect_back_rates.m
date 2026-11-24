function out = defect_back_rates(cfg)
% DEFECT_BACK_RATES  How often do members defect back when a defector defects on them, under the
% runs' encounter rule: against their own defecting single-weight mutants, and against the
% defector that later took over, in the failing populations of the census before failure?
%
% Produces : SI S4.4, members defect back in 0.995 of their own defecting mutants' defections
%            (the pooled 'GRID' row of printed block [1]); prints it itself
% Inputs   : cfg (optional struct; defaults in fill_defaults; 'reportOnly', true reprints the
%              saved report; 'pilot', true runs two populations per architecture, saves nothing)
%            publication/collapse_census/collapse_census.mat   (census_before_failure: per
%              failing population and matched control, the wirings that meet C1 at the restart
%              point, and the defector that later took over)
% Outputs  : publication/later_defector_quality/later_defector_quality.mat and .txt (the printed report)
% Run      : defect_back_rates()   (about 1 min)
% Seed     : all 'twister': rng(20260927 + 1000 a + i) for failing population i of architecture
%            a, rng(20260927 + 1000 a + 500 + i) for control i (the start states);
%            rng(20260927 + 99) for the membership_test sample; rng(20260927 + 7), the bootstrap
%
% Members: a copy of the membership test (encounter_batch, start 0.005, 600 steps, 160
% single-weight mutants), checked against membership_test itself on 60 sampled wirings. Each
% member wiring meets each opponent in K = 20 encounters under the runs' encounter rule, member
% first, the same start pairs for all its opponents. Per population, weighted by individuals:
% defecting back ret = P(member D, opponent D) / P(opponent D); exploited E = P(member C,
% opponent D). Opponents: 'LATER DEFECTOR', the defector that later took over; 'GRID', the member's own
% grid mutants that defect on it from the test's start; 'CTRL', in the matched controls, that
% population's most common defecting wiring. Populations with no member or no defection to
% answer drop out.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
    matFile = fullfile(cfg.outDir, 'later_defector_quality.mat');
    if cfg.reportOnly
        S = load(matFile, 'out'); out = S.out;
    else
        t0 = tic; out = run_all(cfg); out.minutes = toc(t0) / 60;
        if cfg.pilot
            evalc('print_report(out, cfg)'); fprintf('timing run: %.1f min, report code ran\n', out.minutes);
            return;
        end
        save(matFile, 'out', 'cfg');
    end
    txt = fullfile(cfg.outDir, 'later_defector_quality.txt');
    if exist(txt, 'file'), delete(txt); end
    diary(txt); print_report(out, cfg); diary off;
end


% =========================================================================
function cfg = fill_defaults(cfg)
    pub = fullfile(grn_data_file('publication'));
    def = struct('dataDir', pub, 'outDir', '', 'rngSeed', 20260927, 'K', 20, 'nBoot', 2000, ...
        'nCssCheck', 60, 'pilot', false, 'reportOnly', false);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'later_defector_quality'); end
end


% =========================================================================
function out = run_all(cfg)
    L = load(fullfile(cfg.dataDir, 'collapse_census', 'collapse_census.mat'), 'out');
    C = L.out; D = grid_deltas();
    out.arch = {C.name};
    out.inv = []; out.ctrl = [];
    for a = 1:numel(C)
        nC = numel(C(a).coll); nK = numel(C(a).ctrl);
        if cfg.pilot, nC = min(nC, 2); nK = min(nK, 2); end
        for i = 1:nC
            rng(cfg.rngSeed + 1000 * a + i, 'twister');
            e = analyse_event(C(a).coll(i), C(a), D, cfg, true);
            e.arch = a; e.idx = i;
            out.inv = [out.inv, e];
        end
        for i = 1:nK
            rng(cfg.rngSeed + 1000 * a + 500 + i, 'twister');
            e = analyse_event(C(a).ctrl(i), C(a), D, cfg, false);
            e.arch = a; e.idx = i;
            out.ctrl = [out.ctrl, e];
        end
        fprintf('%s done\n', C(a).name);
    end
    out.cssCheck = css_check(out, C, cfg);
end


% =========================================================================
function e = analyse_event(rec, oa, D, cfg, withGrid)
% One population: test its wirings for membership before the fall, then its members against the defector that
% later took over (and against their own grid mutants that defect on them).
    e = empty_event();
    e.found = rec.inv.found;
    e.shareRec = rec.share(1);
    W = rec.resG(:, 1:16) .* rec.resG(:, 17:32);
    [Wu, ~, iw] = unique(W, 'rows'); cnt = accumarray(iw, rec.resCnt);
    X0 = 0.005 * ones(4, 1); nU = size(Wu, 1);
    isMem = false(nU, 1); clsDef = cell(nU, 1);
    for u = 1:nU
        F = grid_pairs(reshape(Wu(u, :), 4, 4), oa.sp, 600, X0, X0, D, 'simpair');
        isMem(u) = F.c1 && ~any(F.aR & ~F.aM);
        clsDef{u} = find(~F.aM);
    end
    e.Wres = Wu; e.cntRes = cnt; e.isMem = isMem;
    e.shareEmu = sum(cnt(isMem)) / 1000;
    e.nMemW = sum(isMem); e.nMemInd = sum(cnt(isMem));
    if ~e.found || ~any(isMem), return; end
    WI = reshape(rec.inv.W, 4, 4);
    mu = find(isMem); w = cnt(mu);
    inv = zeros(numel(mu), 4); grd = zeros(numel(mu), 4); nG = zeros(numel(mu), 1);
    for q = 1:numel(mu)
        Wm = reshape(Wu(mu(q), :), 4, 4);
        X1 = rand(4, cfg.K) * 0.01; X2 = rand(4, cfg.K) * 0.01;
        [a1, a2] = encounter_batch(repmat(Wm, 1, 1, cfg.K), repmat(WI, 1, 1, cfg.K), X1, X2, ...
            oa.spE, 300, 'encounter_rule');
        inv(q, :) = outcome_means(a1, a2);
        if withGrid && ~isempty(clsDef{mu(q)})
            nb = clsDef{mu(q)}; nG(q) = numel(nb);
            Wn = single_weight_mutants(Wm, D(nb, :));
            nn = numel(nb);
            [b1, b2] = encounter_batch(repmat(Wm, 1, 1, cfg.K * nn), Wn(:, :, repelem(1:nn, cfg.K)), ...
                repmat(X1, 1, nn), repmat(X2, 1, nn), oa.spE, 300, 'encounter_rule');
            grd(q, :) = outcome_means(b1, b2);
        else
            grd(q, :) = NaN;
        end
    end
    e.inv = event_stats(inv, w);
    e.grid = event_stats(grd, w);
    e.nGridDef = nG;
    e.memberGrid = member_grid_case(rec, Wu, isMem, cnt, WI, inv, mu, oa, cfg);
end


function e = empty_event()
    s = struct('ret', NaN, 'E', NaN, 'pID', NaN, 'pMD', NaN, 'n', 0);
    e = struct('found', false, 'shareRec', NaN, 'shareEmu', NaN, 'nMemW', 0, 'nMemInd', 0, ...
        'Wres', [], 'cntRes', [], 'isMem', [], 'inv', s, 'grid', s, 'nGridDef', [], ...
        'memberGrid', [], 'arch', 0, 'idx', 0);
end


function x = outcome_means(a1, a2)
% [P(member C & opp D), P(member D & opp D), P(member C & opp C), P(member D & opp C)]
    x = [mean(a1 & ~a2), mean(~a1 & ~a2), mean(a1 & a2), mean(~a1 & a2)];
end


function s = event_stats(x, w)
    ok = all(isfinite(x), 2); x = x(ok, :); w = w(ok);
    s = struct('ret', NaN, 'E', NaN, 'pID', NaN, 'pMD', NaN, 'n', sum(ok));
    if ~any(ok), return; end
    m = (w' * x) / sum(w);
    s.pID = m(1) + m(2); s.E = m(1); s.pMD = m(2) + m(4);
    if s.pID > 0, s.ret = m(2) / s.pID; end
end


% =========================================================================
function g = member_grid_case(rec, Wu, isMem, cnt, WI, inv, mu, oa, cfg)
% A population whose defector that later took over differs from a member wiring in one weight by at
% most 0.6 (block [5] of census_before_failure's report): that member's retaliation against it, whether the
% shift is a grid point (a multiple of 0.12), and both actions under the membership test's start.
    g = [];
    if ~(rec.inv.onGrid == 1 && rec.inv.onGridCls == 1), return; end
    dW = abs(Wu - rec.inv.W);
    onG = isMem & sum(dW > 0, 2) == 1 & max(dW, [], 2) <= 0.6 + 1e-9;
    if ~any(onG), g = struct('found', false); return; end
    [~, j] = max(cnt .* onG);
    q = find(mu == j);
    slot = find(dW(j, :) > 0); delta = rec.inv.W(slot) - Wu(j, slot);
    X0 = 0.005 * ones(4, 1);
    [c1, c2] = encounter_batch(reshape(Wu(j, :), 4, 4), WI, X0, X0, oa.sp, 600, 'simpair');
    x = inv(q, :); pid = x(1) + x(2);
    g = struct('found', true, 'slot', slot, 'delta', delta, ...
        'onGridPoint', abs(delta / 0.12 - round(delta / 0.12)) < 1e-9, ...
        'clsMemberC', c1, 'clsDefectorC', c2, 'ret', x(2) / max(pid, eps), 'E', x(1), 'pID', pid, ...
        'memShare', cnt(j) / 1000, 'nOnGridMembers', sum(onG), 'K', cfg.K);
end


% =========================================================================
function chk = css_check(out, C, cfg)
% Check (b): membership_test itself on a seeded sample of the distinct wirings; it must return the
% emulation's member result on every one.
    rng(cfg.rngSeed + 99, 'twister');
    pool = [];
    for arm = 1:2
        if arm == 1, E = out.inv; else, E = out.ctrl; end
        for k = 1:numel(E)
            for u = 1:size(E(k).Wres, 1), pool = [pool; E(k).arch, arm, k, u]; end %#ok<AGROW>
        end
    end
    pick = pool(randperm(size(pool, 1), min(cfg.nCssCheck, size(pool, 1))), :);
    mm = 0; nMem = 0;
    for p = 1:size(pick, 1)
        if pick(p, 2) == 1, E = out.inv(pick(p, 3)); else, E = out.ctrl(pick(p, 3)); end
        r = membership_test(reshape(E.Wres(pick(p, 4), :), 4, 4), C(pick(p, 1)).sp);
        mm = mm + (r.is_css ~= E.isMem(pick(p, 4))); nMem = nMem + r.is_css;
    end
    chk = struct('n', size(pick, 1), 'mismatch', mm, 'nMember', nMem);
end


% =========================================================================
function D = grid_deltas()
% membership_test's grid: linspace(-3*0.2, 3*0.2, 11) with the exact zero skipped.
    d = linspace(-0.6, 0.6, 11); d = d(d ~= 0);
    [dd, jj, ii] = ndgrid(d, 1:4, 1:4);     % membership_test loop order: ii, jj, delta
    D = [ii(:), jj(:), dd(:)];
end


function Wn = single_weight_mutants(W, D)
    nN = size(D, 1); Wn = repmat(W, 1, 1, nN);
    for q = 1:nN, Wn(D(q,1), D(q,2), q) = Wn(D(q,1), D(q,2), q) + D(q,3); end
end


function F = grid_pairs(W, sp, steps, X1, X2, D, snapMode)
% The wiring (slot 1) against its clone and every grid mutant (slot 2), one start pair.
    nN = size(D, 1);
    W2 = cat(3, W, single_weight_mutants(W, D));
    [a1, a2] = encounter_batch(repmat(W, 1, 1, nN + 1), W2, repmat(X1, 1, nN + 1), ...
        repmat(X2, 1, nN + 1), sp, steps, snapMode);
    F.c1 = a1(1) && a2(1); F.aR = a1(2:end); F.aM = a2(2:end);
end


% =========================================================================
% REPORT
% =========================================================================
function print_report(out, cfg)
    rng(cfg.rngSeed + 7, 'twister');
    fprintf('DEFECT-BACK RATES: do members before the fall defect back against the defector that later took over?\n');
    fprintf('K = %d draws per pairing under the runs'' encounter rule; nBoot %d; %.1f min.\n', ...
        cfg.K, cfg.nBoot, out.minutes);
    print_controls(out);
    [Ri, Rg] = print_arms(out, cfg);
    print_member_grid(out);
    print_pooled(Ri, Rg);
end


% Check (a): the emulated member share of every population equals the census's recorded share to
% 1e-12; check (b) above.
function print_controls(out)
    fprintf('\n[0] CONTROLS\n');
    A = [out.inv, out.ctrl];
    d = abs([A.shareEmu] - [A.shareRec]);
    okA = all(d < 1e-12);
    fprintf('  (a) member-share emulation vs census record: max |diff| %.3g over %d records -> %s\n', ...
        max(d), numel(A), pf(okA));
    c = out.cssCheck; okB = c.mismatch == 0;
    fprintf('  (b) membership_test on %d sampled wirings of the populations (%d members): %d mismatches -> %s\n', ...
        c.n, c.nMember, c.mismatch, pf(okB));
end


function [Ri, Rg] = print_arms(out, cfg)
    fprintf('\n[1] MEMBERS vs OPPONENTS, the runs'' encounter rule, member first. Event = unit, clone-weighted within event\n');
    fprintf('    ret = P(member D | opponent D); E = P(member C, opponent D); 95%% bootstrap CI over events\n');
    fprintf('  %-20s %-6s %4s %5s | %24s | %24s | %7s %7s\n', 'architecture', 'arm', 'ev', 'noMem', ...
        'ret [95% CI]', 'E [95% CI]', 'P(opD)', 'P(memD)');
    for a = [1:numel(out.arch), 0]
        for arm = 1:3
            [S, lab, n0] = arm_sel(out, a, arm);
            r = arm_row(S, cfg);
            fprintf('  %-20s %-6s %4d %5d | %.3f [%.3f, %.3f]     | %.3f [%.3f, %.3f]     | %7.3f %7.3f\n', ...
                aname(out, a), lab, r.n, n0, r.ret, r.retci, r.E, r.Eci, r.pID, r.pMD);
            if a == 0 && arm == 1, Ri = r; end
            if a == 0 && arm == 2, Rg = r; end
        end
        E = events_of(out, a, 1);
        dd = arrayfun(@(e) e.inv.ret - e.grid.ret, E); dd = dd(isfinite(dd));
        b = boot_mean(dd, cfg.nBoot);
        fprintf('  %-20s paired ret(LATER DEFECTOR) - ret(GRID): %+.3f [%+.3f, %+.3f] over %d events\n', aname(out, a), b(1), b(2), b(3), numel(dd));
    end
    E = events_of(out, 0, 1);
    fprintf('  members per collapse event: distinct wirings median %g [%g, %g]; member individuals median %g [%g, %g]; grid mutants per member that defect under the membership test, median %g\n', ...
        median([E.nMemW]), min([E.nMemW]), max([E.nMemW]), median([E.nMemInd]), min([E.nMemInd]), max([E.nMemInd]), median(vertcat(E.nGridDef)));
end


function E = events_of(out, a, arm)
% Populations of an arm (1, 2 = the falls; 3 = the matched controls) with a defector that later
% took over.
    if arm == 3, E = out.ctrl; else, E = out.inv; end
    if a ~= 0, E = E([E.arch] == a); end
    E = E([E.found]);
end


function [S, lab, n0] = arm_sel(out, a, arm)
% Per-population statistics of an arm, restricted to populations with at least one member before
% the fall.
    labs = {'LATER DEFECTOR', 'GRID', 'CTRL'}; lab = labs{arm};
    E = events_of(out, a, arm);
    has = [E.nMemW] > 0; n0 = sum(~has); E = E(has);
    if arm == 2, S = [E.grid]; else, S = [E.inv]; end
    if isempty(S), S = struct('ret', {}, 'E', {}, 'pID', {}, 'pMD', {}, 'n', {}); end
end


function r = arm_row(S, cfg)
    r = struct('n', 0, 'ret', NaN, 'retci', [NaN NaN], 'E', NaN, 'Eci', [NaN NaN], 'pID', NaN, 'pMD', NaN);
    if isempty(S), return; end
    [r.ret, r.retci] = mci([S.ret], cfg.nBoot); [r.E, r.Eci] = mci([S.E], cfg.nBoot);
    r.pID = mean([S.pID], 'omitnan'); r.pMD = mean([S.pMD], 'omitnan'); r.n = sum(isfinite([S.ret]));
end


function [m, ci] = mci(x, nBoot)
    b = boot_mean(x, nBoot); m = b(1); ci = b(2:3);
end


function print_member_grid(out)
    fprintf('\n[2] THE MEMBER-GRID CASES (census block [5]: the defector that took over lies within one slot, |dW| <= 0.6, of a member before the fall)\n');
    fprintf('  %-20s %4s | %5s %7s %6s | %-9s | %6s %6s %6s | %s\n', 'architecture', 'ev', 'slot', 'delta', 'onPt', ...
        'cls M/I', 'ret', 'E', 'P(ID)', 'member share, on-grid members');
    for e = out.inv
        g = e.memberGrid;
        if isempty(g), continue; end
        if ~g.found
            fprintf('  %-20s %4d | on-grid member not recovered\n', out.arch{e.arch}, e.idx); continue;
        end
        fprintf('  %-20s %4d | %5d %+7.3f %6d | %s/%s       | %6.3f %6.3f %6.3f | %.3f, %d\n', out.arch{e.arch}, e.idx, ...
            g.slot, g.delta, g.onGridPoint, cd(g.clsMemberC), cd(g.clsDefectorC), g.ret, g.E, g.pID, g.memShare, g.nOnGridMembers);
    end
    fprintf('  (onPt: delta is exactly a grid point, k * 0.12; cls M/I: actions under the membership test''s dynamics, member/defector)\n');
end


function s = cd(c), if c, s = 'C'; else, s = 'D'; end, end


% Pooled, event = unit: check (c), members defect back against their own defecting grid mutants
% in at least 0.95 of encounters (ret(GRID)), and ret(LATER DEFECTOR) beside it.
function print_pooled(Ri, Rg)
    fprintf('\n[3] POOLED (event = unit)\n');
    fprintf('  control (c): ret(GRID) %.3f >= 0.95 -> %s\n', Rg.ret, pf(Rg.ret >= 0.95));
    fprintf('  ret(LATER DEFECTOR) %.3f [%.3f, %.3f] vs ret(GRID) %.3f [%.3f, %.3f]\n', Ri.ret, Ri.retci, Rg.ret, Rg.retci);
end


function s = aname(out, a)
    if a == 0, s = 'POOLED'; else, s = out.arch{a}; end
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end


function b = boot_mean(x, nBoot)
    x = x(isfinite(x)); n = numel(x);
    if n == 0, b = [NaN NaN NaN]; return; end
    m = zeros(nBoot, 1);
    for k = 1:nBoot, m(k) = mean(x(randi(n, n, 1))); end
    q = quantile(m, [0.025 0.975]); b = [mean(x), q(1), q(2)];
end
