function out = measure_retaliation_fixed_opponent(cfg)
% MEASURE_RETALIATION_FIXED_OPPONENT  For each member of the per-mutation census, which of its 30
% census mutations end its retaliation against its worst-case mutant H, held fixed? Also recovers
% each member's wiring and finds H: the input of Fig. 4A.
%
% Produces : main text, "Conditional cooperation stays within a mutation's reach of the
%            threshold"; Methods, "Mutation tests in the game"; SI S3.3, "The search" (members
%            recovered, members with a defecting mutant); Fig. 4A (all as input: the output is
%            what sample_census_members reads for tipping_input_flag). The printed report's
%            rho_ret statistics are not paper numbers.
% Inputs   : publication/rho_by_architecture.mat   (the census; per_mutation_census)
%            publication/retaliation_margin_vs_rho.mat   (row-aligned readouts against the worst-case
%              mutant; retaliation_margin_vs_rho)
%            publication/<architecture>_100reps.mat   (the runs; run_signaling_architectures)
% Outputs  : publication/retaliation_fixed_opponent/retaliation_fixed_opponent_census.mat (out,
%              cfg), .txt (the printed report) and _census_ck<m>.mat (one checkpoint per
%              architecture m, resumed if present); returns out
% Run      : run_fixed_opponent_census   (57.0 min on 4 workers);
%              measure_retaliation_fixed_opponent(struct('reportOnly', true)) reprints the report
%              from the stored .mat (about 15 min)
% Seed     : no draws in the simulations (fixed start); bootstraps and shuffles use
%            rng(20260928 + k, 'twister'), k = 5, 7 and 9 to 14; the bootstraps that repeat the
%            stored census's cooperation interval use rng(0); `cfg.pilot` samples its members with
%            rng(20260928 + 200 + m, 'twister') for architecture m
%
% Scored: every census row marked ok whose wiring is recovered uniquely. H: of the 160 grid mutants
% that defect on the member, the one against which its readout is highest (nearest zero); the
% retaliation margin is minus that readout, and rho_ret is the predicted shift over that margin.
% Members with no defecting mutant have no H. `cfg.pilot` (a trial on replicates 1 to 3, 20
% members per architecture, which prints a timing projection) feeds no paper number.
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
    matFile = fullfile(cfg.outDir, ['retaliation_fixed_opponent' cfg.fileTag '.mat']);
    if cfg.reportOnly
        S = load(matFile, 'out'); out = S.out;
    else
        t0 = tic; out = run_all(cfg); out.minutes = toc(t0) / 60;
        % `cfg.pilot`: replicates 1 to 3 only, 20 members per architecture; the report is
        % generated and discarded, and only a timing projection is printed
        if cfg.pilot
            evalc('print_report(out, cfg)');      % generated, discarded unread
            out.reportMin = timing_run_report_minutes(out, cfg);
            print_timing_run(out, cfg);
            return;
        end
        save(matFile, 'out', 'cfg', '-v7');
    end

    %% Report
    report = @() print_report(out, cfg); %#ok<NASGU> called inside evalc
    rep = evalc('report();');
    fprintf('%s', rep);

    %% Save
    fid = fopen(fullfile(cfg.outDir, ['retaliation_fixed_opponent' cfg.fileTag '.txt']), 'w');
    fprintf(fid, '%s', rep); fclose(fid);
end


% =========================================================================
function cfg = fill_defaults(cfg)
% sigma: the membership test's grid SD and the one-sided slope step (0.2)
    pub = fullfile(grn_data_file('publication'));
    def = struct('dataDir', pub, 'outDir', '', 'rngSeed', 20260928, 'nTimingRunPerArch', 20, ...
        'nBoot', 2000, 'nShuffle', 200, 'nWorkers', 2, 'sigma', 0.2, 'pilot', false, ...
        'reportOnly', false, 'fileTag', '_census', 'fullCensus', true);
    assert(~isfield(cfg, 'fullCensus') || cfg.fullCensus, 'the census (fullCensus true) is the only run');
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
    if isempty(cfg.outDir), cfg.outDir = fullfile(cfg.dataDir, 'retaliation_fixed_opponent'); end
end


function T = arch_table()
    T = {'graded_separate',     'separate',     Inf
         'graded_contrast',     'contrast', Inf
         'saturating_separate', 'separate',     1
         'saturating_contrast', 'contrast', 1};
end


function sp = params(t)
% simulate_pair parameters, as per_mutation_census's sim_params.
    sp = struct('sensingMode', t{2}, 'hillCoefficient', t{3}, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), ...
        'externalSignal2', zeros(4, 1));
end


% =========================================================================
function out = run_all(cfg)
    R = load(fullfile(cfg.dataDir, 'rho_by_architecture.mat'), 'W_all', 'modeNames', 'metadata');
    Q = load(fullfile(cfg.dataDir, 'retaliation_margin_vs_rho.mat'), 'R_all', 'modeNames');
    rc = R.metadata.cfg; modes = arch_table();
    out.modes = R.modeNames; out.K = rc.K;
    out.census = census_reference(R.W_all, cfg);
    if cfg.nWorkers > 1 && isempty(gcp('nocreate')), parpool('local', cfg.nWorkers); end
    for m = 1:4
        assert(strcmp(R.modeNames{m}, modes{m, 1}) && strcmp(Q.modeNames{m}, modes{m, 1}), 'mode order');
        rows = R.W_all{m}; qrows = Q.R_all{m}; R.W_all{m} = []; Q.R_all{m} = [];
        check_alignment(rows, qrows, modes{m, 1});
        ck = fullfile(cfg.outDir, sprintf('retaliation_fixed_opponent_census_ck%d.mat', m));
        if ~cfg.pilot && isfile(ck)                                    % resume from the checkpoint
            S = load(ck, 'A'); out.arch{m} = S.A;
            fprintf('%s: resumed from checkpoint (%d members)\n', S.A.name, numel(S.A.res));
            continue;
        end
        if cfg.pilot
            keep = [rows.rep] <= 3; rows = rows(keep); qrows = qrows(keep);
        end
        sp = params(modes(m, :)); A = struct();
        tr = tic; [Wp, rec] = recover_parents(rows, modes{m, 1}, rc, cfg); A.recMin = toc(tr) / 60;
        A.recovery = rec; A.name = modes{m, 1};
        [sub, info] = pick_sample(rows, Wp, cfg, m);
        A.sample = info;
        mem = member_inputs(rows, qrows, Wp, sub, m); A.memInputs = mem;
        tm = tic; A.res = run_members(mem, sp, cfg); A.memMin = toc(tm) / 60;
        A.selfCheck = self_check(mem, A.res, sp, cfg);
        out.arch{m} = A;
        if ~cfg.pilot, save(ck, 'A', '-v7'); end
        fprintf('%s: %d members, recovery %.1f min, members %.1f min\n', A.name, numel(sub), A.recMin, A.memMin);
    end
end


function check_alignment(rows, qrows, name)
% The two stored files are row-aligned (row counts, ok flags and parent r asserted equal).
    assert(numel(rows) == numel(qrows), '%s: row counts differ', name);
    assert(isequal(logical([rows.ok]), logical([qrows.ok])), '%s: ok flags differ', name);
    assert(max(abs([rows.r] - [qrows.r])) == 0, '%s: parent r differs', name);
end


function [sub, info] = pick_sample(rows, Wp, cfg, m)
% Every 'ok' recovered row; `cfg.pilot`: nTimingRunPerArch of them at random. info keeps the stored
% census's field names (run1 and overlap are empty and 0).
    ok = find([rows.ok]' & all(isfinite(Wp), 2));
    info = struct('run1', [], 'nPool', numel(ok), 'nOk', sum([rows.ok]), 'overlap', 0);
    sub = ok;
    if cfg.pilot
        rng(cfg.rngSeed + 200 + m, 'twister');
        sub = sort(ok(randperm(numel(ok), min(cfg.nTimingRunPerArch, numel(ok)))));
    end
end


function mem = member_inputs(rows, qrows, Wp, sub, m)
    mem = struct('W', {}, 'arch', {}, 'row', {}, 'rep', {}, 'count', {}, 'r', {}, 'slot', {}, ...
        'dW', {}, 'isToggle', {}, 'rho_dir', {}, 'lost', {}, 'ejected', {}, 'wr0', {}, 'wr_m', {});
    for q = 1:numel(sub)
        i = sub(q); a = rows(i); b = qrows(i);
        mem(q) = struct('W', reshape(Wp(i, :), 4, 4), 'arch', m, 'row', i, 'rep', a.rep, ...
            'count', a.count, 'r', a.r, 'slot', a.slot(:), 'dW', a.dW(:), 'isToggle', a.isToggle(:), ...
            'rho_dir', a.rho_dir(:), 'lost', logical(a.lostReadout(:)), 'ejected', logical(a.ejected(:)), ...
            'wr0', b.wr0, 'wr_m', b.wr_m(:));
    end
end


function res = run_members(mem, sp, cfg)
    n = numel(mem); cells = cell(n, 1);
    if cfg.nWorkers > 1
        parfor q = 1:n, cells{q} = member_metrics(mem(q), sp, cfg); end
    else
        for q = 1:n, cells{q} = member_metrics(mem(q), sp, cfg); end
    end
    res = vertcat(cells{:});
end


function sc = self_check(mem, res, sp, cfg)
% Control (d)(i), the last item in the list of print_controls_census: the instrument run twice on the same
% members must agree exactly.
    n = min(3, numel(mem)); sc.n = n; sc.identical = true;
    for q = 1:n
        again = member_metrics(mem(q), sp, cfg);
        a = rmfield(res(q), 'secs'); b = rmfield(again, 'secs');
        sc.identical = sc.identical && isequaln(a, b);
    end
end


% =========================================================================
% RECOVERY of each census row's parent wiring (as measure_joint_loss.m): among the distinct
% genotypes of the row's snapshot with the row's clone count (searched after the previous row's
% match), the one whose clone readout equals the stored r within 1e-9; ties broken by consistency with the stored draws, then by the
% membership test; still tied or no match: not recovered (counted)
% =========================================================================
function [Wp, rec] = recover_parents(rows, mode, rc, cfg)
    d = load(fullfile(cfg.dataDir, sprintf('%s_100reps.mat', mode)), 'results');
    sp = params_by_name(mode);
    nR = numel(rows); Wp = nan(nR, 16);
    rep = [rows.rep]; gen = [rows.gen];
    nAmb = 0; nMiss = 0; nTie = 0; maxDr = 0;
    for r = unique(rep)
        R = d.results{r};
        sg = double(R.historyGenerations(:))';
        ks = find(sg >= rc.steadyStateFrom);
        ks = ks(round(linspace(1, numel(ks), rc.nSnaps)));
        for k = ks
            idx = find(rep == r & gen == sg(k));
            if isempty(idx), continue; end
            M = double(R.historyM(:,:,:,k)); B = double(R.historyB(:,:,:,k)); N = size(M, 3);
            [Wu, ~, ic] = unique(reshape(M .* B, 16, N)', 'rows');
            counts = accumarray(ic, 1);
            last = 0;
            for i = idx
                cand = find(counts == rows(i).count & (1:numel(counts))' > last);
                dr = inf(numel(cand), 1);
                for c = 1:numel(cand)
                    dr(c) = abs(clone_r(Wu(cand(c), :), sp) - rows(i).r);
                end
                hit = cand(dr < 1e-9);
                if isempty(hit), nMiss = nMiss + 1; continue; end
                maxDr = max(maxDr, min(dr(dr < 1e-9)));
                if numel(hit) > 1
                    nTie = nTie + 1;
                    hit = hit(arrayfun(@(h) draws_consistent(Wu(h, :), rows(i)), hit));
                    if numel(hit) > 1
                        hit = hit(arrayfun(@(h) is_member_emu(Wu(h, :), sp), hit));
                    end
                    if numel(hit) ~= 1
                        nAmb = nAmb + 1; last = min([cand(dr < 1e-9); Inf]); continue;
                    end
                end
                Wp(i, :) = Wu(hit, :); last = hit;
            end
        end
    end
    rec = struct('n', nR, 'missing', nMiss, 'ambiguous', nAmb, 'ties', nTie, 'maxDr', maxDr, ...
        'toggleBad', toggle_check(rows, Wp));
end


function ok = draws_consistent(w, row)
    s = row.slot(:); t = row.isToggle(:); dW = row.dW(:); ws = w(s)'; ws = ws(:);
    okT = ~t | ws == 0 | abs(dW + ws) < 1e-12;
    okW = t | ((ws == 0) == (dW == 0));
    ok = all(okT & okW);
end


function m = is_member_emu(w, sp)
    W = reshape(w, 4, 4); D = grid_deltas_joint(); nN = size(D, 1); X0 = 0.005 * ones(4, 1);
    Wn = repmat(W, 1, 1, nN);
    for q = 1:nN, Wn(D(q,1), D(q,2), q) = W(D(q,1), D(q,2)) + D(q,3); end
    [a1, a2] = encounter_batch(repmat(W, 1, 1, nN + 1), cat(3, W, Wn), repmat(X0, 1, nN + 1), ...
        repmat(X0, 1, nN + 1), sp, 600, 'simpair');
    m = a1(1) && a2(1) && ~any(a1(2:end) & ~a2(2:end));
end


function D = grid_deltas_joint()
    d = linspace(-0.6, 0.6, 11); d = d(d ~= 0);
    [dd, jj, ii] = ndgrid(d, 1:4, 1:4);
    D = [ii(:), jj(:), dd(:)];
end


function sp = params_by_name(mode)
    T = arch_table(); sp = params(T(strcmp(T(:, 1), mode), :));
end


function r = clone_r(w, sp)
    W = reshape(w, 4, 4); E = ones(4);
    [~, fs, ~] = simulate_pair(E, W, [], E, W, [], sp);
    r = fs(1) - fs(2);
end


function nBad = toggle_check(rows, Wp)
    nBad = 0;
    for i = 1:numel(rows)
        if any(isnan(Wp(i, :))), continue; end
        t = rows(i).isToggle; s = rows(i).slot; dW = rows(i).dW;
        live = t & Wp(i, s(:))' ~= 0;
        nBad = nBad + sum(abs(dW(live) + Wp(i, s(live))') > 1e-12);
    end
end


% =========================================================================
% THE INSTRUMENT (one member, its 30 census mutations)
% Every pair is simulate_pair from the fixed start, 600 steps; readout = C - D of agent 1.
% Per member M (clone readout r): its worst-case mutant H on the grid (16 weights x 10 nonzero
% moves of linspace(-0.6, 0.6, 11)) and the retaliation margin c_ret = -readout(M vs H).
% Per mutation M -> M' (dl = W'(slot) - W(slot)): the readout against H held fixed; the loss with
% the opponent fixed is readout(M' vs H) > 0; the predicted shift is g dl, g the one-sided slope of
% readout(. vs H) at M in that slot measured at +/- 0.2 on the side dl went, and rho_ret = g dl /
% c_ret; an analytic line (v_i x_j dl at the coupled fixed point) is carried beside. M' is also
% scored for cooperation with its clone and, if it cooperates, on its own full grid.
% =========================================================================
function o = member_metrics(mem, sp, cfg)
    t0 = tic; E = ones(4); W = mem.W; K = numel(mem.dW); sg = cfg.sigma;
    deltas = linspace(-3 * cfg.sigma, 3 * cfg.sigma, 11);   % the membership test's grid
    [~, fc] = simulate_pair(E, W, [], E, W, [], sp); r = fc(1) - fc(2);
    [gF, gB] = clone_slopes(W, r, sp, sg);
    G0 = grid_scan(W, sp, deltas);
    o = blank_output(K); o.r = r; o.wr0 = G0.wr; o.nDef = G0.nDef;
    o.hSlot = G0.slot; o.hDelta = G0.delta;
    hasH = isfinite(G0.wr);
    if hasH
        H = W; H(G0.slot) = H(G0.slot) + G0.delta;
        [hF, hB] = opp_slopes(W, H, G0.wr, sp, sg);
        [~, fH] = simulate_pair(E, W, [], E, H, [], sp);
        [v, u1] = coupled_v(W, H, fH, sp); x1 = fH(1:4);
        o.cret = -G0.wr;
    end
    for t = 1:K
        q = mem.slot(t); Wm = W; Wm(q) = W(q) + mem.dW(t); dl = Wm(q) - W(q);
        if dl >= 0, g = gF(q); else, g = gB(q); end
        o.rho_dir(t) = -(g * dl) / r;
        [~, fm] = simulate_pair(E, Wm, [], E, Wm, [], sp); rm = fm(1) - fm(2);
        o.lost(t) = ~(rm > 0);
        if ~o.lost(t)
            Gm = grid_scan(Wm, sp, deltas);
            o.ejected(t) = Gm.inv; o.wr_m(t) = Gm.wr; o.wrSlot(t) = Gm.slot; o.wrDelta(t) = Gm.delta;
        else
            o.ejected(t) = true;
        end
        if ~hasH, continue; end
        o = fixed_opponent_draw(o, t, Wm, H, q, dl, hF, hB, v, u1, x1, sp);
        if o.ejected(t) && ~o.lost(t) && ~(o.rH(t) > 0)
            o = new_opponent_draw(o, t, W, Wm, q, dl, sp, sg);
        end
    end
    o.secs = toc(t0);
end


function o = blank_output(K)
    z = nan(K, 1); f = false(K, 1);
    o = struct('r', NaN, 'wr0', NaN, 'nDef', 0, 'hSlot', NaN, 'hDelta', NaN, 'cret', NaN, ...
        'rho_dir', z, 'lost', f, 'ejected', f, 'wr_m', z, 'wrSlot', z, 'wrDelta', z, ...
        'rH', z, 'aH', f, 'measH', z, 'predH', z, 'predHan', z, 'anFallback', f, ...
        'rO', z, 'aO', f, 'predO', z, 'rO2', z, 'aO2', f, 'predO2', z, 'secs', NaN);
end


function o = fixed_opponent_draw(o, t, Wm, H, q, dl, hF, hB, v, u1, x1, sp)
% M' against the SAME H: measured shift, loss, directional-slope and analytic predictions.
    E = ones(4);
    [oc, fh] = simulate_pair(E, Wm, [], E, H, [], sp);
    o.rH(t) = fh(1) - fh(2); o.aH(t) = oc(2) == 'C';
    o.measH(t) = o.rH(t) - o.wr0;
    if dl >= 0, g = hF(q); else, g = hB(q); end
    o.predH(t) = g * dl;
    [i, j] = ind2sub([4 4], q);
    if isfinite(v(i)) && u1(i) > 0
        o.predHan(t) = v(i) * x1(j) * dl;
    else
        o.predHan(t) = g * dl; o.anFallback(t) = true;
    end
end


function o = new_opponent_draw(o, t, W, Wm, q, dl, sp, sg)
% A mutation that ends membership while M' keeps cooperation and still defects back against H
% (it fails against another mutant): its predicted shift
% against M''s worst-case mutant O (and against M's own mutant O2 at the same slot and move),
% computed on M.
    O = Wm; O(o.wrSlot(t)) = O(o.wrSlot(t)) + o.wrDelta(t);
    O2 = W; O2(o.wrSlot(t)) = O2(o.wrSlot(t)) + o.wrDelta(t);
    [o.rO(t), o.aO(t), o.predO(t)] = predicts_against(W, O, q, dl, sp, sg);
    [o.rO2(t), o.aO2(t), o.predO2(t)] = predicts_against(W, O2, q, dl, sp, sg);
end


function [r0, aOpp, pred] = predicts_against(W, O, q, dl, sp, sg)
    E = ones(4);
    [oc, f0] = simulate_pair(E, W, [], E, O, [], sp); r0 = f0(1) - f0(2); aOpp = oc(2) == 'C';
    Ws = W; if dl >= 0, Ws(q) = Ws(q) + sg; else, Ws(q) = Ws(q) - sg; end
    [~, fs] = simulate_pair(E, Ws, [], E, O, [], sp);
    if dl >= 0, g = ((fs(1) - fs(2)) - r0) / sg; else, g = (r0 - (fs(1) - fs(2))) / sg; end
    pred = g * dl;
end


function [gF, gB] = clone_slopes(W, r, sp, sg)
% per_mutation_census's one-sided clone slopes, verbatim.
    E = ones(4); gF = zeros(16, 1); gB = zeros(16, 1);
    for q = 1:16
        Wup = W; Wup(q) = Wup(q) + sg;
        Wdn = W; Wdn(q) = Wdn(q) - sg;
        [~, fu] = simulate_pair(E, Wup, [], E, Wup, [], sp);
        [~, fd] = simulate_pair(E, Wdn, [], E, Wdn, [], sp);
        gF(q) = ((fu(1) - fu(2)) - r) / sg;
        gB(q) = (r - (fd(1) - fd(2))) / sg;
    end
end


function [gF, gB] = opp_slopes(W, H, r0, sp, sg)
% The same one-sided slopes of M's readout, with the opponent H held fixed.
    E = ones(4); gF = zeros(16, 1); gB = zeros(16, 1);
    for q = 1:16
        Wup = W; Wup(q) = Wup(q) + sg;
        Wdn = W; Wdn(q) = Wdn(q) - sg;
        [~, fu] = simulate_pair(E, Wup, [], E, H, [], sp);
        [~, fd] = simulate_pair(E, Wdn, [], E, H, [], sp);
        gF(q) = ((fu(1) - fu(2)) - r0) / sg;
        gB(q) = (r0 - (fd(1) - fd(2))) / sg;
    end
end


function [v, u1] = coupled_v(W1, W2, f, sp)
% Readout sensitivity of agent 1 to an additive input on each of its genes, at the coupled
% fixed point of the pair (agent 2 fixed): v' = e_r' (mu I - D_f J_u)^{-1} D_f, agent-1 block.
    x1 = f(1:4); x2 = f(5:8); K = sp.michaelisConstant; mu = sp.degradationRate;
    [s1, D1] = sense(x2, sp); [s2, D2] = sense(x1, sp);
    u = [W1 * x1 + s1; W2 * x2 + s2]; u1 = u(1:4);
    fp = zeros(8, 1); pos = u > 0; fp(pos) = K ./ (K + u(pos)).^2;
    A = mu * eye(8) - diag(fp) * [W1, D1; D2, W2];
    if rcond(A) < 1e-12, v = nan(4, 1); return; end
    row = ([1 -1 0 0 0 0 0 0] / A) * diag(fp);
    v = row(1:4)';
end


function [s, D] = sense(x, sp)
% What a player's sensors receive from partner state x, and its derivative (simulate_pair).
    ir = sp.sensorIntegrationRate; s = zeros(4, 1); D = zeros(4);
    if isfinite(sp.hillCoefficient)
        assert(sp.hillCoefficient == 1, 'only n = 1 implemented');
        den = x(1) + x(2) + 1e-12; oC = x(1) / den; oD = x(2) / den;
        dO = [(x(2) + 1e-12), -x(1); -x(2), (x(1) + 1e-12)] / den^2;
    else
        oC = x(1); oD = x(2); dO = eye(2);
    end
    if strcmp(stored_value('sensingMode', sp.sensingMode), 'separate')
        s(3) = ir * oC; s(4) = ir * oD; D(3:4, 1:2) = ir * dO;
    else
        sig = ir * (oC - oD); s(3) = sig; s(4) = -sig;
        rw = ir * (dO(1, :) - dO(2, :)); D(3, 1:2) = rw; D(4, 1:2) = -rw;
    end
end


% =========================================================================
% CENSUS REFERENCE (stored data only): the cooperation-side shares below and above rho = 1 and
% their interval, and the half-loss point, pooled and per architecture. Each bootstrap seeds its
% own stream.
% reshape(vertcat(...), K, [])' lays the stored draws out as one row per member, K = 30 columns.
% =========================================================================
function C = census_reference(W_all, cfg)
    RHO = []; LOST = []; w = []; arch = [];
    for m = 1:4
        R = W_all{m}; Rk = R(logical([R.ok])'); K = numel(Rk(1).rho_dir);
        RHO = [RHO; reshape(vertcat(Rk.rho_dir), K, [])'];                 %#ok<AGROW>
        LOST = [LOST; reshape(double(vertcat(Rk.lostReadout)), K, [])'];   %#ok<AGROW>
        w = [w; double([Rk.count])']; arch = [arch; m * ones(numel(Rk), 1)]; %#ok<AGROW>
    end
    C.nWir = numel(w); C.nDraw = numel(RHO);
    C.step = step_stats(RHO, LOST, w, arch, true(size(RHO)), cfg.nBoot, cfg.rngSeed + 5);
    [C.x50, C.x50ci] = half_loss(RHO, LOST, w, arch, 1000, 0, false);
    for m = 1:4                                           % per-architecture census values
        ix = find(arch == m); C.nOkArch(m) = numel(ix);
        C.archStep{m} = step_stats(RHO(ix, :), LOST(ix, :), w(ix), arch(ix), true(numel(ix), size(RHO, 2)), ...
            cfg.nBoot, cfg.rngSeed + 5);
        [C.archX50(m), C.archX50ci(m, :)] = half_loss(RHO(ix, :), LOST(ix, :), w(ix), arch(ix), 1000, 0, false);
    end
end


% =========================================================================
% STATISTICS
% =========================================================================
function s = step_point(X, Y, w, keep)
    ww = repmat(w(:), 1, size(X, 2)); k = keep & isfinite(X) & ww > 0;
    lo = k & X < 1; hi = k & X >= 1;
    s.below = sum(ww(lo) .* Y(lo)) / sum(ww(lo));
    s.above = sum(ww(hi) .* Y(hi)) / sum(ww(hi));
    s.ratio = s.above / s.below; s.nBelow = sum(lo(:)); s.nAbove = sum(hi(:));
end


function S = step_stats(X, Y, w, strata, keep, nBoot, seed)
% Clone-weighted P(Y | X < 1) and P(Y | X >= 1); member bootstrap stratified by architecture.
    ww = repmat(w(:), 1, size(X, 2)); k = keep & isfinite(X) & ww > 0;
    lo = k & X < 1; hi = k & X >= 1; Y = double(Y); Y(~k) = 0;
    P = [sum(ww .* lo, 2), sum(ww .* Y .* lo, 2), sum(ww .* hi, 2), sum(ww .* Y .* hi, 2)];
    S = step_point(X, Y, w, keep);
    rng(seed, 'twister'); bb = nan(nBoot, 3);
    for b = 1:nBoot
        c = resample_counts(strata);
        t = c' * P; bb(b, :) = [t(2) / t(1), t(4) / t(3), (t(4) / t(3)) / (t(2) / t(1))];
    end
    S.belowCI = prctile(bb(:, 1), [2.5 97.5]); S.aboveCI = prctile(bb(:, 2), [2.5 97.5]);
    S.ratioCI = prctile(bb(:, 3), [2.5 97.5]);
end


function c = resample_counts(strata)
    n = numel(strata); c = zeros(n, 1);
    for s = unique(strata(:))'
        ix = find(strata == s); pick = ix(randi(numel(ix), numel(ix), 1));
        c = c + accumarray(pick, 1, [n 1]);
    end
end


function [x50, ci] = half_loss(X, Y, w, strata, nBoot, seed, stratified)
% The census's half-loss point: bins of the fixed grid -2.5:0.125:2.5 with >= 40 raw draws,
% clone-weighted loss share, linear crossing of 0.5; bootstrap over members.
    e = -2.5:0.125:2.5; nb = numel(e) - 1; bC = 0.5 * (e(1:end-1) + e(2:end));
    [nW, K] = size(X); b = discretize(X, e);
    wid = repmat((1:nW)', 1, K); ww = repmat(w(:), 1, K);
    v = ~isnan(b) & isfinite(X) & ww > 0; Y = double(Y);
    Sw = sparse(wid(v), b(v), ww(v), nW, nb);
    Swy = sparse(wid(v), b(v), ww(v) .* Y(v), nW, nb);
    Sn = sparse(wid(v), b(v), 1, nW, nb);
    x50 = cross50(bC, boot_curve(ones(nW, 1), Sw, Swy, Sn));
    rng(seed, 'twister'); bs = nan(nBoot, 1);
    for t = 1:nBoot
        if stratified, cnt = resample_counts(strata);
        else, cnt = accumarray(randi(nW, nW, 1), 1, [nW 1]); end
        bs(t) = cross50(bC, boot_curve(cnt, Sw, Swy, Sn));
    end
    bs = sort(bs(~isnan(bs)));
    if isempty(bs), ci = [NaN NaN]; return; end
    ci = [bs(max(1, floor(0.025 * numel(bs)))), bs(min(numel(bs), ceil(0.975 * numel(bs))))];
end


function p = boot_curve(cnt, Sw, Swy, Sn)
    sw = full(cnt' * Sw)'; swy = full(cnt' * Swy)'; n = full(cnt' * Sn)';
    p = swy ./ max(sw, eps); p(n < 40 | sw <= 0) = NaN;
end


function c = cross50(bC, p)
    c = NaN; ok = find(~isnan(p));
    for i = 1:numel(ok) - 1
        a = ok(i); b = ok(i + 1);
        if p(a) < 0.5 && p(b) >= 0.5
            c = bC(a) + (0.5 - p(a)) * (bC(b) - bC(a)) / (p(b) - p(a)); return;
        end
    end
end


function C = corr_stats(P, Mv, w, strata, nBoot, seed)
% (b): Spearman (unweighted) and Pearson (clone-weighted) of predicted vs measured, per draw.
    ww = repmat(w(:), 1, size(P, 2)); k = isfinite(P) & isfinite(Mv);
    C.n = sum(k(:));
    C.spear = corr(P(k), Mv(k), 'Type', 'Spearman'); C.pear = wpear(P(k), Mv(k), ww(k));
    rng(seed, 'twister'); bs = nan(nBoot, 2);
    for b = 1:nBoot
        c = resample_counts(strata); idx = repelem((1:numel(c))', c);
        Pb = P(idx, :); Mb = Mv(idx, :); wb = ww(idx, :); kb = isfinite(Pb) & isfinite(Mb);
        bs(b, :) = [corr(Pb(kb), Mb(kb), 'Type', 'Spearman'), wpear(Pb(kb), Mb(kb), wb(kb))];
    end
    C.spearCI = prctile(bs(:, 1), [2.5 97.5]); C.pearCI = prctile(bs(:, 2), [2.5 97.5]);
end


function r = wpear(x, y, w)
    w = w / sum(w); mx = sum(w .* x); my = sum(w .* y);
    r = sum(w .* (x - mx) .* (y - my)) / sqrt(sum(w .* (x - mx).^2) * sum(w .* (y - my).^2));
end


function sh = shuffle_control(X, Y, w, keep, nPerm, seed)
% (c): permute X across a member's draws, Y fixed; the step must vanish.
    rng(seed, 'twister'); [n, K] = size(X); d = nan(nPerm, 1);
    for p = 1:nPerm
        Xs = X;
        for i = 1:n, Xs(i, :) = X(i, randperm(K)); end
        s = step_point(Xs, Y, w, keep); d(p) = s.above - s.below;
    end
    sh.mean = mean(d); sh.range = prctile(d, [2.5 97.5]); sh.pass = abs(sh.mean) <= 0.02;
end


% =========================================================================
% ASSEMBLY OF THE DRAW MATRICES
% =========================================================================
function D = assemble(out)
    D = struct('X', [], 'Xan', [], 'Y', [], 'Ystrict', [], 'meas', [], 'pred', [], 'predAn', [], ...
        'lost', [], 'endsMem', [], 'rH', [], 'aH', [], 'w', [], 'arch', [], 'hasH', [], 'rhoC', [], ...
        'lostStored', [], 'endsMembershipStored', [], 'rhoStored', [], 'wr0', [], 'wr0Stored', [], ...
        'wrm', [], 'wrmStored', [], 'rO', [], 'aO', [], 'predO', [], 'rO2', [], 'aO2', [], ...
        'predO2', [], 'r', [], 'rStored', [], 'fallback', [], 'secs', []);
    for m = 1:numel(out.arch)
        A = out.arch{m}; R = A.res(:); mem = A.memInputs(:);
        c = @(f) cell2mat(arrayfun(@(x) x.(f)(:)', R, 'uni', 0));
        cret = [R.cret]';
        D.X = [D.X; c('predH') ./ cret]; D.Xan = [D.Xan; c('predHan') ./ cret];
        D.Y = [D.Y; c('rH') > 0]; D.Ystrict = [D.Ystrict; c('rH') > 0 & ~c('aH')];
        D.meas = [D.meas; c('measH')]; D.pred = [D.pred; c('predH')]; D.predAn = [D.predAn; c('predHan')];
        D.lost = [D.lost; c('lost')]; D.endsMem = [D.endsMem; c('ejected')];
        D.rH = [D.rH; c('rH')]; D.aH = [D.aH; c('aH')];
        D.rhoC = [D.rhoC; c('rho_dir')]; D.fallback = [D.fallback; c('anFallback')];
        D.rO = [D.rO; c('rO')]; D.aO = [D.aO; c('aO')]; D.predO = [D.predO; c('predO')];
        D.rO2 = [D.rO2; c('rO2')]; D.aO2 = [D.aO2; c('aO2')]; D.predO2 = [D.predO2; c('predO2')];
        D.wrm = [D.wrm; c('wr_m')]; D.wr0 = [D.wr0; [R.wr0]']; D.r = [D.r; [R.r]'];
        D.secs = [D.secs; [R.secs]'];
        D.lostStored = [D.lostStored; cell2mat(arrayfun(@(x) x.lost(:)', mem, 'uni', 0))];
        D.endsMembershipStored = [D.endsMembershipStored; cell2mat(arrayfun(@(x) x.ejected(:)', mem, 'uni', 0))];
        D.rhoStored = [D.rhoStored; cell2mat(arrayfun(@(x) x.rho_dir(:)', mem, 'uni', 0))];
        D.wrmStored = [D.wrmStored; cell2mat(arrayfun(@(x) x.wr_m(:)', mem, 'uni', 0))];
        D.wr0Stored = [D.wr0Stored; [mem.wr0]']; D.rStored = [D.rStored; [mem.r]'];
        D.w = [D.w; double([mem.count])']; D.arch = [D.arch; m * ones(numel(R), 1)];
        D.hasH = [D.hasH; isfinite(cret)];
    end
    D.Y = double(D.Y); D.Ystrict = double(D.Ystrict);
    lg = {'lost', 'endsMem', 'aH', 'aO', 'aO2', 'fallback', 'lostStored', 'endsMembershipStored', 'hasH'};
    for i = 1:numel(lg), D.(lg{i}) = logical(D.(lg{i})); end
end


% =========================================================================
% REPORT
% =========================================================================
function print_timing_run(out, cfg)
% Census trial run: timing projection only.
    n = 0; secs = [];
    for m = 1:4, R = out.arch{m}.res; n = n + numel(R); secs = [secs; [R.secs]']; end %#ok<AGROW>
    recFull = sum(cellfun(@(A) A.recMin, out.arch)) * 100 / 3;
    nT = out.census.nWir; memMin = mean(secs) * nT / cfg.nWorkers / 60;
    fprintf('TIMING RUN: %.2f min total; %d members; %.2f s per member (range %.2f-%.2f)\n', ...
        out.minutes, n, mean(secs), min(secs), max(secs));
    fprintf(['TIMING RUN: projected full census (%d ok rows, %d workers): members %.1f min + recovery %.1f ' ...
        '+ report %.1f (10%% replica x10) = %.1f min (%.2f h)\n'], nT, cfg.nWorkers, memMin, recFull, ...
        out.reportMin, memMin + recFull + out.reportMin, (memMin + recFull + out.reportMin) / 60);
end


function m = timing_run_report_minutes(out, cfg)
% Census trial run: time the report (discarded) on its members replicated to 10% of each
% architecture's census 'ok' rows, and scale x10.
    o2 = out;
    for a = 1:4
        A = out.arch{a}; k = numel(A.res); t = max(1, round(0.1 * out.census.nOkArch(a)));
        ix = repmat(1:k, 1, ceil(t / k)); ix = ix(1:t);
        A.res = A.res(ix); A.memInputs = A.memInputs(ix); o2.arch{a} = A;
    end
    t0 = tic; evalc('print_report(o2, cfg)'); m = 10 * toc(t0) / 60;
end


function print_report(out, cfg)
    D = assemble(out);
    fprintf('RETALIATION CROSSING WITH THE OPPONENT FIXED (analysis/s3_one_mutation/measure_retaliation_fixed_opponent.m)\n');
    fprintf('Full census, nBoot %d, %.1f min.\n', ...
        cfg.nBoot, out.minutes);
    print_controls_census(out, D);
    P = print_main_measure(out, D, cfg);
    print_other_measures(D, cfg);
end


function ctl = print_controls_census(out, D)
% Census controls. (a1) recovery: no
% row unmatched, unresolved
% ties <= 1% per architecture, no toggle mismatch; (a2) recomputed rho_dir (within 1e-9), lostReadout
% and 'ejected' identical to the stored census on >= 99.9% of draws; (a3') the cooperation-side
% shares from this code equal those from the stored values to 1e-9; (a4') the cooperation
% half-loss is 0.903 [0.871, 0.994] on the full 'ok' set (else within 0.02 of 0.903); (d)(i) the
% instrument run twice agrees.
% Stored field names: 'ejected' = the mutation ends membership, 'lostReadout' = ends cooperation.
    fprintf('\n[0] CONTROLS (full census)\n');
    okRec = true; okSelf = true;
    for m = 1:numel(out.arch)
        A = out.arch{m}; r = A.recovery; s = A.sample; nS = numel(A.res); nH = sum(isfinite([A.res.cret]));
        o1 = r.missing == 0 && r.ambiguous <= 0.01 * r.n && r.toggleBad == 0;
        okRec = okRec && o1; okSelf = okSelf && A.selfCheck.identical;
        fprintf('  (a1) %-20s recovered %d/%d (ties %d, unresolved %d), toggle mismatches %d -> %s | (d)(i) rerun of %d members identical: %d\n', ...
            A.name, r.n - r.missing - r.ambiguous, r.n, r.ties, r.ambiguous, r.toggleBad, pf(o1), ...
            A.selfCheck.n, A.selfCheck.identical);
        fprintf('  (rows) %-20s census ok rows %d; excluded unrecovered (tie/missing) %d; scored %d; no defecting mutant %d; primary members %d\n', ...
            A.name, s.nOk, s.nOk - nS, nS, nS - nH, nH);
    end
    rel = abs(D.rhoC - D.rhoStored) ./ max(1, abs(D.rhoStored));
    same = rel <= 1e-9 & D.lost == D.lostStored & D.endsMem == D.endsMembershipStored; fr = mean(same(:));
    nL = sum(D.lost(:) ~= D.lostStored(:)); nE = sum(D.endsMem(:) ~= D.endsMembershipStored(:));
    fin = isfinite(D.wrm) & isfinite(D.wrmStored); o2 = fr >= 0.999;
    fprintf(['  (a2) replay on %d members x %d draws: rho_dir and the ends-cooperation and ends-membership flags all identical on %.5f of draws ' ...
        '(rho within 1e-9 on %.5f, max rel %.2e; ends-cooperation mismatches %d; ends-membership mismatches %d) -> %s\n'], ...
        numel(D.w), size(D.X, 2), fr, mean(rel(:) <= 1e-9), max(rel(:)), nL, nE, pf(o2));
    fprintf('       reported: parent wr0 identical %d; parent r max|d| %.1e; mutant wr_m max|d| %.1e over %d\n', ...
        isequaln(D.wr0, D.wr0Stored), max(abs(D.r - D.rStored)), max(abs(D.wrm(fin) - D.wrmStored(fin))), sum(fin(:)));
    all1 = true(size(D.rhoC));
    sT = step_point(D.rhoC, double(D.lost), D.w, all1);
    sC = step_point(D.rhoStored, double(D.lostStored), D.w, all1);
    o3 = abs(sT.below - sC.below) <= 1e-9 && abs(sT.above - sC.above) <= 1e-9;
    fprintf('  (a3'') cooperation shares, this code vs stored census, same %d members: below %.10f vs %.10f, above %.10f vs %.10f (|d| %.1e, %.1e) -> %s\n', ...
        numel(D.w), sT.below, sC.below, sT.above, sC.above, abs(sT.below - sC.below), abs(sT.above - sC.above), pf(o3));
    C = out.census; [x50, ci] = half_loss(D.rhoC, double(D.lost), D.w, D.arch, 1000, 0, false);
    sameSet = numel(D.w) == C.nWir;
    if sameSet
        o4 = isequal(round([x50 ci] * 1000), [903 871 994]);
    else
        o4 = abs(x50 - 0.903) <= 0.02;
    end
    fprintf('  (a4'') cooperation half-loss, this code, census bootstrap (rng(0), 1000, unstratified), %d included of %d census ok: %.3f [%.3f, %.3f]; rule: %s -> %s\n', ...
        numel(D.w), C.nWir, x50, ci, ternary(sameSet, 'same set, must be 0.903 [0.871, 0.994]', ...
        'sets differ, point within 0.02 of 0.903'), pf(o4));
    fprintf('       census by this function (all ok rows, stored values): %.3f [%.3f, %.3f] (stored census 0.903 [0.871, 0.994])\n', C.x50, C.x50ci);
    ctl.a = okRec && o2 && o3 && o4; ctl.d = okSelf;
end


function P = print_main_measure(out, D, cfg)
    h = logical(D.hasH); X = D.X(h, :); Y = D.Y(h, :); w = D.w(h); ar = D.arch(h);
    keep = true(size(X)); keepC1 = ~D.lost(h, :);
    fprintf('\n[1] Fixed-opponent loss = readout(M'' vs H) > 0 against rho_ret = slope-predicted shift / c_ret\n');
    fprintf('    members with a defecting mutant %d of %d (no H: %d); draws %d; clone-weighted; 95%% member bootstrap by architecture\n', ...
        sum(h), numel(h), sum(~h), numel(X));
    P = step_stats(X, Y, w, ar, keep, cfg.nBoot, cfg.rngSeed + 7);
    print_step('all draws (slope at +/- 0.2)', P);
    print_step('cooperation-keeping draws only', step_stats(X, Y, w, ar, keepC1, cfg.nBoot, cfg.rngSeed + 7));
    print_step('unweighted', step_stats(X, Y, ones(size(w)), ar, keep, cfg.nBoot, cfg.rngSeed + 7));
    print_step('strict loss (H still defects)', step_stats(X, D.Ystrict(h, :), w, ar, keep, cfg.nBoot, cfg.rngSeed + 7));
    Pan = step_stats(D.Xan(h, :), Y, w, ar, keep, cfg.nBoot, cfg.rngSeed + 7);
    print_step('ROBUSTNESS analytic v_i x_j dl', Pan);
    fprintf('    analytic line used the measured slope on %.4f of draws (rectified-off target or singular)\n', mean(D.fallback(h, :), 'all'));
    Pan.excess = excess_stats(D.Xan(h, :), Y, w, ar, keep, cfg);
    fprintf('    analytic excess step %+.4f [%+.4f, %+.4f]\n', Pan.excess.point, Pan.excess.ci);
    for m = 1:4
        s = ar == m;
        print_step(sprintf('  %s', out.arch{m}.name), step_stats(X(s, :), Y(s, :), w(s), ar(s), keep(s, :), cfg.nBoot, cfg.rngSeed + 7));
    end
    [P.x50, P.x50ci] = half_loss(X, Y, w, ar, cfg.nBoot, cfg.rngSeed + 9, true);
    fprintf('  half-loss point of rho_ret: %.3f [%.3f, %.3f]\n', P.x50, P.x50ci);
    [xa, xaci] = half_loss(D.Xan(h, :), Y, w, ar, cfg.nBoot, cfg.rngSeed + 9, true);
    fprintf('  half-loss point, analytic line: %.3f [%.3f, %.3f]\n', xa, xaci);
    print_bins(X, Y, w);
    P.corr = corr_stats(D.pred(h, :), D.meas(h, :), w, ar, cfg.nBoot, cfg.rngSeed + 10);
    Ca = corr_stats(D.predAn(h, :), D.meas(h, :), w, ar, cfg.nBoot, cfg.rngSeed + 10);
    fprintf('\n  (b) predicted vs measured fixed-opponent shift, %d draws: slope-predicted Spearman %.3f [%.3f, %.3f], Pearson(w) %.3f [%.3f, %.3f]%s\n', ...
        P.corr.n, P.corr.spear, P.corr.spearCI, P.corr.pear, P.corr.pearCI, flag_b(P.corr.spear));
    fprintf('      analytic Spearman %.3f [%.3f, %.3f], Pearson(w) %.3f [%.3f, %.3f]%s\n', ...
        Ca.spear, Ca.spearCI, Ca.pear, Ca.pearCI, flag_b(Ca.spear));
    P.shuffle = shuffle_control(X, Y, w, keep, cfg.nShuffle, cfg.rngSeed + 11);
    fprintf('  (c) within-member shuffle of rho_ret (%d perms): mean(above - below) %+.4f [%+.4f, %+.4f] (not a gate)\n', ...
        cfg.nShuffle, P.shuffle.mean, P.shuffle.range);
    P.excess = excess_stats(X, Y, w, ar, keep, cfg);
    P.codeOK = abs(P.excess.nullMean - P.shuffle.mean) <= 1e-12;
    fprintf('  code check: shuffle mean by excess_stats %+.6f vs shuffle_control %+.6f, |d| %.1e -> %s\n', ...
        P.excess.nullMean, P.shuffle.mean, abs(P.excess.nullMean - P.shuffle.mean), pf(P.codeOK));
    fprintf('  EXCESS STEP: (above - below) %+.4f - shuffle mean %+.4f = %+.4f [%+.4f, %+.4f]\n', ...
        P.excess.step, P.excess.nullMean, P.excess.point, P.excess.ci);
    P.paired = paired_stats(X, Y, w, ar, keep, cfg);
    fprintf('  PAIRED WITHIN-MEMBER STEP (reported): %+.4f [%+.4f, %+.4f] over %d members with draws on both sides of 1 (of %d)\n', ...
        P.paired.point, P.paired.ci, P.paired.n, numel(w));
    print_side_by_side(out, D, X, Y, w, ar, P, cfg);
end


function print_side_by_side(out, D, X, Y, w, ar, P, cfg)
% Retaliation side (opponent fixed, this run) beside the cooperation-side census (stored values,
% all 'ok' rows) and the cooperation side on the included members, pooled and per architecture.
    C = out.census; all1 = true(size(D.rhoC));
    Sc = step_stats(D.rhoC, double(D.lost), D.w, D.arch, all1, cfg.nBoot, cfg.rngSeed + 5);
    [xc, xcci] = half_loss(D.rhoC, double(D.lost), D.w, D.arch, 1000, 0, false);
    fprintf('\n  SIDE BY SIDE: below | above | ratio | half-loss, 95%% CIs\n');
    sbs_row('RETALIATION, opponent fixed (pooled)', P, P.x50, P.x50ci);
    sbs_row('cooperation, stored census (pooled)', C.step, C.x50, C.x50ci);
    sbs_row('cooperation, this code, included members', Sc, xc, xcci);
    for m = 1:4
        s = ar == m;
        Sr = step_stats(X(s, :), Y(s, :), w(s), ar(s), true(sum(s), size(X, 2)), cfg.nBoot, cfg.rngSeed + 7);
        [xr, xrci] = half_loss(X(s, :), Y(s, :), w(s), ar(s), cfg.nBoot, cfg.rngSeed + 9, true);
        sbs_row(sprintf('RETALIATION  %s', out.arch{m}.name), Sr, xr, xrci);
        sbs_row(sprintf('cooperation  %s (census)', out.arch{m}.name), C.archStep{m}, C.archX50(m), C.archX50ci(m, :));
    end
end


function sbs_row(lab, S, x, xci)
    fprintf('  %-44s %.4f [%.4f, %.4f] | %.4f [%.4f, %.4f] | %6.2f [%.2f, %.2f] | %.3f [%.3f, %.3f]\n', ...
        lab, S.below, S.belowCI, S.above, S.aboveCI, S.ratio, S.ratioCI, x, xci);
end


function Q = step_parts(X, Y, ww, keep)
% Per-member clone-weighted sums [w(<1), wY(<1), w(>=1), wY(>=1)], pooled as step_point pools.
    k = keep & isfinite(X) & ww > 0; lo = k & X < 1; hi = k & X >= 1;
    Q = [sum(ww .* lo, 2), sum(ww .* Y .* lo, 2), sum(ww .* hi, 2), sum(ww .* Y .* hi, 2)];
end


function d = dstep(t)
% above - below from pooled sums, one column per replicate.
    d = t(4, :) ./ t(3, :) - t(2, :) ./ t(1, :);
end


function E = excess_stats(X, Y, w, strata, keep, cfg)
% Within-member excess step: (above - below) minus the mean of the same under within-member
% shuffles of X (Y fixed; shuffle_control's seed and draw order); member bootstrap by stratum,
% each member's permutations held fixed.
    [n, K] = size(X); nP = cfg.nShuffle; ww = repmat(w(:), 1, K); Y = double(Y);
    P0 = step_parts(X, Y, ww, keep); Pp = zeros(n, 4, nP);
    rng(cfg.rngSeed + 11, 'twister');
    for p = 1:nP
        Xs = X;
        for i = 1:n, Xs(i, :) = X(i, randperm(K)); end
        Pp(:, :, p) = step_parts(Xs, Y, ww, keep);
    end
    Pr = reshape(Pp, n, 4 * nP);
    E.step = dstep(sum(P0, 1)'); E.nullMean = mean(dstep(reshape(sum(Pr, 1), 4, nP)));
    E.point = E.step - E.nullMean;
    rng(cfg.rngSeed + 13, 'twister'); bb = nan(cfg.nBoot, 1);
    for b = 1:cfg.nBoot
        c = resample_counts(strata);
        bb(b) = dstep((c' * P0)') - mean(dstep(reshape(c' * Pr, 4, nP)));
    end
    E.ci = prctile(bb, [2.5 97.5]);
end


function S = paired_stats(X, Y, w, strata, keep, cfg)
% Paired within-member step (reported): clone-weighted mean over members with draws on both
% sides of 1 of (loss share among draws >= 1) - (loss share among draws < 1).
    k = keep & isfinite(X); lo = k & X < 1; hi = k & X >= 1; Y = double(Y);
    both = any(lo, 2) & any(hi, 2);
    d = sum(Y .* hi, 2) ./ max(sum(hi, 2), 1) - sum(Y .* lo, 2) ./ max(sum(lo, 2), 1);
    d(~both) = 0; ew = w(:) .* both;
    S.n = sum(both); S.point = sum(ew .* d) / sum(ew);
    rng(cfg.rngSeed + 14, 'twister'); bb = nan(cfg.nBoot, 1);
    for b = 1:cfg.nBoot
        c = resample_counts(strata); bb(b) = sum(c .* ew .* d) / sum(c .* ew);
    end
    S.ci = prctile(bb, [2.5 97.5]);
end


function s = flag_b(sp)
    if sp < 0.5, s = '  FLAG: predictor does not track displacement'; else, s = ''; end
end


function print_step(lab, S)
    fprintf('  %-34s below %.4f [%.4f, %.4f] (n %d) | above %.4f [%.4f, %.4f] (n %d) | ratio %.2f [%.2f, %.2f]\n', ...
        lab, S.below, S.belowCI, S.nBelow, S.above, S.aboveCI, S.nAbove, S.ratio, S.ratioCI);
end


function print_bins(X, Y, w)
    ed = [-Inf 0 0 0.25 0.5 0.75 1 1.5 2 Inf];
    lab = {'< 0', '= 0 (dead)', '(0, 0.25)', '[0.25, 0.5)', '[0.5, 0.75)', '[0.75, 1)', '[1, 1.5)', '[1.5, 2)', '>= 2'};
    ww = repmat(w(:), 1, size(X, 2));
    fprintf('  bins of rho_ret: %-12s %8s %10s\n', 'bin', 'n draws', 'P(loss) w');
    for b = 1:numel(lab)
        switch b
            case 1, s = X < 0;
            case 2, s = X == 0;
            case 3, s = X > 0 & X < 0.25;
            otherwise, s = X >= ed(b) & X < ed(b + 1);
        end
        fprintf('                   %-12s %8d %10.4f\n', lab{b}, sum(s(:)), sum(ww(s) .* Y(s)) / sum(ww(s)));
    end
end


function print_other_measures(D, cfg)
    h = logical(D.hasH); w = D.w(h); ar = D.arch(h);
    E = D.endsMem(h, :) & ~D.lost(h, :); same = E & D.rH(h, :) > 0; sameS = same & ~D.aH(h, :);
    fprintf('\n[2] also reported: mutations that end only retaliation (M'' keeps C1, fails the full C2 grid), clone-weighted\n');
    fprintf('    mutations that end only retaliation: %d draws; share of all draws %.4f; of all mutations ending membership %.4f\n', ...
        sum(E(:)), wshare(E, true(size(E)), w), wshare(E, D.endsMem(h, :), w));
    [s1, c1] = boot_share(same, E, w, ar, cfg);
    [s2, c2] = boot_share(sameS, E, w, ar, cfg);
    fprintf('    SAME-OPPONENT (M'' fails against H itself)        %.4f [%.4f, %.4f]; NEW-OPPONENT %.4f\n', s1, c1, 1 - s1);
    fprintf('    strict SAME (M'' cooperates and H still defects)  %.4f [%.4f, %.4f]\n', s2, c2);
    X = D.X(h, :);
    fprintf('    SAME-OPPONENT mutations that end only retaliation, with rho_ret >= 1: %.4f of them\n', wshare(same & X >= 1, same, w));
    newE = E & ~same;
    for k = 1:2
        if k == 1, rO = D.rO(h, :); pr = D.predO(h, :); aO = D.aO(h, :); nm = 'O  = M''s new worst-case mutant';
        else, rO = D.rO2(h, :); pr = D.predO2(h, :); aO = D.aO2(h, :); nm = 'O2 = M''s own single-weight mutant, same slot/delta'; end
        pre = rO < 0; rhoN = pr ./ -rO;
        fprintf('    NEW-OPPONENT, %s: M already fails to defect back (readout >= 0) %.4f (of which O cooperates %.4f); rho_new >= 1 %.4f; rho_new < 1 %.4f; n %d\n', ...
            nm, wshare(newE & ~pre, newE, w), wshare(newE & ~pre & aO, newE & ~pre, w), ...
            wshare(newE & pre & rhoN >= 1, newE, w), wshare(newE & pre & rhoN < 1, newE, w), sum(newE(:)));
    end
end


function s = wshare(num, den, w)
    ww = repmat(w(:), 1, size(num, 2));
    s = sum(ww(num & den)) / sum(ww(den));
end


function [s, ci] = boot_share(num, den, w, strata, cfg)
    ww = repmat(w(:), 1, size(num, 2));
    P = [sum(ww .* den, 2), sum(ww .* (num & den), 2)];
    s = sum(P(:, 2)) / sum(P(:, 1));
    rng(cfg.rngSeed + 12, 'twister'); bb = nan(cfg.nBoot, 1);
    for b = 1:cfg.nBoot
        c = resample_counts(strata); t = c' * P; bb(b) = t(2) / t(1);
    end
    ci = prctile(bb, [2.5 97.5]);
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
