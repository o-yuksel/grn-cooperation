function out = safety_search_losses(stage)
% SAFETY_SEARCH_LOSSES  Which of the 1,000 test mutations end cooperation or retaliation, for each
% start and each final wiring of the safety search: the per-network shares of Fig. 4B.
%
% Produces : main text, "Conditional cooperation stays within a mutation's reach of the
%            threshold"; Fig. 4B; SI S3.4; Table S7 (through long_flags.mat, read by
%            fig4_tipping_input, run_safety_search('longreport'), fig4b_retaliation_from_rest and
%            safety_search_cooperate_only)
% Inputs   : publication/safety_search/, all written by run_safety_search: search_setup.mat (the 30
%              starts), search_main_starts_1_4.mat and search_main_starts_5_12.mat (the 300-step searches on starts
%              1 to 12), long_jobs/long_<k>_<S|C>.mat (the 1,000-step searches, and the 100-step
%              need-only-cooperate searches on starts 13 to 30)
% Outputs  : publication/safety_search/long_flags.mat   (per network: group, start, per test
%              mutation whether it ends cooperation (lostC) or retaliation (lostR), and the
%              breaking shares: uS, mutants that are not members cooperating reliably with their
%              clone; uC, mutants not cooperating reliably)
% Run      : safety_search_losses('long')   (4.4 min on 4 workers; the report prints to the console and is not saved)
% Seed     : the test mutations: RandStream('mt19937ar', 'Seed', 99991), the search's own sample;
%            coop_clone_batch's 100 random starts: rng(20261022, 'twister'), caller's restored
%
% Three networks per start: the start (group 1), the final wiring of the search that must also
% retaliate (group 2, 1,000 steps) and of the search that need only cooperate (group 3, 300 steps
% on starts 1 to 12, 100 on the rest). A mutation ends cooperation when the mutant's readout with
% its own clone is at or below zero, and retaliation when the mutant cooperates (readout above
% zero) against the parent's worst-case mutant H; retaliation is not scored for group 3. The
% recomputed breaking shares are checked against the stored ones (the control line); the flags are
% saved only when they pass.

    pub = grn_data_file('publication');
    odir = fullfile(pub, 'safety_search');
    switch stage
        case 'long', out = do_long(odir);
        otherwise, error('safety_search_losses:stage', 'unknown stage %s', stage);
    end
end


function sp = sim_params()
% tipping_input_flag sim_params(4): saturating contrast.
    sp = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), ...
        'externalSignal2', zeros(4, 1));
end


function o = one_net(j, fr)
% One network j against the test mutations fr: its worst-case mutant H (grid_scan), its tipping
% inputs, and per mutation the input it adds (dW, x_j), lostC (the mutant's clone readout <= 0),
% lostR (the mutant still cooperates against H; groups 1 and 2 only), membership and the breaking
% shares uS and uC as run_safety_search's fresh_eval counts them.
    t0 = tic; sp = sim_params(); E = ones(4); W = double(j.M) .* double(j.B); needRet = j.grp ~= 3;
    G = grid_scan(W, sp, linspace(-0.6, 0.6, 11)); o.nDef = G.nDef; o.hSlot = G.slot; o.hDelta = G.delta;
    hasH = isfinite(G.wr); H = W; if hasH, H(G.slot) = H(G.slot) + G.delta; end
    [o.thrC, o.thrR] = thresholds(W, H, sp, needRet && hasH);
    [~, fc] = simulate_pair(E, W, [], E, W, [], sp); xc = fc(1:4); o.rClone = fc(1) - fc(2);
    xh = nan(4, 1); o.rH0 = NaN; if hasH, [~, fh] = simulate_pair(E, W, [], E, H, [], sp); xh = fh(1:4); o.rH0 = fh(1) - fh(2); end
    [Mm, Bm] = apply_muts(j.M, j.B, fr); Wm = double(Mm) .* double(Bm); K = size(Wm, 3);
    o.slot = sub2ind([4 4], fr.g1, fr.g2); o.isToggle = fr.tog;
    o.dW = arrayfun(@(t) Wm(o.slot(t) + 16 * (t - 1)) - W(o.slot(t)), (1:K)');
    o.xjC = xc(fr.g2); o.xjR = xh(fr.g2);
    o.lostC = false(K, 1); o.lostR = false(K, 1);
    for t = 1:K
        [~, fm] = simulate_pair(E, Wm(:, :, t), [], E, Wm(:, :, t), [], sp); o.lostC(t) = ~((fm(1) - fm(2)) > 0);
        if needRet && hasH
            [~, fmh] = simulate_pair(E, Wm(:, :, t), [], E, H, [], sp); o.lostR(t) = (fmh(1) - fmh(2)) > 0;
        end
    end
    [o.mem, info] = member_batch(Wm, true); o.weak = info.weak;
    base = info.weak & info.rClone >= 1e-6; cc = nan(K, 1); cc(base) = coop_clone_batch(Wm(:, :, base));
    o.rel = base & cc >= 0.9; o.uS = mean(~(o.mem & o.rel)); o.uC = mean(~o.rel);
    o.sec = toc(t0);
end


function [thrC, thrR] = thresholds(W, H, sp, doRet)
    thrC = nan(4, 2); thrR = nan(4, 2);
    for i = 1:4
        for sgn = 1:2
            sv = 3 - 2 * sgn;
            thrC(i, sgn) = flip_threshold(@(d) flip_coop(W, sp, i, sv * d));
            if doRet, thrR(i, sgn) = flip_threshold(@(d) flip_ret(W, H, sp, i, sv * d)); end
        end
    end
end


function O = do_long(odir)
% The flags of one_net for every start with a 1,000-step search: group 1 start, 2 the final
% wiring that must also retaliate, 3 the final wiring that need only cooperate.
    s = load(fullfile(odir, 'search_setup.mat'));
    V = [getfield(load(fullfile(odir, 'search_main_starts_1_4.mat'), 'R'), 'R'); ...
         getfield(load(fullfile(odir, 'search_main_starts_5_12.mat'), 'R'), 'R')]; %#ok<GFLD>
    L = {}; d = dir(fullfile(odir, 'long_jobs', 'long_*.mat'));
    for i = 1:numel(d), L{end + 1} = getfield(load(fullfile(d(i).folder, d(i).name)), 'r'); end %#ok<AGROW,GFLD>
    J = struct('grp', {}, 'k', {}, 'M', {}, 'B', {}, 'uS', {}, 'uC', {}, 'step', {});
    for k = 1:numel(s.reps)
        iS = find(cellfun(@(r) r.k == k && r.mode == 'S' && size(r.trace, 1) == 1000, L));
        if isempty(iS), continue; end
        rS = L{iS}; iC = find(cellfun(@(r) r.k == k && r.mode == 'C', L));
        if isempty(iC), rC = V{cellfun(@(r) r.k == k && r.mode == 'C', V)}; else, rC = L{iC}; end
        assert(isequal(rS.fresh.uS(1), rC.fresh.uS(1)) && isequal(rS.fresh.uC(1), rC.fresh.uC(1)));
        J(end + 1) = struct('grp', 1, 'k', k, 'M', s.S0M{k}, 'B', s.S0B{k}, 'uS', rS.fresh.uS(1), 'uC', rS.fresh.uC(1), 'step', 0); %#ok<AGROW>
        J(end + 1) = struct('grp', 2, 'k', k, 'M', rS.M1, 'B', rS.B1, 'uS', rS.fresh.uS(end), 'uC', rS.fresh.uC(end), 'step', size(rS.trace, 1)); %#ok<AGROW>
        J(end + 1) = struct('grp', 3, 'k', k, 'M', rC.M1, 'B', rC.B1, 'uS', rC.fresh.uS(end), 'uC', rC.fresh.uC(end), 'step', size(rC.trace, 1)); %#ok<AGROW>
    end
    fr = draw_muts(RandStream('mt19937ar', 'Seed', 99991), 1000);
    if isempty(gcp('nocreate')), parpool('local', 4); end
    t0 = tic; N = cell(numel(J), 1);
    parfor n = 1:numel(J), N{n} = one_net(J(n), fr); end
    grp = [J.grp]'; kk = [J.k]';
    lostC = cellfun(@(o) mean(o.lostC), N); retOnly = cellfun(@(o) mean(o.lostR & ~o.lostC), N);
    O = struct('J', J, 'N', {N}, 'wallSec', toc(t0), 'reps', s.reps, 'grp', grp, 'k', kk, ...
        'lostC', lostC, 'retOnly', retOnly);
    p = @(varargin) fprintf(varargin{:});
    O.ctrl = control_fresh(N, J, p);
    if ~O.ctrl, error('safety_search_losses:control', 'CONTROL FAILED: long_flags.mat not saved'); end
    nm = {'start', 'S end', 'C end'};
    for g = 1:3
        m = grp == g;
        p('%-8s n %2d: lostC median %.3f [%.3f, %.3f], zeros %d; lostR-only median %.3f [%.3f, %.3f], zeros %d\n', nm{g}, sum(m), ...
            median(lostC(m)), min(lostC(m)), max(lostC(m)), sum(lostC(m) == 0), ...
            median(retOnly(m)), min(retOnly(m)), max(retOnly(m)), sum(retOnly(m) == 0));
    end
    for k = unique(kk)'
        a = grp == 1 & kk == k; b = grp == 2 & kk == k; c = grp == 3 & kk == k;
        p('start %2d (rep %3d): lostC %.3f -> S %.3f / C %.3f; lostR-only %.3f -> S %.3f\n', k, s.reps(k), ...
            lostC(a), lostC(b), lostC(c), retOnly(a), retOnly(b));
    end
    p('long flags: %d networks, wall %.1f min\n', numel(J), O.wallSec / 60);
    save(fullfile(odir, 'long_flags.mat'), '-struct', 'O');
end

function ok = control_fresh(N, J, p)
    d = zeros(numel(N), 2);
    for n = 1:numel(N), d(n, :) = abs([N{n}.uS N{n}.uC] - [J(n).uS J(n).uC]); end
    ok = all(d(:) < 5e-4);
    p('[CONTROL fresh shares] max |uS - stored| %.4g, max |uC - stored| %.4g over %d networks (starts at step 0, ends at step 300) -> %s\n', ...
        max(d(:, 1)), max(d(:, 2)), numel(N), pf(ok));
    bw = arrayfun(@(n) sum(N{n}.lostC ~= ~N{n}.weak), 1:numel(N));
    p('  cooperation-lost flag (simulate_pair) vs member_batch C1 flag: %d disagreements of %d draws\n', sum(bw), 1000 * numel(N));
end


function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end
