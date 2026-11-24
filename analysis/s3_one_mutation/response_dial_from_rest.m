function response_dial_from_rest(nPilot, nWorkers)
% RESPONSE_DIAL_FROM_REST  As a network's response to its partner is scaled by s, does it still
% cooperate with its clone and retaliate when every pair starts at the resting state x* and runs to
% rest? Scores the 267 networks kept from rest (Table S8) against their 160 single-weight mutants at
% each s, and the change L(2) - L(1) among those conditional (cooperating and retaliating) at both
% s = 1 and s = 2; L(s) is the share of a network's mutants whose readout from rest is not positive.
%
% Produces : main text, "Conditional cooperation stays within a mutation's reach of the
%            threshold" (the 191 networks, losses 0.036 to 0.069, printed here); Methods, "Mutation
%            tests in the game"; Fig. 4 (panel C: response_dial_rest.mat and .txt, read by
%            fig4_tipping_input); SI S3.5 (the counts retaliating at each s); Table S8 (last row)
% Inputs   : nPilot (default 0, all 267; n > 0 is a trial run on the first n: checks and timing);
%            nWorkers (default 3), the pool size
%            publication/conditionality_score/response_dial/response_dial_inputs.mat and chunks/
%              (written by response_dial)
% Outputs  : under publication/conditionality_score/response_dial/: response_dial_rest.txt
%            (response_dial_rest_timing_run.txt for a trial run; the lines it prints) and, in the full
%            run once checks 1 to 3 pass, response_dial_rest.mat (S, the scan per network and
%            setting; id, rep, L)
% Run      : response_dial_from_rest   (under 1 min on 3 workers, stored report)
% Seed     : the scan draws nothing; every sign-flip test draws under rng(K.seed + 54000,
%            'twister'), K.seed = 20261028 (conditional_cooperation_measures)

    if nargin < 1, nPilot = 0; end
    if nargin < 2, nWorkers = 3; end
    t0 = tic;
    [K, Hd] = response_dial('handles');
    assert(isequal(K.S, [0 0.5 1 1.5 2]), 'settings differ');
    assert(K.rest.tol == 1e-10 && K.rest.cap == 20000, 'K.rest tol / cap differ from 1e-10 / 20000');
    D = data(K, Hd, nPilot);
    tag = 'response_dial_rest'; if nPilot > 0, tag = [tag '_timing_run']; end
    f = fullfile(K.rd, [tag '.txt']);

    %% Report
    % The scan runs inside: its pool and check lines are part of the report.
    report = @() rest_report(K, Hd, D, nPilot, nWorkers, f, t0); %#ok<NASGU> called inside evalc
    rep = evalc('report()');
    fprintf('%s', rep);

    %% Save
    fid = fopen(f, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function rest_report(K, Hd, D, nPilot, nWorkers, f, t0)
% The scan, its checks and the report's lines; in the full run with the checks passed, also
% response_dial_rest.mat.
    fprintf('===== RESPONSE DIAL REST %s UTC =====\n', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('chunks %s\n', D.dir);
    fprintf('from-rest set %d (%d replicates); from-rest networks run %d%s\n', ...
        D.nSet, numel(unique(D.rep)), numel(D.id), ternary(nPilot > 0, ' (TIMING RUN)', ''));
    fprintf('pair runs from rest, both networks from x*, tol %g, cap %d; membership_test grid %s (delta 0 skipped); settings %s\n', ...
        K.rest.tol, K.rest.cap, mat2str(D.grid, 4), mat2str(D.Z.S));
    nW = Hd.open_workers(nWorkers);
    tr = tic; S = run_all(Hd, D.Z, D.id, D.sp, D.grid, 1:5, K.rest); wall = toc(tr);
    delete(gcp('nocreate'));
    fprintf('pool closed (%d workers); from-rest scan wall %.1f s (%.2f s per network)\n', nW, wall, wall / numel(D.id));
    ok = [check1(D, S), check2(D), check3(S)];
    if nPilot > 0
        timing_run_end(Hd, D, ok, wall); fprintf('\ncompute time %.1f min; report %s\n', toc(t0) / 60, f); return;
    end
    if ~all(ok)
        fprintf('\nresult: not computed (check FAIL: %s)\n', mat2str(ok)); return;
    end
    save(fullfile(K.rd, 'response_dial_rest.mat'), 'S', '-v7.3');
    save(fullfile(K.rd, 'response_dial_rest.mat'), '-struct', 'D', 'id', 'rep', 'L', '-append');
    reading(K, D, S);
    secondary(K, D, S);
    reported(K, D, S);
    fprintf('\nper-network results %s\ncompute time %.1f min; report %s\n', fullfile(K.rd, 'response_dial_rest.mat'), toc(t0) / 60, f);
end


function D = data(K, Hd, nPilot)
% The from-rest set: kept at s = 1, return test (response_dial's check_rest) passed at all five s,
% reaction fifths 2 to 5 (asserted 267). L(n, s): the share of network n's kept mutants (toggles of
% the four dialed connections left out) whose readout from rest is not positive.
    L0 = load(fullfile(K.rd, 'response_dial_inputs.mat'), 'I'); I = L0.I; dirC = fullfile(K.rd, 'chunks');
    Z = Hd.load_set(K, I, (1:500)', dirC); GR = Hd.check_rest(K, Z);
    rest = GR.kept & all(GR.pass(:, 1:5), 2) & I.fifth >= 2; assert(sum(rest) == 267, 'from-rest set is not 267');
    keep = ~(ismember(Z.q, Hd.dialed_slots()) & Z.tog);
    id = find(rest); kp = keep(id, :); lost = ~(Z.rM(id, :, 1:5) > 0);
    assert(all(sum(kp, 2) > 0), 'a network has no kept mutant');
    L = squeeze(sum(lost & kp, 2)) ./ sum(kp, 2);
    [cp, sp] = Hd.census_params(); Zsp = Z.sp; Zsp.sensingMode = stored_value('sensingMode', Zsp.sensingMode);   % stored spelling
    assert(isequal(Zsp, sp), 'stored sp differs from the census sp');
    assert(isequal(size(Z.xstar), [500 4]) && isequal(size(Z.rr0), [500 5]), 'stored xstar / rr0 shape differs');
    nR = numel(id); if nPilot > 0, nR = nPilot; end
    D = struct('Z', Z, 'idAll', id, 'id', id(1:nR), 'rep', I.rep(id), 'L', L, 'nSet', numel(id), ...
        'dir', dirC, 'sp', sp, 'grid', linspace(-3 * cp.sigma_mut, 3 * cp.sigma_mut, 11));
end


function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end


% =========================================================================
% THE SCAN FROM REST (rest_pair: response_dial's runs to rest)
% =========================================================================
function o = scan(dial, rest_pair, M, B, x, s, sp, grid, Rc)
% The network (M, B) dialed to setting s (W = W_s, b = b_s), then the clone pair and the 160
% single-weight mutants (each slot (ii, jj) plus each nonzero grid value, in membership_test's
% order), every pair with b_s on both agents, started at x* and run to rest; the network is agent 1,
% the mutant agent 2. A mutant defects if its C - D is not positive; a defecting
% mutant is met with cooperation if the network cooperates with it (code flag inv). ret
% (retaliates): at least one mutant defects and none is met with cooperation; coop: clone readout
% > 0; cond (conditional): both. worst: the network's
% highest C - D against a defecting mutant (-inf if none); slot and restWorst: that mutant's
% position and rest flag; nNoRest: mutant runs not at rest by the cap.
    [W, b] = dial(M, B, M .* B, x, s);
    [x1, ~, ~, o.restClone] = rest_pair(W, b, x, W, b, x, sp, Rc.tol, Rc.cap); o.rClone = x1(1) - x1(2);
    worst = -inf; slot = NaN; restW = NaN; inv = false; nDef = 0; nNo = 0; k = 0;
    for ii = 1:4
        for jj = 1:4
            for delta = grid
                if delta == 0, continue; end
                Wm = W; Wm(ii, jj) = Wm(ii, jj) + delta;
                [x1, x2, ~, ok] = rest_pair(W, b, x, Wm, b, x, sp, Rc.tol, Rc.cap); k = k + 1; nNo = nNo + ~ok;
                fC = x1(1) > x1(2); nD = ~(x2(1) > x2(2));
                if nD
                    nDef = nDef + 1; r = x1(1) - x1(2);
                    if isfinite(r) && r > worst, worst = r; slot = ii + 4 * (jj - 1); restW = ok; end
                end
                if fC && nD, inv = true; end
            end
        end
    end
    o.nRun = k; o.nDef = nDef; o.nNoRest = nNo; o.inv = inv; o.worst = worst; o.slot = slot;
    o.restWorst = restW; o.ret = nDef > 0 && ~inv; o.coop = o.rClone > 0; o.cond = o.coop && o.ret;
end


function S = run_all(Hd, Z, id, sp, grid, ks, Rc)
% scan() at the settings ks, one parfor job per network; each field stacked as n x numel(ks).
    n = numel(id); C = cell(n, 1); nk = numel(ks); dial = Hd.dial; rest_pair = Hd.rest_pair;
    M = Z.M(:, :, id); B = Z.B(:, :, id); X = Z.xstar(id, :); sv = Z.S(ks);
    parfor j = 1:n
        Mj = M(:, :, j); Bj = B(:, :, j); x = X(j, :)'; r = cell(1, nk);
        for k = 1:nk
            r{k} = scan(dial, rest_pair, Mj, Bj, x, sv(k), sp, grid, Rc);
        end
        C{j} = [r{:}];
    end
    S = struct();
    for fn = fields_()
        S.(fn{1}) = cell2mat(cellfun(@(c) [c.(fn{1})], C, 'uni', 0));               % n x nk
    end
end


function f = fields_()
    f = {'rClone', 'restClone', 'nRun', 'nDef', 'nNoRest', 'inv', 'worst', 'slot', 'restWorst', 'ret', 'coop', 'cond'};
end


% =========================================================================
% CHECKS
% =========================================================================
function ok = check1(D, S)
% Check 1: the clone run from x* reproduces the stored readout rr0 to 1e-12 and is positive, at
% every setting.
    st = D.Z.rr0(D.id, :); d = abs(S.rClone - st); pos = S.rClone > 0;
    ok = all(d(:) <= 1e-12) && all(pos(:));
    fprintf('\nCHECK 1 (clone from rest reproduces stored rr0): %d networks x 5 settings; max |rClone - rr0| per setting %s; > 1e-12 in %d; rClone > 0 in %s of %d; clone runs not at rest by the cap %s -> %s\n', ...
        numel(D.id), mat2str(max(d, [], 1), 3), sum(d(:) > 1e-12), mat2str(sum(pos, 1)), numel(D.id), ...
        mat2str(sum(~S.restClone, 1)), pf(ok));
end


function ok = check2(D)
% Check 2: the mean loss shares print as the fixed string of dial_loss_shares's from-rest row.
    a = mat2str(mean(D.L, 1), 4); want = '[0.01271 0.01904 0.03526 0.05701 0.06579]'; ok = strcmp(a, want);
    fprintf('CHECK 2 (losses reproduce fig4c_dial.txt from rest): mean L over the %d at s = 0 / 0.5 / 1 / 1.5 / 2 %s (want %s) -> %s\n', ...
        size(D.L, 1), a, want, pf(ok));
end


function ok = check3(S)
% Check 3: no network retaliates at s = 0, where the sensor-to-output connections are zero.
    n0 = sum(S.ret(:, 1)); ok = n0 == 0;
    fprintf('CHECK 3 (nothing to answer at s = 0): networks retaliating at s = 0: %d of %d -> %s\n', n0, size(S.ret, 1), pf(ok));
end


% =========================================================================
% TRIAL RUN (nPilot > 0): check 4 and the projected time
% =========================================================================
function timing_run_end(Hd, D, ok, wall)
% Check 4 and the projected full scan (scan time per network x the set, pool start excluded),
% against 2 h.
    ok4 = check4(Hd, D); h = wall / numel(D.id) * D.nSet / 3600;
    fprintf('projection: %d from-rest networks x %.2f s = %.3f h (pool start excluded) -> %s\n', ...
        D.nSet, wall / numel(D.id), h, ternary(h <= 2, 'within 2 h', 'OVER 2 h: STOP'));
    fprintf('checks 1-3 %s; check 4 %s\n', mat2str(ok), pf(ok4));
end


function ok = check4(Hd, D)
% Check 4: rest_pair with no stopping for 600 steps equals simulate_pair bit for bit (both agents
% from x*, constant input b_s, stability threshold 0) for the clone and all 160 mutants of each
% trial-run network at each setting.
    Z = D.Z; sp = D.sp; dial = Hd.dial; rest_pair = Hd.rest_pair; readout = Hd.readout; sp0 = sp; sp0.stabilityThreshold = 0; nP = 0; nBad = 0;
    for j = 1:numel(D.id)
        i = D.id(j); M = Z.M(:, :, i); B = Z.B(:, :, i); x = Z.xstar(i, :)';
        for k = 1:5
            [Ws, b] = dial(M, B, M .* B, x, Z.S(k)); s1 = sp0; s1.initialState1 = x; s1.initialState2 = x;
            for t = 0:160
                Wm = Ws; if t > 0, [ii, jj, dl] = grid_pos(t, D.grid); Wm(ii, jj) = Wm(ii, jj) + dl; end
                f = readout(Ws, Wm, b, b, s1); [a1, a2] = rest_pair(Ws, b, x, Wm, b, x, sp, -1, 600);
                nP = nP + 1; nBad = nBad + ~isequal(f, [a1; a2]);
            end
        end
    end
    ok = nBad == 0;
    fprintf('\nCHECK 4 (timing run, %d from-rest networks x 5 settings x 161 runs): rest_pair with no stopping for 600 steps equals simulate_pair (initialState1/2 = x*, externalSignal1/2 = b_s, stabilityThreshold 0) bit for bit in %d of %d runs -> %s\n', ...
        numel(D.id), nP - nBad, nP, pf(ok));
end


function [ii, jj, dl] = grid_pos(t, grid)
% The t-th mutant in membership_test's order: slot (ii, jj), then the nonzero grid values.
    g = grid(grid ~= 0); m = numel(g); a = floor((t - 1) / m); dl = g(mod(t - 1, m) + 1);
    ii = floor(a / 4) + 1; jj = mod(a, 4) + 1;
end


% =========================================================================
% RESULT, ALSO REPORTED
% =========================================================================
function reading(K, D, S)
% K: networks conditional from rest at s = 1 and s = 2. TOO FEW if n_K < 60 or below half the
% number conditional at s = 1; else, for L(2) - L(1) on K, CONDITIONAL RISE if its interval lies
% above 0 and p < 0.05, NO RISE if the interval holds 0, FALL if it lies below 0.
    n = size(S.cond, 1); sh = @(v) strjoin(arrayfun(@(k) sprintf('s = %g %d (%.4f)', D.Z.S(k), sum(v(:, k)), ...
        sum(v(:, k)) / n), 1:5, 'uni', 0), ', ');
    fprintf('\nC(s) = share of the %d cooperating with their clone from rest: %s\n', n, sh(S.coop));
    fprintf('R(s) = share retaliating from rest: %s\n', sh(S.ret));
    fprintf('Cond(s) = share conditional from rest (C and R): %s\n', sh(S.cond));
    nC = sum(S.cond, 1); Kc = S.cond(:, 3) & S.cond(:, 5); nK = sum(Kc); few = nK < 60 || nK < nC(3) / 2;
    fprintf('K = conditional at s = 1 AND s = 2: n_K = %d (conditional at s = 1: %d; TOO FEW if n_K < 60 or < %.1f)\n', nK, nC(3), nC(3) / 2);
    x = show(K, 'dL = L(2) - L(1) on K', D.L(Kc, 5) - D.L(Kc, 3), D.rep(Kc));
    fprintf('  on K: mean L(1) %.5f, mean L(2) %.5f\n', mean(D.L(Kc, 3)), mean(D.L(Kc, 5)));
    if few, t = 'TOO FEW';
    elseif isempty(x), t = 'not computed (contrast not computable)';
    elseif x.ci(1) > 0 && x.p < 0.05, t = 'CONDITIONAL RISE';
    elseif x.ci(1) <= 0 && x.ci(2) >= 0, t = 'NO RISE';
    elseif x.ci(2) < 0, t = 'FALL';
    else, t = 'NO CATEGORY (interval above 0, p >= 0.05)';
    end
    fprintf('result: %s\n', t);
end


function x = show(K, name, d, rp)
% Cluster-t interval and sign-flip p over replicates, or only n and mean when there are fewer than
% 2 networks or 2 replicates.
    G = numel(unique(rp)); x = [];
    if numel(d) < 2 || G < 2
        fprintf('  %s: n %d (%d replicates), mean %.4f; interval not computable\n', name, numel(d), G, mean(d)); return;
    end
    [x.e, x.ci] = cluster_t(d, rp); x.p = signflip(K, d, rp);
    fprintf('  %s: n %d (%d replicates) %.4f [%.4f, %.4f]; sign-flip p %.4f\n', name, numel(d), G, x.e, x.ci, x.p);
end


function secondary(K, D, S)
% Beside the reading: the rise on two other subsets and on all 267; the last is compared with the
% fixed string of dial_loss_shares's from-rest UPPER line.
    fprintf('\n==================== also reported: the rise on other subsets ====================\n');
    k15 = S.cond(:, 3) & S.cond(:, 4);
    show(K, 'L(1.5) - L(1) on conditional at s = 1 AND 1.5', D.L(k15, 4) - D.L(k15, 3), D.rep(k15));
    k2n1 = S.cond(:, 5) & ~S.cond(:, 3);
    show(K, 'L(2) - L(1) on conditional at s = 2, not at s = 1', D.L(k2n1, 5) - D.L(k2n1, 3), D.rep(k2n1));
    x = show(K, 'L(2) - L(1) on all 267', D.L(:, 5) - D.L(:, 3), D.rep);
    b = sprintf('%.4f [%.4f, %.4f] p %.4f', x.e, x.ci, x.p); want = '0.0305 [0.0253, 0.0358] p 0.0001';
    fprintf('    control against fig4c_dial.txt UPPER: %s (want %s) -> %s\n', b, want, pf(strcmp(b, want)));
end


function reported(K, D, S)
% Counts per setting. Stored rDr <= 0: the network defects back against the default defector from
% rest.
    fprintf('\n==================== also reported: counts per setting ====================\n');
    n = size(S.ret, 1); m = @(v) mat2str(sum(v, 1));
    fprintf('settings s = %s, counts of %d\n', mat2str(D.Z.S), n);
    fprintf('no defecting mutant %s; worst-case readout exactly 0 %s; any mutant invades %s\n', m(S.nDef == 0), m(S.worst == 0), m(S.inv));
    fprintf('not at rest by the cap (%d steps): clone runs %s; mutant runs %s of %d; worst-case mutant runs %s\n', ...
        K.rest.cap, m(~S.restClone), m(S.nNoRest), sum(S.nRun(:, 1)), m(S.restWorst == 0));
    w = isfinite(S.worst);
    fprintf('defecting mutants per network, median: %s; worst-case readout, median: %s\n', ...
        mat2str(median(S.nDef, 1)), mat2str(arrayfun(@(k) median(S.worst(w(:, k), k)), 1:5), 4));
    fprintf('network defects back against the built-in defector from rest (stored rDr <= 0): %s\n', m(D.Z.rDr(D.id, :) <= 0));
    pats = S.cond * [16; 8; 4; 2; 1];
    fprintf('conditional status across settings (s = 0, 0.5, 1, 1.5, 2): %s\n', strjoin(arrayfun(@(p) sprintf('%s: %d', ...
        dec2bin(p, 5), sum(pats == p)), unique(pats)', 'uni', 0), ', '));
end
