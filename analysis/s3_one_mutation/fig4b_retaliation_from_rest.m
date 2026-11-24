function T = fig4b_retaliation_from_rest(pubDir)
% FIG4B_RETALIATION_FROM_REST  Of the 30 + 30 final wirings of the two safety searches, how many
% cooperate and retaliate when scored from rest (Fig. 4C's rule, response unscaled: s = 1)?
%
% Produces : Fig. 4B, bottom (the share that retaliate; fig4_tipping_input reads the report);
%            SI S3.4, "Only the final wirings that had to stay members retaliate from rest"
% Inputs   : pubDir   the publication folder (default grn_data_file('publication'))
%            publication/safety_search/long_flags.mat   (which wirings; safety_search_losses)
%            publication/safety_search/search_main_starts_1_4.mat, search_main_starts_5_12.mat and
%              long_jobs/long_<start>_<S|C>.mat   (the search records; run_safety_search)
%            publication/conditionality_score/response_dial/response_dial_inputs.mat and chunks/ (response_dial),
%              response_dial_rest.mat (response_dial_from_rest): the opening check only
% Outputs  : publication/safety_search/fig4b_retaliation_rest.txt (the lines this run prints);
%            returns T, the per-wiring counts of each search
% Run      : fig4b_retaliation_from_rest   (under 1 min)
% Seed     : none (no random draws)
%
% Rest x*: the clone pair's state after 600 steps from the fixed start, iterated until no gene
% changes by 1e-10 per step. From x*, a wiring cooperates if its clone pair, run to rest, has
% x_C > x_D; it retaliates if at least one of its 160 single-weight mutants defects and each one
% that defects meets defection. Pair runs not settled by 20,000 steps are read at that step and
% counted. The opening check must reproduce 5 of the 267 networks of response_dial_rest.mat.
    if nargin < 1, pubDir = grn_data_file('publication'); end
    d = fullfile(pubDir, 'safety_search'); f = fullfile(d, 'fig4b_retaliation_rest.txt');

    %% Report
    report = @() counts_report(pubDir, d, f); %#ok<NASGU> called inside evalc
    [rep, T] = evalc('report()');
    fprintf('%s', rep);

    %% Save
    fid = fopen(f, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function T = counts_report(pubDir, d, f)
% The check, then the counts of each group; prints the report's lines.
    fprintf('===== FIG 4B RETALIATION FROM REST %s UTC =====\n', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    sp = census_sp(); grid = linspace(-3 * 0.2, 3 * 0.2, 11); Rc = struct('tol', 1e-10, 'cap', 20000); t0 = tic;
    fprintf('rule: panel C from-rest rule at s = 1, no dial; rest_pair tol %g cap %d; grid %s (delta 0 skipped)\n', Rc.tol, Rc.cap, mat2str(grid, 4));
    T = struct();
    if ~check_dial5(pubDir, sp, grid, Rc), fprintf('\nNot computed (CHECK failed)\n'); return; end
    G = groups(d);
    for g = 1:2, T.(G(g).tag) = count(G(g), sp, grid, Rc); end
    fprintf('\nFROM REST, s = 1, no dial\n');
    for g = 1:2
        S = T.(G(g).tag); n = numel(S.ret);
        fprintf('  %s: n %d; cooperate %d; retaliate %d; conditional %d; refine not converged %d; clone not at rest %d; single-weight mutant runs not at rest %d of %d\n', ...
            G(g).name, n, sum(S.coop), sum(S.ret), sum(S.cond), sum(~S.refOk), sum(~S.restClone), sum(S.nNoRest), sum(S.nRun));
        fprintf('    conditional share %.4f (%d of %d); cooperate but not retaliate %d; retaliate but not cooperate %d; no defecting mutant %d; invaded %d\n', ...
            mean(S.cond), sum(S.cond), n, sum(S.coop & ~S.ret), sum(~S.coop & S.ret), sum(S.nDef == 0), sum(S.inv));
    end
    fprintf('\ncompute time %.1f min; report %s\n', toc(t0) / 60, f);
end


function S = count(g, sp, grid, Rc)
% Per wiring: x* as response_dial finds it, then scan() from rest at s = 1.
    n = size(g.W, 3); z4 = zeros(4, 1); C = cell(n, 1); refOk = false(n, 1);
    for i = 1:n
        W = g.W(:, :, i); f0 = readout(W, W, z4, z4, sp); [xstar, ~, refOk(i)] = refine(W, z4, f0(1:4), sp);
        C{i} = scan(W, z4, xstar, sp, grid, Rc);
    end
    S = struct('refOk', refOk);
    for fn = {'rClone', 'restClone', 'nRun', 'nDef', 'nNoRest', 'inv', 'ret', 'coop', 'cond'}
        S.(fn{1}) = cellfun(@(o) double(o.(fn{1})), C);
    end
    S.coop = S.coop == 1; S.ret = S.ret == 1; S.cond = S.cond == 1; S.restClone = S.restClone == 1; S.inv = S.inv == 1;
end


function ok = check_dial5(pubDir, sp, grid, Rc)
% Opening check: this file's scan at s = 1 on 5 of the 267 from-rest dial networks must equal
% response_dial_rest.mat's S.ret and S.rClone bit for bit; the counts are not printed otherwise.
    rd = fullfile(pubDir, 'conditionality_score', 'response_dial');
    [K, Hd] = response_dial('handles'); K.rd = rd;
    L0 = load(fullfile(rd, 'response_dial_inputs.mat'), 'I'); Z = Hd.load_set(K, L0.I, (1:500)', fullfile(rd, 'chunks'));
    R = load(fullfile(rd, 'response_dial_rest.mat')); ix = [1 60 120 200 267]; eq = true(1, 5); sx = zeros(1, 5);
    for k = 1:5
        id = R.id(ix(k)); M = Z.M(:, :, id); B = Z.B(:, :, id); W = M .* B;
        [Ws, bias] = dial(M, B, W, Z.xstar(id, :)', 1);
        o = scan(Ws, bias, Z.xstar(id, :)', sp, grid, Rc);
        eq(k) = isequal(o.ret, logical(R.S.ret(ix(k), 3))) && isequal(o.rClone, R.S.rClone(ix(k), 3));
        f0 = readout(W, W, zeros(4, 1), zeros(4, 1), sp); xr = refine(W, zeros(4, 1), f0(1:4), sp); sx(k) = isequal(xr, Z.xstar(id, :)');
    end
    ok = all(eq);
    fprintf('CHECK (scan at s = 1 on networks %s of the 267 = S.ret and S.rClone of response_dial_rest.mat, bit for bit): %s; refine reproduces stored x* bit for bit: %s -> %s\n', ...
        mat2str(ix), mat2str(eq), mat2str(logical(sx)), pf(ok));
end


function sp = census_sp()
% simulate_pair parameters of the membership test (saturating contrast; 600 development steps),
% as census_params in response_dial_from_rest.m.
    sp = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), 'externalSignal2', zeros(4, 1));
end


function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end


% =========================================================================
% The two groups of final wirings, as long_flags.mat lists them
% =========================================================================
function G = groups(d)
% Panel B's wirings (long_flags.mat J; grp 2 = must stay conditional, grp 3 = need only cooperate,
% 30 each) with their stored search records.
    L = load(fullfile(d, 'long_flags.mat'), 'J'); J = L.J;
    V = [getfield(load(fullfile(d, 'search_main_starts_1_4.mat'), 'R'), 'R'); ...
         getfield(load(fullfile(d, 'search_main_starts_5_12.mat'), 'R'), 'R')]; %#ok<GFLD>
    spec = {2, 'S', 'mustStay', 'must stay conditional'; 3, 'C', 'needCoop', 'need only cooperate'};
    for g = 1:2
        j = J([J.grp] == spec{g, 1}); assert(numel(j) == 30, 'group %d: %d networks, not 30', spec{g, 1}, numel(j));
        n = numel(j); mem = nan(n, 1); nDef = nan(n, 1); src = cell(n, 1); W = zeros(4, 4, n);
        for i = 1:n
            [r, src{i}] = record(d, V, j(i).k, spec{g, 2});
            assert(isequal(r.M1, j(i).M) && isequal(r.B1, j(i).B), 'group %d start %d: J wiring differs from the record', spec{g, 1}, j(i).k);
            assert(r.ckSteps(end) == size(r.trace, 1) && j(i).step == size(r.trace, 1), 'last checkpoint is not the final step');
            mem(i) = r.fresh.member(end); nDef(i) = r.fresh.nDefMut(end); W(:, :, i) = double(j(i).M) .* double(j(i).B);
        end
        G(g) = struct('grp', spec{g, 1}, 'mode', spec{g, 2}, 'tag', spec{g, 3}, 'name', spec{g, 4}, 'k', [j.k]', ...
            'W', W, 'mem', mem, 'nDef', nDef, 'src', {src}); %#ok<AGROW>
    end
end


function [r, src] = record(d, V, k, mode)
    f = fullfile(d, 'long_jobs', sprintf('long_%02d_%s.mat', k, mode));
    if isfile(f), r = getfield(load(f), 'r'); src = f; return; end %#ok<GFLD>
    assert(mode == 'C', 'no long record for start %d mode S', k);
    i = find(cellfun(@(x) x.k == k && x.mode == 'C', V)); assert(numel(i) == 1, 'no stored C record for start %d', k);
    r = V{i}; src = 'search_main';
end



% =========================================================================
% The from-rest rule, as in response_dial.m and response_dial_from_rest.m
% =========================================================================
function [x, it, ok] = refine(W, bias, x, sp)
% The clone pair's Euler map with the partner equal to itself, until the largest step change is
% below 1e-10 (at most 1e6 steps; ok = false if it does not settle).
    km = sp.michaelisConstant; mu = sp.degradationRate; dt = sp.timeStep; ok = false;
    for it = 1:1e6
        u = max(W * x + partner_signal(x, sp) + bias, 0);
        xn = max(x + dt * (u ./ (km + u) - mu * x), 0);
        ch = max(abs(xn - x)); x = xn;
        if ch < 1e-10, ok = true; return; end
    end
end


function s = partner_signal(x, sp)
% simulate_pair's sensing (Hill, contrast) with the partner equal to x.
    n = sp.hillCoefficient; cn = x(1)^n; dn = x(2)^n; den = cn + dn + 1e-12;
    g = sp.sensorIntegrationRate * (cn / den - dn / den); s = [0; 0; g; -g];
end


function o = scan(W, b, x, sp, grid, Rc)
% The clone and the 160 single-weight mutants, every pair from x*, run to rest; focal = agent 1,
% mutant = agent 2. Retaliates: at least one mutant defects and none that defects meets cooperation.
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
