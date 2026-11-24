function F = dial_loss_shares(reportFile)
% DIAL_LOSS_SHARES  How does the share of a network's mutations that end cooperation with its clone
% change as its response to its partner is scaled? From the stored outcomes of response_dial (no
% dynamics are run), prints under each of the three ways of counting a loss the mean share at each
% setting s, with its interval, and the contrasts between settings.
%
% Produces : SI S3.5; Table S8 (its first three rows, one per counting rule), printed here
% Inputs   : reportFile (default publication/conditionality_score/response_dial/fig4c_dial.txt)
%            publication/conditionality_score/response_dial/response_dial_inputs.mat and chunks/chunk_<k>.mat
%              (written by response_dial, stages 'checks' and 'full')
% Outputs  : reportFile (the lines this run prints); returns F, one struct per rule
% Run      : dial_loss_shares   (about 30 s); the control line compares the shares and the
%            L(1) - L(0.5) contrasts with those response_dial_recount printed
% Seed     : K.seed = 20261028, Mersenne twister: the interval of rule r at its j-th setting
%            resamples replicates under K.seed + 91000 + 10 r + j; sign-flip tests K.seed + 54000
%
% Rules (stored outcome per mutant and setting): the per-mutation census rule (cMut; s <= 1), the
% Fig. 5 rule (eMut, code tag 'encounter_rule'; s <= 1) and from rest (rM; all five settings), each defined
% in SI S3.5 and in response_dial's dial_focal.

    if nargin < 1
        d = fullfile(grn_data_file('publication'), 'conditionality_score', 'response_dial');
        if ~isfolder(d), mkdir(d); end
        reportFile = fullfile(d, 'fig4c_dial.txt');
    end
    [K, Hd] = response_dial('handles');
    [~, ~, NP] = conditional_cooperation_measures('np_handles');
    L0 = load(fullfile(K.rd, 'response_dial_inputs.mat'), 'I', 'cls'); I = L0.I; cls = L0.cls;
    Z = Hd.load_set(K, I, (1:500)', fullfile(K.rd, 'chunks'));
    keep = ~(ismember(Z.q, Hd.dialed_slots()) & Z.tog);
    R = rules(K, Hd, Z, I);

    %% Report
    report = @() shares_report(K, NP, R, keep, I, cls, reportFile); %#ok<NASGU> called inside evalc
    [rep, F] = evalc('report()');
    fprintf('%s', rep);

    %% Save
    fid = fopen(reportFile, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function F = shares_report(K, NP, R, keep, I, cls, reportFile)
% Each rule's shares and contrasts, the control and the description; prints the report's lines.
    fprintf('===== FIG 4C RESPONSE DIAL %s UTC =====\n', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('chunks %s; toggles on the four dialed connections excluded\n', fullfile(K.rd, 'chunks'));
    F = struct();
    for r = 1:numel(R)
        F.(R(r).tag) = one_rule(K, NP, R(r), keep, I.rep, r);
    end
    control(F);
    describe(F, R, cls, keep);
    fprintf('\nReport %s\n', reportFile);
end


function R = rules(K, Hd, Z, I)
% The three rules: a mutant is lost if its census readout is not positive (s <= 1), if its own-copy
% cooperation is below 0.5 (s <= 1), or if its readout from rest is not positive (all five s).
% Sets, all within reaction fifths 2 to 5: the developmental rules keep networks kept at s = 1
% (census readout > 0, or own-copy cooperation >= 0.99) that pass the development test at s = 0,
% 0.5 and 1; from rest keeps networks kept at s = 1 that pass the return test at all five s (the
% two tests: check_dev and check_rest in response_dial).
    GD = Hd.check_dev(K, Z, I.self); GR = Hd.check_rest(K, Z);
    devSet = @(kept) kept & all(GD.pass(:, K.dev), 2) & I.fifth >= 2;
    rest = GR.kept & all(GR.pass(:, 1:5), 2) & I.fifth >= 2;
    R = struct('tag', {'census', 'encounter_rule', 'rest'}, ...
        'name', {'census (fixed start, 600 steps, C - D > 0)', 'encounter rule (100 encounters, own copy < 0.5)', 'from rest (C - D > 0 at rest)'}, ...
        'lost', {~(Z.cMut > 0), Z.eMut < 0.5, ~(Z.rM > 0)}, ...
        'set', {devSet(GD.keptC), devSet(GD.keptE), rest}, 'k', {K.dev, K.dev, 1:5});
    assert(sum(R(1).set) == 199 && sum(R(2).set) == 199 && sum(R(3).set) == 267, 'set sizes differ from the recount');
    for r = 1:3, R(r).s = K.S(R(r).k); end
end


function O = one_rule(K, NP, R, keep, rep, r)
% L(n, s): the share of network n's kept mutants lost at setting s (toggles of the four dialed
% connections left out, since their size scales with s). Per setting, the mean over networks with
% a 95% percentile interval over resampled replicates; LOWER L(1) - L(0.5), STEP L(1) - L(0) and,
% from rest, UPPER L(2) - L(1), each with a cluster-t interval and a sign-flip p over replicates.
    p = R.set; kp = keep(p, :); lost = R.lost(p, :, R.k);
    L = squeeze(sum(lost & kp, 2)) ./ sum(kp, 2); rp = rep(p);
    assert(all(sum(kp, 2) > 0), 'a network has no kept mutant');
    O = struct('s', R.s, 'L', mean(L, 1), 'ci', nan(2, numel(R.k)), 'n', sum(p), 'nrep', numel(unique(rp)), ...
        'nmut', sum(kp(:)), 'Lnet', L, 'rep', rp, 'idx', find(p));
    for j = 1:numel(R.k)
        [~, O.ci(:, j)] = NP.boot(K, rp, @(i) mean(L(i, j)), K.seed + 91000 + 10 * r + j);
    end
    j1 = find(R.s == 1); jh = find(R.s == 0.5); j0 = find(R.s == 0);
    O.lower = contrast(K, L(:, j1) - L(:, jh), rp);
    O.step = contrast(K, L(:, j1) - L(:, j0), rp);
    if any(R.s == 2), O.upper = contrast(K, L(:, R.s == 2) - L(:, j1), rp); else, O.upper = []; end
    fprintf('\n%s: %d networks (%d replicates), %d kept mutants\n', R.name, O.n, O.nrep, O.nmut);
    for j = 1:numel(R.k)
        fprintf('  s = %.1f  L %.5f [%.5f, %.5f] (bootstrap over replicates)\n', R.s(j), O.L(j), O.ci(1, j), O.ci(2, j));
    end
    pr = @(t, x) fprintf('  %-22s %.4f [%.4f, %.4f]; sign-flip p %.4f\n', t, x.e, x.ci, x.p);
    pr('LOWER L(1) - L(0.5)', O.lower); pr('STEP  L(1) - L(0)', O.step);
    if ~isempty(O.upper), pr('UPPER L(2) - L(1)', O.upper); end
end


function x = contrast(K, d, rp)
    [x.e, x.ci] = cluster_t(d, rp); x.p = signflip(K, d, rp);
end


function control(F)
% Compares each rule's L by setting and LOWER contrast with the strings below, the values
% response_dial_recount prints for the same sets (stored response_dial_recount.txt); reports only.
    want = {'census', '[0.01225 0.01734 0.08144]', '0.0641 [0.0529, 0.0753] p 0.0001'; ...
            'encounter_rule', '[0.02276 0.02553 0.06339]', '0.0379 [0.0270, 0.0487] p 0.0001'; ...
            'rest',   '[0.01271 0.01904 0.03526 0.05701 0.06579]', '0.0162 [0.0115, 0.0209] p 0.0001'};
    ok = true;
    fprintf('\n[CONTROL] against response_dial_recount.txt (excluded, primary sets)\n');
    for r = 1:3
        O = F.(want{r, 1}); a = mat2str(O.L, 4);
        b = sprintf('%.4f [%.4f, %.4f] p %.4f', O.lower.e, O.lower.ci, O.lower.p);
        okr = strcmp(a, want{r, 2}) && strcmp(b, want{r, 3}); ok = ok && okr;
        fprintf('  %-14s L %s (want %s); LOWER %s (want %s) -> %s\n', want{r, 1}, a, want{r, 2}, b, want{r, 3}, pf(okr));
    end
    if ~ok, fprintf('CONTROL FAILED: stored outcomes do not reproduce the recount\n'); end
end


function describe(F, R, cls, keep)
% The census set's membership classes (1 member, 2 cooperates with its clone but does not
% retaliate, 3 other), its overlap with the other sets, and its kept and left-out mutants.
    c = F.census; nm = {'member', 'cooperate_only', 'other'};
    fprintf('\n[DESCRIBE] census set: classes %s\n', strjoin(arrayfun(@(k) sprintf('%s %d', nm{k}, sum(cls(c.idx) == k)), 1:3, 'uni', 0), ', '));
    fprintf('  census set within the from-rest set: %d of %d (from-rest set %d)\n', numel(intersect(c.idx, F.rest.idx)), c.n, F.rest.n);
    fprintf('  census set equals the encounter-rule set: %s\n', pf(isequal(c.idx, F.encounter_rule.idx)));
    fprintf('  kept mutants per network (census set): median %d, range [%d, %d]; excluded %d of %d\n', ...
        median(sum(keep(c.idx, :), 2)), min(sum(keep(c.idx, :), 2)), max(sum(keep(c.idx, :), 2)), sum(~keep(c.idx, :), 'all'), 30 * c.n);
    fprintf('  census L ratio L(1) / L(0.5) %.2f, L(1) / L(0) %.2f\n', c.L(3) / c.L(2), c.L(3) / c.L(1));
    fprintf('  census, networks losing none of their kept mutants at s = 0 / 0.5 / 1: %s of %d\n', mat2str(sum(c.Lnet == 0, 1)), c.n);
    for r = 1:numel(R), fprintf('  rule %d: %s\n', r, R(r).name); end
end


function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end
