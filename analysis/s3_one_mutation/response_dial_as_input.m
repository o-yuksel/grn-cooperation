function response_dial_as_input(nWorkers)
% RESPONSE_DIAL_AS_INPUT  Does a mutation act on a network along the response dial as the constant
% input it adds at rest? For 100 of the 267 networks kept from rest, each kept mutation is replaced by
% that input, du = dW_s x* (dW_s, the change the mutation makes to the scaled wiring), added to
% both copies with the wiring unchanged; the pair runs from x* to rest at s = 0.5 and 2, and the
% outcome is compared with the stored outcome of the mutant itself.
%
% Produces : SI S3.5, printed here: the rise in losses from s = 0.5 to 2 with each mutation replaced
%            by its input, against the mutants' own rise; the share of mutant losses the inputs catch
% Inputs   : nWorkers (default 3), the pool size
%            publication/conditionality_score/response_dial/response_dial_inputs.mat and chunks/
%              (written by response_dial)
% Outputs  : publication/conditionality_score/response_dial/response_dial_signal.txt
% Run      : response_dial_as_input   (under 1 min on 3 workers, stored report); the first 50 runs
%            time the rest and the projection for the remaining ones is printed
% Seed     : K.seed = 20261028, Mersenne twister: the 100 networks are sort(randperm(267, 100))
%            under K.seed + 90000; every sign-flip test K.seed + 54000; the runs draw nothing
%
% Per network, L_mut and L_sig: the share of its kept mutations whose readout at rest is not
% positive, for the mutant (stored rM) and for its input du. SLOPE = mean over networks of
% L(2) - L(0.5); R = SLOPE_sig / SLOPE_mut. catch: the share of mutant losses du also loses.

    if nargin < 1, nWorkers = 3; end
    t0 = tic;
    [K, Hd] = response_dial('handles');
    L = load(fullfile(K.rd, 'response_dial_inputs.mat'), 'I');
    D = Hd.rest_data(K, L.I, fullfile(K.rd, 'chunks')); E = Hd.exclusion(D);   % the from-rest set (338), primary 267
    assert(K.rest.tol == 1e-10 && K.rest.cap == 20000, 'rest tolerance / cap differ from 1e-10 / 20000');
    f = fullfile(K.rd, 'response_dial_signal.txt');

    %% Report
    % The runs go inside: their pool, timing and check lines are part of the report.
    report = @() signal_report(K, Hd, D, E, nWorkers, f, t0); %#ok<NASGU> called inside evalc
    rep = evalc('report()');
    fprintf('%s', rep);

    %% Save
    fid = fopen(f, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function signal_report(K, Hd, D, E, nWorkers, f, t0)
% The sample, the signal runs, the checks, the measures and the descriptive tables.
    fprintf('===== RESPONSE DIAL SIGNAL %s UTC =====\n', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('chunks %s; set %d, primary %d; tol %g, cap %d\n', ...
        D.dir, sum(D.set), sum(D.prim), K.rest.tol, K.rest.cap);
    S = sample(K, Hd, D, E);
    nW = Hd.open_workers(nWorkers); cp = onCleanup(@() delete(gcp('nocreate')));
    S = signal_runs(K, Hd, D, S, t0);
    ok = [check1(K, Hd, D, S), check2(K, Hd, D, S)];
    delete(gcp('nocreate'));
    fprintf('pool closed (%d workers used)\n', nW);
    A = measures(K, S);
    cs = agreement(S);
    reading(ok, A, cs);
    fprintf('\n==================== also reported ====================\n');
    fprintf('SLOPE_mut on the sample %.4f\n', A.slm);
    fprintf('not at rest by the cap: signal runs at s = 0.5 / 2: %d / %d of %d each; stored mutant runs not at rest among the same mutants: %d / %d\n', ...
        sum(~S.okS(:, 1)), sum(~S.okS(:, 2)), S.n, sum(~S.okM(:, 1)), sum(~S.okM(:, 2)));
    split(Hd, S);
    du_sizes(S);
    fprintf('\ncompute time %.1f min; report %s\n', toc(t0) / 60, f);
end


% =========================================================================
% SAMPLE AND JOBS
% =========================================================================
function S = sample(K, Hd, D, E)
% The 100 sampled networks and their kept mutants, one row per (network, mutant), at s = 0.5
% (k = 2) and s = 2 (k = 5): du = (dW_s) x*, the change the mutation makes to W_s times x*.
    Z = D.Z; dial = Hd.dial; apply_mut = Hd.apply_mut; pid = find(D.prim); rs = rng; rng(K.seed + 90000, 'twister'); ix = sort(randperm(numel(pid), 100)); rng(rs);
    id = pid(ix); [t, j] = find(E.keep(id, :)'); net = id(j); ks = [2 5]; n = numel(net);
    S = struct('id', id, 'net', net, 't', t, 'n', n, 'ks', ks, 'du', zeros(n, 4, 2), 'Ws', {cell(n, 2)}, 'b', zeros(n, 4, 2));
    for r = 1:n
        i = net(r); M = Z.M(:, :, i); B = Z.B(:, :, i); W = M .* B; x = Z.xstar(i, :)';
        for a = 1:2
            [Ws, bias, Bs] = dial(M, B, W, x, Z.S(ks(a)));
            [Mm, Bm] = apply_mut(M, Bs, Z.q(i, t(r)), Z.tog(i, t(r)), Z.d(i, t(r)));
            S.du(r, :, a) = (Mm .* Bm - M .* Bs) * x; S.Ws{r, a} = Ws; S.b(r, :, a) = bias;
        end
    end
    li = sub2ind(size(Z.q), net, t);
    S.q = Z.q(li); S.tog = Z.tog(li); S.pres = E.pres(li); S.rep = D.rep(net);
    for a = 1:2
        S.lostM(:, a) = D.lost(sub2ind(size(D.lost), net, t, ks(a) * ones(n, 1)));
        S.okM(:, a) = Z.rMOk(sub2ind(size(Z.rMOk), net, t, ks(a) * ones(n, 1)));
    end
    fprintf('sample: %d networks (%d replicates), %d kept mutants; settings s = %s; signal runs %d\n', ...
        numel(id), numel(unique(D.rep(id))), n, mat2str(Z.S(ks)), 2 * n);
end


function S = signal_runs(K, Hd, D, S, t0)
% Each du run: W_s on both copies, constant input b_s + du on both, from x* to rest; lost if the
% readout there is not positive. The first 50 runs time the rest; the projection for the remaining
% runs is printed.
    sp = D.Z.sp; x = num2cell(D.Z.xstar(S.net, :)', 1)'; J = struct('W', {}, 'b', {}, 'x', {});
    for a = 1:2
        for r = 1:S.n, J(end + 1, 1) = struct('W', S.Ws{r, a}, 'b', S.b(r, :, a)' + S.du(r, :, a)', 'x', x{r}); end %#ok<AGROW>
    end
    tt = tic; R1 = run_jobs(Hd, J(1:50), sp, K.rest); t50 = toc(tt);
    proj = t50 / 50 * (numel(J) - 50); used = toc(t0);
    fprintf('timing: 50 signal runs %.1f s on the pool -> projection %.1f min for the remaining %d; used so far %.1f min\n', ...
        t50, proj / 60, numel(J) - 50, used / 60);
    R = [R1; run_jobs(Hd, J(51:end), sp, K.rest)];
    S.rS = reshape([R.r]', S.n, 2); S.okS = reshape([R.ok]', S.n, 2); S.lostS = ~(S.rS > 0);
end


function R = run_jobs(Hd, J, sp, Rc)
% Each job run to rest with response_dial's rest_pair, one parfor job per run.
    R = repmat(struct('r', NaN, 'ok', false), numel(J), 1); tol = Rc.tol; cap = Rc.cap; rest_pair = Hd.rest_pair;
    parfor j = 1:numel(J)
        R(j) = run_one(rest_pair, J(j), sp, tol, cap);
    end
end


function r = run_one(rest_pair, j, sp, tol, cap)
% One job, both copies from x* to rest: the readout there and whether the run came to rest.
    [x1, ~, ~, ok] = rest_pair(j.W, j.b, j.x, j.W, j.b, j.x, sp, tol, cap); r = struct('r', x1(1) - x1(2), 'ok', ok);
end


% =========================================================================
% CHECKS
% =========================================================================
function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end


function ok = check1(K, Hd, D, S)
% Check 1: du = 0 reproduces the stored parent readout rr0 to 1e-12, both settings.
    Z = D.Z; J = struct('W', {}, 'b', {}, 'x', {}); want = []; dial = Hd.dial;
    for i = S.id'
        M = Z.M(:, :, i); B = Z.B(:, :, i); x = Z.xstar(i, :)';
        for k = S.ks, [Ws, b] = dial(M, B, M .* B, x, Z.S(k)); J(end + 1, 1) = struct('W', Ws, 'b', b, 'x', x); want(end + 1, 1) = Z.rr0(i, k); end %#ok<AGROW>
    end
    R = run_jobs(Hd, J, Z.sp, K.rest); e = max(abs([R.r]' - want)); ok = e <= 1e-12;
    fprintf('\nCHECK 1 (du = 0 reproduces the stored parent readout rr0, %d networks x 2 settings): max |difference| %.3g (<= 1e-12) -> %s\n', ...
        numel(S.id), e, pf(ok));
end


function ok = check2(K, Hd, D, S)
% Check 2: the first 200 sampled kept mutants, rerun as mutants at both settings, reproduce the
% stored rM to 1e-12.
    Z = D.Z; J = struct('W', {}, 'b', {}, 'x', {}); want = []; dial = Hd.dial; apply_mut = Hd.apply_mut;
    for r = 1:min(200, S.n)
        i = S.net(r); M = Z.M(:, :, i); B = Z.B(:, :, i); x = Z.xstar(i, :)';
        for k = S.ks
            [~, b, Bs] = dial(M, B, M .* B, x, Z.S(k)); [Mm, Bm] = apply_mut(M, Bs, Z.q(i, S.t(r)), Z.tog(i, S.t(r)), Z.d(i, S.t(r)));
            J(end + 1, 1) = struct('W', Mm .* Bm, 'b', b, 'x', x); want(end + 1, 1) = Z.rM(i, S.t(r), k); %#ok<AGROW>
        end
    end
    R = run_jobs(Hd, J, Z.sp, K.rest); e = max(abs([R.r]' - want)); ok = e <= 1e-12;
    fprintf('CHECK 2 (the first %d sampled kept mutants rerun as mutants at both settings, %d runs): max |difference| to stored rM %.3g (<= 1e-12) -> %s\n', ...
        min(200, S.n), numel(J), e, pf(ok));
end


% =========================================================================
% MEASURES, AGREEMENT, RESULT
% =========================================================================
function A = measures(K, S)
% Per network L at s = 0.5 and 2 for mutants and for du; SLOPE and its sign-flip p over replicates.
    [u, ~, g] = unique(S.net); rp = accumarray(g, S.rep, [], @(v) v(1));
    Lm = [accumarray(g, S.lostM(:, 1), [], @mean), accumarray(g, S.lostM(:, 2), [], @mean)];
    Ls = [accumarray(g, S.lostS(:, 1), [], @mean), accumarray(g, S.lostS(:, 2), [], @mean)];
    dm = Lm(:, 2) - Lm(:, 1); ds = Ls(:, 2) - Ls(:, 1);
    A = struct('slm', mean(dm), 'sls', mean(ds), 'pm', signflip(K, dm, rp), 'ps', signflip(K, ds, rp));
    A.R = A.sls / A.slm;
    fprintf('\nMEASURES (sample, %d networks; s = 0.5 / 2):\n', numel(u));
    fprintf('  L_mut by setting %s; L_sig by setting %s\n', mat2str(mean(Lm, 1), 4), mat2str(mean(Ls, 1), 4));
    fprintf('  SLOPE_mut %.4f (sign-flip p %.4f); SLOPE_sig %.4f (sign-flip p %.4f); R = SLOPE_sig / SLOPE_mut %.4f\n', A.slm, A.pm, A.sls, A.ps, A.R);
end


function c = agreement(S)
    sl = {'0.5', '2'}; c = zeros(1, 2);
    fprintf('\nAGREEMENT (all kept mutants of the sample)\n');
    for a = 1:2
        m = S.lostM(:, a); s = S.lostS(:, a); T = [sum(m & s), sum(m & ~s); sum(~m & s), sum(~m & ~s)];
        c(a) = T(1, 1) / sum(m); pr = T(1, 1) / sum(s);
        fprintf('  s = %s: mutant lost & signal lost %d | mutant lost & signal kept %d | mutant kept & signal lost %d | mutant kept & signal kept %d; catch %.4f; precision %.4f\n', ...
            sl{a}, T(1, 1), T(1, 2), T(2, 1), T(2, 2), c(a), pr);
    end
end


function reading(ok, A, c)
% SIGNAL AT FULL SIZE if R >= 0.75 and the lower catch >= 0.75; NOT A SIGNAL AT FULL SIZE if
% R < 0.25 or that catch < 0.50; else PARTLY. No reading if check 1 or 2 fails.
    cm = min(c);
    fprintf('\nInputs: R %.4f; c = min catch %.4f; checks 1 %s, 2 %s\n', A.R, cm, pf(ok(1)), pf(ok(2)));
    if ~all(ok), r = sprintf('not computed (check %s FAIL)', strjoin(arrayfun(@num2str, find(~ok), 'uni', 0), ', '));
    elseif A.R >= 0.75 && cm >= 0.75, r = 'SIGNAL AT FULL SIZE';
    elseif A.R < 0.25 || cm < 0.50, r = 'NOT A SIGNAL AT FULL SIZE';
    else, r = 'PARTLY'; end
    fprintf('result: %s\n', r);
end


% =========================================================================
% ALSO REPORTED
% =========================================================================
function split(Hd, S)
% Agreement by slot group and by type (toggle off: a toggle of a present edge; toggle on: of an
% absent edge; weight change).
    G = {'inputs to sensors (rows 3-4)', ismember(S.q, Hd.sensor_rows()); 'the four dialed connections', ismember(S.q, Hd.dialed_slots()); ...
        'inputs among C and D', ismember(S.q, Hd.cd_slots()); 'toggle off', S.tog & S.pres; 'toggle on', S.tog & ~S.pres; 'weight change', ~S.tog};
    sl = {'0.5', '2'};
    fprintf('agreement split:\n  %-30s | s   | n    | mut lost | sig lost | both | catch  | precision\n', 'group');
    for g = 1:size(G, 1)
        for a = 1:2
            w = G{g, 2}; m = S.lostM(w, a); s = S.lostS(w, a); b = sum(m & s);
            fprintf('  %-30s | %-3s | %4d | %8d | %8d | %4d | %6.3f | %9.3f\n', G{g, 1}, sl{a}, sum(w), sum(m), sum(s), b, b / sum(m), b / sum(s));
        end
    end
end


function du_sizes(S)
% Median Euclidean norm of du among mutant losses that du misses and among those it catches.
    sl = {'0.5', '2'};
    for a = 1:2
        n = vecnorm(S.du(:, :, a), 2, 2); m = S.lostM(:, a); s = S.lostS(:, a);
        fprintf('median |du| among mutant losses, s = %s: signal misses %.4g (n %d); signal catches %.4g (n %d)\n', ...
            sl{a}, median(n(m & ~s)), sum(m & ~s), median(n(m & s)), sum(m & s));
    end
end
