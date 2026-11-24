function varargout = response_dial_networks(mode, K, H, varargin)
% RESPONSE_DIAL_NETWORKS  Draws the 500 saturating-contrast networks of the response dial (SI S3.5)
% and 30 random mutations of each, and measures each network's reaction to 80 partners from the
% same populations: its probability of cooperating when the partner cooperated minus that when
% the partner defected.
%
% Produces : as input only, for main text, "Conditional cooperation stays within a mutation's
%            reach of the threshold"; Methods, "Mutation tests in the game"; Fig. 4; SI S3.4;
%            SI S3.5; Table S8 (response_dial reads pool.mat and chunks/)
% Inputs   : mode (the phase), K and H (constants and helpers of conditional_cooperation_measures);
%            for 'np_timing_run' and 'np_score' optionally nWorkers (default 3), the pool size
%            publication/<architecture>_100reps.mat   (the run configs, and the stored runs of
%              saturating contrast; run_signaling_architectures)
%            'np_check2' also: publication/rho_by_architecture.mat (per_mutation_census) and
%              retaliation_fixed_opponent/retaliation_fixed_opponent_census.mat
%              (measure_retaliation_fixed_opponent)
% Outputs  : under publication/conditionality_score/natural_partners/: pool.mat (partners Pt,
%            networks F with their mutants, K), checks.mat (pass flags), timing_run_manifest.mat and
%            timing_run/ (the trial on networks 1 to 50), chunks/chunk_<k>.mat (25 networks each),
%            natural_partners[_timing_run].txt and log_<phase>.txt; 'np_handles' returns K and helpers
% Run      : through conditional_cooperation_measures, whose header says what each phase does
% Seed     : K.seed = 20261028, Mersenne twister: partners K.seed + 9101, networks K.seed + 9102,
%            mutations of network f K.seed + 10000 + f; checks and bootstraps further fixed offsets
%            (9401 to 9614; calibration 30000 on); encounters under mutate_with_refusal's rng(20261022)
%
% The pool is every individual at generations 100, 200, 300, 376 and 500 of the 100 replicates
% (500,000). Networks are drawn by individual, so a common genotype is more likely to be drawn,
% and each genotype at most once; partners are 80 individuals drawn without replacement. What
% 'np_analyze' prints (checks 3 to 5, measures 1 to 3 defined in the function readings, and the
% descriptive tables) feeds no paper item.

    K.np = fullfile(K.out, 'natural_partners'); if ~isfolder(K.np), mkdir(K.np); end
    K.npChunk = 25; K.nFocal = 500; K.nPart = 80; K.nMut = 30; K.nPilot = 50; K.encSeed = 20261022;
    K.gens = [100 200 300 376 500];
    switch mode
        case 'np_pool',    logged(K, 'pool', @() do_pool(K, H));
        case 'np_check1',  logged(K, 'check1', @() do_check1(K, H));
        case 'np_check2',  logged(K, 'check2', @() do_check2(K));
        case 'np_timing_run',   logged(K, 'timing_run', @() do_timing_run(K, H, opt(varargin, 1, 3)));
        case 'np_score',   logged(K, 'score', @() do_score(K, H, opt(varargin, 1, 3)));
        case 'np_analyze', do_analyze(K, H, opt(varargin, 1, 'full'));
        case 'np_handles'    % read-only handles for the response dial (no evaluation, no output file)
            varargout = {K, struct('load_chunks', @load_chunks, 'measures', @measures, 'boot', @boot)};
        otherwise, error('unknown mode %s', mode);
    end
end


function v = opt(c, i, d)
    if numel(c) >= i, v = c{i}; else, v = d; end
end


function logged(K, tag, fn)
    diary(fullfile(K.np, ['log_' tag '.txt'])); c = onCleanup(@() diary('off'));
    fprintf('\n===== NATURAL PARTNERS %s %s UTC =====\n', tag, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fn();
end


function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end


function c = satcon_config(K, H)
    cfg = H.load_configs(K); c = cfg{4};
    assert(strcmp(stored_value('sensingMode', c.sensingMode), 'contrast') && c.hillCoefficient == 1, 'config 4 is not saturating contrast');
    assert(c.numOpponents == 1 && strcmp(c.payoffMode, 'focal') && ~c.requireStability, 'coop_pair assumptions');
end


% =========================================================================
% POOL and mutants. Partners: 80 individuals drawn without replacement from the pool (half A the
% first 40 drawn, half B the last 40); they are not excluded from the network draw.
% =========================================================================
function do_pool(K, H)
    P = build_snapshots(K);
    [Wg, ~, gid] = unique(P.Wind, 'rows');
    fprintf('pool: %d replicates x %d snapshots x %d individuals = %d individuals; %d distinct genotypes\n', ...
        P.nRep, numel(K.gens), P.N, size(P.Wind, 1), size(Wg, 1));
    rng(K.seed + 9101, 'twister'); part = randperm(size(P.Wind, 1), K.nPart)';
    fprintf('seed partners K.seed + 9101 = %d (randperm, without replacement)\n', K.seed + 9101);
    [foc, nRe] = draw_focals(K, gid, size(P.Wind, 1));
    fprintf('seed networks K.seed + 9102 = %d; %d networks, %d draws redrawn (genotype already chosen)\n', K.seed + 9102, numel(foc), nRe);
    F = focal_table(K, P, foc); Pw = P.Wind(part, :);
    Pt = struct('ind', part, 'W', Pw, 'rep', P.rep(part), 'gen', P.gen(part));
    config = satcon_config(K, H); Pt.self = nan(K.nPart, 1);
    for j = 1:K.nPart
        mutate_with_refusal('reset'); Wj = reshape(Pw(j, :), 4, 4); Pt.self(j) = mutate_with_refusal('coop', Wj, Wj, config);
    end
    mutate_with_refusal('reset');
    report_pool(K, Pt, F, gid);
    save(fullfile(K.np, 'pool.mat'), 'Pt', 'F', 'K', '-v7.3');
    fprintf('saved %s\n', fullfile(K.np, 'pool.mat'));
end


function P = build_snapshots(K)
% Every individual of the census snapshots (gen >= 100, 5 per replicate, asserted to be generations
% 100, 200, 300, 376 and 500 with 1,000 individuals each), with its replicate, generation and the
% M and B of the first individual carrying its W in that snapshot (the mutation parent, as in the
% per-mutation census). Genotype = W = M .* B.
    d = load(fullfile(K.pub, 'saturating_contrast_100reps.mat'), 'results'); R = d.results; clear d;
    assert(numel(R) == 100, 'expected 100 replicates, found %d', numel(R));
    nR = numel(R); N = 1000; nS = numel(K.gens); n = nR * nS * N;
    P = struct('nRep', nR, 'N', N, 'Wind', zeros(n, 16), 'rep', zeros(n, 1), 'gen', zeros(n, 1), ...
        'Mp', zeros(4, 4, n), 'Bp', zeros(4, 4, n));
    o = 0;
    for r = 1:nR
        rep = R{r}; sg = double(rep.historyGenerations(:))';
        ks = find(sg >= 100); ks = ks(round(linspace(1, numel(ks), 5)));
        assert(isequal(sg(ks), K.gens), 'replicate %d snapshots %s', r, mat2str(sg(ks)));
        for k = ks
            M = double(rep.historyM(:, :, :, k)); B = double(rep.historyB(:, :, :, k));
            assert(size(M, 3) == N, 'replicate %d gen %d has %d individuals', r, sg(k), size(M, 3));
            Wall = reshape(M .* B, 16, N)'; [~, ~, ic] = unique(Wall, 'rows');
            firstIdx = accumarray(ic, (1:N)', [], @min); fi = firstIdx(ic);
            id = o + (1:N); P.Wind(id, :) = Wall; P.rep(id) = r; P.gen(id) = sg(k);
            P.Mp(:, :, id) = M(:, :, fi); P.Bp(:, :, id) = B(:, :, fi); o = o + N;
        end
    end
    fprintf('snapshots: generations %s in all %d replicates (asserted), %d individuals each\n', mat2str(K.gens), nR, N);
end


function [foc, nRe] = draw_focals(K, gid, n)
% Individuals drawn uniformly with replacement; a draw whose genotype is already chosen is redrawn,
% until K.nFocal. A network's replicate (the bootstrap cluster) is that of the drawn individual.
    rng(K.seed + 9102, 'twister'); foc = zeros(K.nFocal, 1); seen = false(max(gid), 1); nRe = 0; f = 0;
    while f < K.nFocal
        g = randi(n);
        if seen(gid(g)), nRe = nRe + 1; continue; end
        f = f + 1; foc(f) = g; seen(gid(g)) = true;
    end
end


function F = focal_table(K, P, foc)
    n = numel(foc); F = struct('ind', foc, 'W', P.Wind(foc, :), 'rep', P.rep(foc), 'gen', P.gen(foc), ...
        'Wmut', zeros(n, K.nMut, 16), 'slot', zeros(n, K.nMut), 'dW', zeros(n, K.nMut), ...
        'isTog', false(n, K.nMut), 'mutSeed', K.seed + 10000 + (1:n)');
    cfgK = struct('pToggle', 0.2, 'sdWeight', 0.2, 'K', K.nMut);
    for f = 1:n
        rng(F.mutSeed(f), 'twister');
        Mi = P.Mp(:, :, foc(f)); Bi = P.Bp(:, :, foc(f));
        assert(isequal(reshape(Mi .* Bi, 1, 16), F.W(f, :)), 'Mfirst .* Bfirst differs from the network''s W');
        [Wm, slot, dl, tog] = census_mutations(Mi, Bi, cfgK);
        F.Wmut(f, :, :) = reshape(Wm, 1, K.nMut, 16); F.slot(f, :) = slot; F.dW(f, :) = dl; F.isTog(f, :) = tog;
    end
end


function [Wout, slot, dWv, isTog] = census_mutations(Mi, Bi, cfg)
% The mutation distribution of the per-mutation census (per_mutation_census.m,
% wiring_metrics): one connection drawn uniformly, toggled with probability pToggle, else its weight
% changed by a normal draw of SD sdWeight (0.2 and 0.2 here); the draw lines are copied, and check 2
% tests the copy. Draws that leave W unchanged (a weight change on an absent connection) are kept,
% as in the census. Returns each mutant's W (K x 16), slot, dW, toggle.
    nG = 4; W = Mi .* Bi; K = cfg.K;
    Wout = zeros(K, 16); dWv = nan(K,1); slot = nan(K,1); isTog = false(K,1);
    for t = 1:K
        Mm = Mi; Bm = Bi;
        g1 = randi(nG); g2 = randi(nG);
        tog = rand() < cfg.pToggle;
        if tog
            Mm(g1,g2) = ~Mm(g1,g2);
        else
            Bm(g1,g2) = Bm(g1,g2) + randn()*cfg.sdWeight;
        end
        Wm = Mm .* Bm;
        q = sub2ind([nG nG], g1, g2);
        dl = Wm(g1,g2) - W(g1,g2);
        dWv(t) = dl; slot(t) = q; isTog(t) = tog; Wout(t, :) = reshape(Wm, 1, 16);
    end
end


function report_pool(K, Pt, F, gid)
    pg = gid(Pt.ind);
    fprintf('partners: %d individuals, %d distinct genotypes, from %d replicates; generations %s\n', K.nPart, ...
        numel(unique(pg)), numel(unique(Pt.rep)), mat2str(histcounts(Pt.gen, [K.gens Inf])));
    fprintf('  partners own-copy cooperation: >= 0.99 %d, in [0.5, 0.99) %d, < 0.5 %d; mean %.3f\n', ...
        sum(Pt.self >= 0.99), sum(Pt.self >= 0.5 & Pt.self < 0.99), sum(Pt.self < 0.5), mean(Pt.self));
    fprintf('  half A own-copy >= 0.99 / < 0.5: %d / %d; half B %d / %d\n', sum(Pt.self(1:40) >= 0.99), ...
        sum(Pt.self(1:40) < 0.5), sum(Pt.self(41:80) >= 0.99), sum(Pt.self(41:80) < 0.5));
    neutral = squeeze(all(F.Wmut == reshape(F.W, [], 1, 16), 3));
    fprintf('networks: %d distinct genotypes from %d replicates; per generation %s\n', numel(F.ind), ...
        numel(unique(F.rep)), mat2str(histcounts(F.gen, [K.gens Inf])));
    fprintf('mutants: %d (toggle share %.4f; W unchanged %d); mutation seeds K.seed + 10000 + f = %d..%d\n', ...
        numel(F.slot), mean(F.isTog(:)), sum(neutral(:)), F.mutSeed(1), F.mutSeed(end));
end


% =========================================================================
% CHECK 1 (instrument) and CHECK 2 (mutations)
% =========================================================================
function do_check1(K, H)
    S = load(fullfile(K.np, 'pool.mat'), 'Pt', 'F'); config = satcon_config(K, H);
    rng(K.seed + 9401, 'twister'); fi = randi(K.nFocal, 200, 1); pj = randi(K.nPart, 200, 1);
    fprintf('seed K.seed + 9401 = %d: 200 (network, partner) pairs\n', K.seed + 9401);
    v = nan(200, 4); same = false(200, 1);
    for q = 1:200
        W = reshape(S.F.W(fi(q), :), 4, 4); P = reshape(S.Pt.W(pj(q), :), 4, 4);
        mutate_with_refusal('reset'); [c1, b, fa] = mutate_with_refusal('coop_pair', W, P, config);
        mutate_with_refusal('reset'); c0 = mutate_with_refusal('coop', W, P, config);
        v(q, :) = [c0 c1 b isequal(W, P)]; same(q) = isequal(c0, c1) && isequal(mean(fa), c1);
    end
    mutate_with_refusal('reset');
    fprintf('  pairs with W = P (own-clone branch): %d; network mean coop %.4f; partner mean coop %.4f\n', ...
        sum(v(:, 4)), mean(v(:, 1)), mean(v(:, 3)));
    fprintf('  coop_pair network side equal to coop exactly: %d of 200; max |diff| %.3g\n', sum(same), max(abs(v(:, 1) - v(:, 2))));
    ok = all(same); fprintf('CHECK 1: %s\n', pf(ok)); save_check(K, 'check1', ok);
end


function save_check(K, name, ok)
    f = fullfile(K.np, 'checks.mat'); C = struct();
    if isfile(f), C = load(f); end
    C.(name) = ok; C.([name '_utc']) = char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss');
    save(f, '-struct', 'C');
end


function do_check2(K)
% Check 2: the copied mutation lines, replayed under the census's own seed, must reproduce 200
% stored census draws (slot, change and toggle). The census (rho_by_architecture.mat, architecture
% 4) seeded replicate r by rng(rngSeed + 1000*4 + r) in process_rep; membership_test saves and
% restores the RNG, and simulate_pair at its default start and alpha_and_beta draw nothing, so the
% mutation draws are the stream's only consumer: rows in stored order, 30 draws each. Replayed on
% the client; if that fails, on a pool worker (the census ran under parfor, whose stream type differs).
    Q = load(fullfile(K.pub, 'rho_by_architecture.mat'), 'W_all', 'metadata'); c = Q.metadata.cfg;
    assert(c.pToggle == 0.2 && c.sdWeight == 0.2 && c.K == 30 && c.rngSeed == 20260904, 'census cfg differs');
    rho = Q.W_all{4}; clear Q;
    S = load(fullfile(K.pub, 'retaliation_fixed_opponent', 'retaliation_fixed_opponent_census.mat'), 'out');
    mi = S.out.arch{4}.memInputs(:); clear S;
    rng(K.seed + 9402, 'twister'); pick = randperm(numel(mi) * 30, 200)'; [iq, jq] = ind2sub([numel(mi) 30], pick);
    fprintf('seed K.seed + 9402 = %d: 200 stored draws from %d member rows (W known); census seed %d + 4000 + r\n', ...
        K.seed + 9402, numel(mi), c.rngSeed);
    T = struct('row', [mi(iq).row]', 'draw', jq, 'W', {arrayfun(@(x) reshape(x.W, 1, 16), mi(iq), 'uni', 0)});
    par = parents_for(K, rho, T);
    [ok, n] = replay_check(rho, T, par, c, 'client');
    if ~ok
        p = gcp('nocreate'); if isempty(p), p = parpool('Processes', 1); end
        f = parfeval(p, @replay_check, 2, rho, T, par, c, 'worker'); [ok, n] = fetchOutputs(f);
        fprintf('%s', f.Diary);    % the worker's own lines (generator type, counts)
    end
    fprintf('CHECK 2: %s (%d of 200 draws equal in slot, change and toggle)\n', pf(ok), n);
    save_check(K, 'check2', ok);
end


function par = parents_for(K, rho, T)
% Mfirst / Bfirst of each picked row's W in its (replicate, generation) snapshot (census convention).
    d = load(fullfile(K.pub, 'saturating_contrast_100reps.mat'), 'results'); par = cell(numel(T.row), 2);
    for q = 1:numel(T.row)
        rr = rho(T.row(q)); rep = d.results{rr.rep}; k = find(double(rep.historyGenerations) == rr.gen, 1);
        M = double(rep.historyM(:, :, :, k)); B = double(rep.historyB(:, :, :, k));
        Wall = reshape(M .* B, 16, size(M, 3))'; i = find(ismember(Wall, T.W{q}, 'rows'), 1);
        assert(~isempty(i), 'row %d: W not in its snapshot', T.row(q));
        par(q, :) = {M(:, :, i), B(:, :, i)};
    end
end


function [ok, nOk] = replay_check(rho, T, par, c, where)
    reps = unique([rho(T.row).rep]); hit = false(numel(T.row), 1); nAll = 0; okAll = 0;
    cfgK = struct('pToggle', c.pToggle, 'sdWeight', c.sdWeight, 'K', c.K); s0 = rng;
    for r = reps
        rng(c.rngSeed + 1000 * 4 + r); rows = find([rho.rep] == r);
        for i = rows
            q = find(T.row == i)'; Mi = zeros(4); Bi = zeros(4);
            if ~isempty(q), Mi = par{q(1), 1}; Bi = par{q(1), 2}; end
            [~, slot, dl, tog] = census_mutations(Mi, Bi, cfgK);
            okST = isequal(slot, double(rho(i).slot(:))) && isequal(tog, logical(rho(i).isToggle(:)));
            nAll = nAll + 1; okAll = okAll + okST;
            for qq = q, j = T.draw(qq); hit(qq) = slot(j) == rho(i).slot(j) && tog(j) == rho(i).isToggle(j) && dl(j) == rho(i).dW(j); end
        end
    end
    g = RandStream.getGlobalStream; rng(s0);
    nOk = sum(hit); ok = nOk == numel(T.row);
    fprintf('  replay on the %s (generator %s): %d replicates, %d rows; slot and toggle of all 30 draws equal in %d rows\n', ...
        where, g.Type, numel(reps), nAll, okAll);
    fprintf('  picked draws equal (slot, dW exactly, toggle): %d of %d; toggles among them %d\n', nOk, numel(T.row), sum(arrayfun(@(i, j) rho(i).isToggle(j), T.row, T.draw)));
end


% =========================================================================
% SCORING ('np_timing_run' and 'np_score'): chunks of 25 networks, resumable
% =========================================================================
function [Z, ev, ws, wall] = score_focals(K, H, fset, dirC, nW)
    S = load(fullfile(K.np, 'pool.mat'), 'Pt', 'F'); config = satcon_config(K, H);
    if ~isfolder(dirC), mkdir(dirC); end
    nW = open_workers(nW); nC = ceil(numel(fset) / K.npChunk); Z = cell(nC, 1); t0 = tic;
    fprintf('encounter seed %d (mutate_with_refusal coop / coop_pair); %d networks in %d chunks\n', K.encSeed, numel(fset), nC);
    for k = 1:nC
        id = fset((k - 1) * K.npChunk + 1:min(k * K.npChunk, numel(fset)));
        f = fullfile(dirC, sprintf('chunk_%03d.mat', k)); want = chunk_inputs(K, S, id);
        if isfile(f)
            Zk = load(f); assert_chunk(Zk, want, f);
        else
            Zk = run_chunk(want, S.Pt.W, config); save(f, '-struct', 'Zk');
        end
        Z{k} = Zk;
        fprintf('  chunk %d/%d (%d networks) done, elapsed %.1f min\n', k, nC, numel(id), toc(t0) / 60);
    end
    wall = toc(t0); Z = cat_chunks(Z); ws = sum(Z.sec);
    ev = numel(fset) * (K.nPart + 1 + K.nMut * (K.nPart / 2 + 1));
    fprintf('  %d networks, %d evaluations, wall %.1f min on %d workers, worker-seconds %.0f; non-finite %d\n', ...
        numel(fset), ev, wall / 60, nW, ws, sum(~isfinite([Z.a(:); Z.am(:); Z.selfM(:); Z.self(:)])));
end


function in = chunk_inputs(K, S, id)
    in = struct('focal', id(:), 'Wf', S.F.W(id, :), 'Wmut', S.F.Wmut(id, :, :), 'Wpart', S.Pt.W, ...
        'mutSeed', S.F.mutSeed(id), 'encSeed', K.encSeed);
end


function assert_chunk(Z, want, f)
% A stored chunk is used only if its networks, mutants, partners and seeds match the pool.
    for fn = {'focal', 'Wf', 'Wmut', 'Wpart', 'mutSeed', 'encSeed'}
        assert(isfield(Z, fn{1}) && isequal(Z.(fn{1}), want.(fn{1})), 'chunk %s: %s does not match the pool (stale)', f, fn{1});
    end
end


function Z = run_chunk(in, Wpart, config)
    n = numel(in.focal); nP = size(Wpart, 1); nA = nP / 2; nM = size(in.Wmut, 2);
    [a, b] = deal(nan(n, nP)); [fa, pa] = deal(false(n, 100, nP)); self = nan(n, 1);
    selfM = nan(n, nM); am = nan(n, nM, nA); sec = nan(n, 1);
    P = arrayfun(@(j) reshape(Wpart(j, :), 4, 4), 1:nP, 'uni', 0); Wf = in.Wf; Wmut = in.Wmut;
    parfor i = 1:n
        t = tic; R = eval_focal(reshape(Wf(i, :), 4, 4), reshape(Wmut(i, :, :), nM, 16), P, nA, config);
        a(i, :) = R.a; b(i, :) = R.b; fa(i, :, :) = R.fa; pa(i, :, :) = R.pa; self(i) = R.self;
        selfM(i, :) = R.selfM; am(i, :, :) = R.am; sec(i) = toc(t);
    end
    Z = in; Z.a = a; Z.b = b; Z.fa = fa; Z.pa = pa; Z.self = self; Z.selfM = selfM; Z.am = am; Z.sec = sec;
end


function R = eval_focal(W, Wm, P, nA, config)
% 80 coop_pair (parent vs every partner, with the focal and the partner action of the same
% encounter), parent own copy, 30 mutants x (half A + own copy); 100 encounters per pair under the
% runs' encounter rule (mutate_with_refusal).
    mutate_with_refusal('reset'); nP = numel(P); nM = size(Wm, 1);
    R = struct('a', nan(1, nP), 'b', nan(1, nP), 'fa', false(100, nP), 'pa', false(100, nP), ...
        'self', NaN, 'selfM', nan(1, nM), 'am', nan(nM, nA));
    for j = 1:nP
        [R.a(j), R.b(j), fa, pa] = mutate_with_refusal('coop_pair', W, P{j}, config);
        R.fa(:, j) = fa(:); R.pa(:, j) = pa(:);
    end
    R.self = mutate_with_refusal('coop', W, W, config);
    for t = 1:nM
        Wt = reshape(Wm(t, :), 4, 4); R.selfM(t) = mutate_with_refusal('coop', Wt, Wt, config);
        for j = 1:nA, R.am(t, j) = mutate_with_refusal('coop', Wt, P{j}, config); end
    end
    mutate_with_refusal('reset');
end


function Z = cat_chunks(C)
    Z = C{1}; f = fieldnames(Z);
    for k = 2:numel(C)
        for i = 1:numel(f)
            if ismember(f{i}, {'Wpart', 'encSeed'}), continue; end
            Z.(f{i}) = cat(1, Z.(f{i}), C{k}.(f{i}));
        end
    end
end


function do_timing_run(K, H, nW)
    [~, ev, ws, wall] = score_focals(K, H, (1:K.nPilot)', fullfile(K.np, 'timing_run'), nW);
    hW = wall / K.nPilot * K.nFocal / 3600; hS = ws / nW / K.nPilot * K.nFocal / 3600; h = max(hW, hS);
    fprintf('  seconds per evaluation: wall %.4f, per worker %.4f\n', wall / ev, ws / ev);
    fprintf('  projection for %d networks on %d workers: %.2f h (wall-scaled), %.2f h (worker-seconds / %d); read %.2f h\n', ...
        K.nFocal, nW, hW, hS, nW, h);
    fprintf('  on 4 workers (scaled): %.2f h; on 3: %.2f h\n', h * nW / 4, h * nW / 3);
    fprintf('  projected time: %.2f h\n', h);
    T = struct('wall', wall, 'ev', ev, 'ws', ws, 'nW', nW, 'h', h);
    save(fullfile(K.np, 'timing_run_manifest.mat'), 'T');
end


function do_score(K, H, nW)
    C = load(fullfile(K.np, 'checks.mat'));
    assert(isfield(C, 'check1') && C.check1 && isfield(C, 'check2') && C.check2, ...
        'check 1 or check 2 failed; not scoring');
    score_focals(K, H, (1:K.nFocal)', fullfile(K.np, 'chunks'), nW);
end


% =========================================================================
% ANALYSIS: measures
% =========================================================================
function A = measures(K, Z, F)
% Per network: reaction R (all 80 partners; the dial keeps fifths 2 to 5 of it), R_A / R_B (the
% two halves), Q (mean cooperation over the partners), M (share of its 30 mutants x 40 half-A
% partners whose cooperation falls on the other side of 0.5 from the parent's), mean |da|, own-copy
% cooperation, L (mutants whose own-copy cooperation is below 0.5), R_pair (R from pair means).
    A.rep = F.rep(Z.focal); A.focal = Z.focal; n = numel(A.focal); nA = K.nPart / 2;
    [A.R, A.Rdef] = reaction(Z.fa, Z.pa, 1:K.nPart);
    [A.RA, A.RAdef] = reaction(Z.fa, Z.pa, 1:nA);
    [A.RB, A.RBdef] = reaction(Z.fa, Z.pa, nA + 1:K.nPart);
    A.Q = mean(Z.a, 2);
    aP = reshape(Z.a(:, 1:nA), n, 1, nA);                          % parent on half A
    flip = (Z.am >= 0.5) ~= (aP >= 0.5);                           % n x 30 x 40, side of 0.5 changed
    A.M = mean(mean(flip, 3), 2);
    A.dabs = mean(mean(abs(Z.am - aP), 3), 2);
    A.self = Z.self; A.lost = Z.selfM < 0.5; A.L = sum(A.lost, 2); A.coopSelf = Z.self >= 0.99;
    a = Z.a; b = Z.b;
    A.Rpair = sum(a .* b, 2) ./ sum(b, 2) - sum(a .* (1 - b), 2) ./ sum(1 - b, 2);
    A.vary = reshape(any(Z.fa, 2) & ~all(Z.fa, 2), n, []);        % n x 80: focal action varies in the pair
    A.a = a; A.b = b;
end


function [R, def] = reaction(fa, pa, js)
% P(i coop | partner coop in that encounter) - P(i coop | partner defected); defined if each
% condition has >= 100 encounters among partners js.
    fa = fa(:, :, js); pa = pa(:, :, js);
    nC = sum(sum(pa, 2), 3); nD = sum(sum(~pa, 2), 3);
    R = sum(sum(fa & pa, 2), 3) ./ nC - sum(sum(fa & ~pa, 2), 3) ./ nD;
    def = nC >= 100 & nD >= 100; R(~def) = NaN;
end


% =========================================================================
% ANALYSIS: statistics
% =========================================================================
function r = sp(x, y)
    r = corr(x(:), y(:), 'Type', 'Spearman');
end


function r = sp_partial(x, y, z)
    rxy = sp(x, y); rxz = sp(x, z); ryz = sp(y, z);
    r = (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2));
end


function [est, ci, nd] = boot(K, rep, fn, seed)
% fn(idx) on the network set; K.nBoot resamples of the replicates in the set (networks carry their
% replicate's multiplicity, each network unweighted); percentile 95%; nd = resamples dropped as
% non-finite.
    est = fn((1:numel(rep))'); [u, ~, loc] = unique(rep); v = nan(K.nBoot, 1);
    rs = rng; rng(seed, 'twister');
    for b = 1:K.nBoot
        mult = accumarray(randi(numel(u), numel(u), 1), 1, [numel(u) 1]);
        v(b) = fn(repelem((1:numel(rep))', mult(loc)));
    end
    rng(rs); ci = prctile(v(isfinite(v)), [2.5 97.5]); nd = sum(~isfinite(v));
end


function T = n3_rows(A, s)
% Rows for H.boot_fit: one per network, weight 1, k of 30 own-copy losses, R standardized once
% over the set (mean / SD).
    R = A.R(s); T = struct('arch', 4 * ones(sum(s), 1), 'rep', A.rep(s), 'count', ones(sum(s), 1), ...
        'lost', A.lost(s, :), 'score', (R - mean(R)) / std(R));
end


function [o, ci, nd] = or_sd(K, H, T, seed)
    Fi = H.boot_fit(T, {@(t) [ones(numel(t.arch), 1), t.score]}, K, seed);
    v = exp(Fi.B{1}(:, 2)); nd = sum(~isfinite(v));
    o = exp(Fi.b{1}(2)); ci = prctile(v(isfinite(v)), [2.5 97.5]);
end


% =========================================================================
% ANALYSIS: checks 5, 3, 4, measures 1 to 3 and the descriptive tables
% =========================================================================
function do_analyze(K, H, src)
    if strcmp(src, 'pilot'), src = 'timing_run'; end   % the earlier name, still accepted
    isPilot = strcmp(src, 'timing_run'); S = load(fullfile(K.np, 'pool.mat'), 'Pt', 'F');
    if isPilot, fset = (1:K.nPilot)'; dirC = fullfile(K.np, 'timing_run'); name = 'natural_partners_timing_run.txt';
    else, fset = (1:K.nFocal)'; dirC = fullfile(K.np, 'chunks'); name = 'natural_partners.txt'; end
    Z = load_chunks(K, S, fset, dirC); A = measures(K, Z, S.F);
    C = load(fullfile(K.np, 'checks.mat')); f = fullfile(K.np, name);

    %% Report
    report = @() analyze_report(K, H, A, S.Pt, C, src, numel(fset), isPilot, f); %#ok<NASGU> called inside evalc
    rep = evalc('report()');
    fprintf('%s', rep);

    %% Save
    fid = fopen(f, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function analyze_report(K, H, A, Pt, C, src, n, isPilot, f)
% Checks 1 and 2 (stored), 5, 4 and 3, then measures 1 to 3, the descriptive tables and the calibration of
% check 5; prints the report's lines.
    fprintf('===== NATURAL PARTNERS analyze (%s) %s UTC =====\n', src, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('saturating contrast; %d networks, %d partners\n', n, K.nPart);
    if isPilot, fprintf('trial run: the first %d networks\n', K.nPilot); end
    ok1 = isfield(C, 'check1') && C.check1; ok2 = isfield(C, 'check2') && C.check2;
    fprintf('Checks 1 and 2 (stored in checks.mat): %s, %s\n', pf(ok1), pf(ok2));
    s1 = A.Rdef; s3 = A.coopSelf & A.Rdef;
    check5(K, H, A, s1, s3);
    checks34(A, isPilot);
    readings(K, H, A, s1, s3);
    not_read(K, A, s1, Pt);
    calibrate(K, H, A, s1, s3);
    fprintf('\nReport written to %s\n', f);
end


function Z = load_chunks(K, S, fset, dirC)
    nC = ceil(numel(fset) / K.npChunk); C = cell(nC, 1);
    for k = 1:nC
        id = fset((k - 1) * K.npChunk + 1:min(k * K.npChunk, numel(fset)));
        f = fullfile(dirC, sprintf('chunk_%03d.mat', k)); assert(isfile(f), 'missing %s', f);
        C{k} = load(f); assert_chunk(C{k}, chunk_inputs(K, S, id), f);
    end
    Z = cat_chunks(C);
end


function check5(K, H, A, s1, s3)
% Can measures 1 and 3 return the other answer: (a) M permuted across all networks, and (a') M
% shuffled within replicate (the measure 1 interval must hold 0); (b) M planted at Spearman 0.3 (0.3 must lie in
% the interval); (c) losses planted at OR 1.5 per SD of R.
    fprintf('\nCHECK 5 (can measures 1 and 3 return the other answer)\n');
    R = A.R(s1); M = A.M(s1); rp = A.rep(s1);
    Ma = shuffle_all(M, K.seed + 9504);
    [e, ci, nd] = boot(K, rp, @(i) sp(R(i), Ma(i)), K.seed + 9614); okA = ci(1) <= 0 && ci(2) >= 0;
    fprintf('  (a) M permuted across all networks (seed %d): measure 1 %.3f [%.3f, %.3f], %d resamples dropped -> %s (pass needs 0 inside the interval)\n', ...
        K.seed + 9504, e, ci, nd, pf(okA));
    Mw = shuffle_within(M, rp, K.seed + 9501);
    [e, ci, nd] = boot(K, rp, @(i) sp(R(i), Mw(i)), K.seed + 9611); okW = ci(1) <= 0 && ci(2) >= 0;
    fprintf('  (a'') M shuffled within replicate (seed %d): measure 1 %.3f [%.3f, %.3f], %d resamples dropped -> %s (pass needs 0 inside the interval)\n', ...
        K.seed + 9501, e, ci, nd, pf(okW));
    Mp = plant_m(R, K.seed + 9502);
    [e, ci, nd] = boot(K, rp, @(i) sp(R(i), Mp(i)), K.seed + 9612); okB = ci(1) <= 0.3 && ci(2) >= 0.3;
    fprintf('  (b) planted Spearman 0.3 (seed %d, copula rho %.4f): measure 1 %.3f [%.3f, %.3f], %d resamples dropped -> %s (pass needs 0.3 inside the interval)\n', ...
        K.seed + 9502, 2 * sin(pi * 0.3 / 6), e, ci, nd, pf(okB));
    okC = false;
    if sum(s3) >= 2
        T = H.plant(n3_rows(A, s3), 1.5, K.seed + 9503, 1); [o, ci, nd] = or_sd(K, H, T, K.seed + 9613);
        okC = ci(1) <= 1.5 && ci(2) >= 1.5;
        fprintf('  (c) planted OR 1.5 per SD (seed %d) on %d networks: OR %.3f [%.3f, %.3f], %d resamples dropped -> %s (pass needs 1.5 inside the interval)\n', ...
            K.seed + 9503, sum(s3), o, ci, nd, pf(okC));
    else, fprintf('  (c) fewer than 2 networks for measure 3: FAIL\n');
    end
    fprintf('  CHECK 5 with (a): %s; with (a'') in place of (a): %s\n', pf(okA && okB && okC), pf(okW && okB && okC));
end


function x = shuffle_all(M, seed)
    rs = rng; rng(seed, 'twister'); x = M(randperm(numel(M))); rng(rs);
end


function x = shuffle_within(M, rp, seed)
% One randperm per replicate, replicates in sorted order.
    rs = rng; rng(seed, 'twister'); x = M;
    for r = unique(rp)', i = find(rp == r); x(i) = M(i(randperm(numel(i)))); end
    rng(rs);
end


function Mp = plant_m(R, seed)
% Gaussian copula with Spearman 0.3 on the normal scores of tiedrank(R):
% M = rho_g z + sqrt(1 - rho_g^2) e, rho_g = 2 sin(pi 0.3 / 6).
    rs = rng; rng(seed, 'twister'); rg = 2 * sin(pi * 0.3 / 6);
    z = norminv(tiedrank(R) / (numel(R) + 1)); Mp = rg * z + sqrt(1 - rg^2) * randn(size(z)); rng(rs);
end


function calibrate(K, H, A, s1, s3)
% Calibration of check 5, printed only (no paper number uses it): each part over 100 seeds x 500 resamples (seeds
% K.seed + 30000 + offset + s), the share of seeds whose interval excludes 0 or holds the planted value.
    Kc = K; Kc.nBoot = 500; nS = 100; b0 = K.seed + 30000; R = A.R(s1); M = A.M(s1); rp = A.rep(s1);
    fprintf('\nCALIBRATION OF CHECK 5 (reported, not read; %d seeds x %d resamples; seeds %d..%d)\n', nS, Kc.nBoot, b0 + 1, b0 + 7000 + nS);
    t = tic; x = false(nS, 4); ndT = 0;
    for s = 1:nS
        Ma = shuffle_all(M, b0 + s); Mw = shuffle_within(M, rp, b0 + 2000 + s); Mp = plant_m(R, b0 + 4000 + s);
        [~, ci, n1] = boot(Kc, rp, @(i) sp(R(i), Ma(i)), b0 + 1000 + s); x(s, 1) = ci(1) > 0 || ci(2) < 0;
        [~, ci, n2] = boot(Kc, rp, @(i) sp(R(i), Mw(i)), b0 + 3000 + s); x(s, 2) = ci(1) > 0 || ci(2) < 0;
        [~, ci, n3] = boot(Kc, rp, @(i) sp(R(i), Mp(i)), b0 + 5000 + s); x(s, 3) = ci(1) <= 0.3 && ci(2) >= 0.3;
        ndT = ndT + n1 + n2 + n3;
    end
    fprintf('  (a) across-all shuffle: measure 1 interval excludes 0 in %d of %d seeds (%.2f)\n', sum(x(:, 1)), nS, mean(x(:, 1)));
    fprintf('  (a) within-replicate shuffle: measure 1 interval excludes 0 in %d of %d seeds (%.2f)\n', sum(x(:, 2)), nS, mean(x(:, 2)));
    fprintf('  (b) planted Spearman 0.3: interval holds 0.3 in %d of %d seeds (%.2f)\n', sum(x(:, 3)), nS, mean(x(:, 3)));
    fprintf('  (a), (b): %d of %d resamples dropped as non-finite; %.1f s\n', ndT, 3 * nS * Kc.nBoot, toc(t));
    if sum(s3) < 2, fprintf('  (c) fewer than 2 networks for measure 3: not run\n'); return; end
    t = tic; ndC = 0; T0 = n3_rows(A, s3);
    for s = 1:nS
        [~, ci, n4] = or_sd(Kc, H, H.plant(T0, 1.5, b0 + 6000 + s, 1), b0 + 7000 + s);
        x(s, 4) = ci(1) <= 1.5 && ci(2) >= 1.5; ndC = ndC + n4;
    end
    fprintf('  (c) planted OR 1.5 per SD: interval holds 1.5 in %d of %d seeds (%.2f); %d of %d resamples dropped; %.1f s\n', ...
        sum(x(:, 4)), nS, mean(x(:, 4)), ndC, nS * Kc.nBoot, toc(t));
end


function checks34(A, isPilot)
    n = numel(A.R); nd = sum(A.Rdef);
    if isPilot, ok4 = nd >= 0.1 * n; rule = sprintf('>= 10%% of %d (trial run, pro rata)', n);
    else, ok4 = nd >= 50; rule = '>= 50 of 500'; end
    fprintf('\nCHECK 4 (variation): R defined in %d of %d networks (%.3f); rule %s -> %s\n', nd, n, nd / n, rule, pf(ok4));
    b = A.RAdef & A.RBdef; r = NaN; if sum(b) >= 3, r = sp(A.RA(b), A.RB(b)); end
    ok3 = isfinite(r) && r >= 0.5;
    fprintf('CHECK 3 (reliability): Spearman R(half A) vs R(half B) = %.3f over %d networks (defined: half A %d, half B %d) -> %s (>= 0.5)\n', ...
        r, sum(b), sum(A.RAdef), sum(A.RBdef), pf(ok3));
end


function readings(K, H, A, s1, s3)
% Measure 1 Spearman(R, M), measure 2 its rank partial on Q, measure 3 the odds ratio of own-copy
% losses per SD of R.
    R = A.R(s1); M = A.M(s1); Q = A.Q(s1); rp = A.rep(s1);
    fprintf('\nMEASURES 1 TO 3 (%d networks with R defined; intervals: %d resamples of %d replicates)\n', sum(s1), K.nBoot, numel(unique(rp)));
    [e1, c1, d1] = boot(K, rp, @(i) sp(R(i), M(i)), K.seed + 9601);
    fprintf('  measure 1 Spearman(R, M):     %.3f [%.3f, %.3f], %d resamples dropped\n', e1, c1, d1);
    [e2, c2, d2] = boot(K, rp, @(i) sp_partial(R(i), M(i), Q(i)), K.seed + 9601);
    fprintf('  measure 2 rank partial r(R, M | Q):  %.3f [%.3f, %.3f], %d resamples dropped\n', e2, c2, d2);
    if sum(s3) >= 2
        [o, ci, d3] = or_sd(K, H, n3_rows(A, s3), K.seed + 9603);
        fprintf('  measure 3 logit L (k of 30) on R per SD, %d networks (own copy >= 0.99, R defined): OR %.3f [%.3f, %.3f], %d resamples dropped\n', ...
            sum(s3), o, ci, d3);
    else, fprintf('  measure 3: fewer than 2 networks\n'); end
end


% =========================================================================
% ANALYSIS: printed only, no paper number uses it
% =========================================================================
function not_read(K, A, s1, Pt)
    fprintf('\nalso reported\n');
    fprintf('  R defined: %d of %d (%.3f)\n', sum(s1), numel(s1), mean(s1));
    dist('R (defined)', A.R(s1)); dist('M (R defined)', A.M(s1)); dist('M (all)', A.M);
    dist('Q (all)', A.Q); dist('Q (R defined)', A.Q(s1)); dist('mean |a(mut) - a(par)| (all)', A.dabs);
    dist('own-copy cooperation (all)', A.self); dist('L of 30, own copy >= 0.99', A.L(A.coopSelf));
    fprintf('  networks with own copy >= 0.99: %d; with R defined too: %d\n', sum(A.coopSelf), sum(A.coopSelf & s1));
    fprintf('  partners'' mean cooperation against the networks: overall %.3f\n', mean(A.b(:)));
    dist('  per partner (mean over networks)', mean(A.b, 1)');
    dist('  networks'' cooperation per partner (mean over networks)', mean(A.a, 1)');
    fprintf('  partners'' own-copy cooperation: mean %.3f, >= 0.99 %d, < 0.5 %d\n', mean(Pt.self), sum(Pt.self >= 0.99), sum(Pt.self < 0.5));
    quintiles(A, s1);
    within_coop(K, A, s1);
    within_pair(K, A, s1);
    R = A.R(s1); d = A.dabs(s1); [e, ci, nd] = boot(K, A.rep(s1), @(i) sp(R(i), d(i)), K.seed + 9604);
    fprintf('  Spearman(R, mean |a(mut) - a(par)|): %.3f [%.3f, %.3f], %d resamples dropped\n', e, ci, nd);
end


function dist(name, x)
    x = x(isfinite(x));
    if isempty(x), fprintf('  %-46s none\n', name); return; end
    q = prctile(x, [0 10 25 50 75 90 100]);
    fprintf('  %-46s n %4d mean %.4f | min %.4f p10 %.4f p25 %.4f med %.4f p75 %.4f p90 %.4f max %.4f\n', name, numel(x), mean(x), q);
end


function quintiles(A, s1)
    R = A.R(s1); M = A.M(s1); Qd = A.Q(s1); ed = unique(prctile(R, 0:20:100)); ed(end) = Inf;
    fprintf('  M by R quintile (R defined; edges by prctile, ties merged):\n');
    for b = 1:numel(ed) - 1
        k = R >= ed(b) & R < ed(b + 1);
        fprintf('    R [%7.4f, %7.4f)  n %4d  mean M %.4f  mean Q %.4f\n', ed(b), ed(b + 1), sum(k), mean(M(k)), mean(Qd(k)));
    end
end


function within_coop(K, A, s1)
    s = s1 & A.coopSelf; R = A.R(s); M = A.M(s);
    if sum(s) < 3, fprintf('  measure 1 within own copy >= 0.99: fewer than 3 networks\n'); return; end
    [e, ci, nd] = boot(K, A.rep(s), @(i) sp(R(i), M(i)), K.seed + 9605);
    fprintf('  measure 1 within networks with own copy >= 0.99 (%d): %.3f [%.3f, %.3f], %d resamples dropped\n', sum(s), e, ci, nd);
end


function within_pair(K, A, s1)
% The within-pair part of R: R_pair from pair means, beside R; "network action varies" = a pair's
% 100 network actions are not all equal.
    R = A.R(s1); Rp = A.Rpair(s1); M = A.M(s1); rp = A.rep(s1);
    fprintf('  Within-pair part of R (R_pair from pair means, %d networks with R defined):\n', sum(s1));
    dist('    R_pair', Rp); fprintf('    mean R - R_pair: %.4f\n', mean(R - Rp));
    [e, ci, nd] = boot(K, rp, @(i) sp(R(i), Rp(i)), K.seed + 9606);
    fprintf('    Spearman(R, R_pair): %.3f [%.3f, %.3f], %d resamples dropped\n', e, ci, nd);
    [e, ci, nd] = boot(K, rp, @(i) sp(Rp(i), M(i)), K.seed + 9607);
    fprintf('    measure 1 with R_pair in place of R: %.3f [%.3f, %.3f], %d resamples dropped\n', e, ci, nd);
    v = A.vary(s1, :);
    fprintf('    pairs whose network action varies across the 100 encounters: %d of %d (%.4f); all networks %.4f\n', ...
        sum(v(:)), numel(v), mean(v(:)), mean(A.vary(:)));
end
