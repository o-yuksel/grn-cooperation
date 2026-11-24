function out = measure_defector_encounter_frequency(cfg)
% MEASURE_DEFECTOR_ENCOUNTER_FREQUENCY  At generation 500 of the four signaling architectures, how
% often does a member meet a partner that defects (f_D), and how does the selection on
% retaliation it implies, s_ret = f_D (P - S) / wbar with wbar the mean payoff, compare with the
% rate u at which mutation removes retaliation?
%
% Produces : main text, "How cooperating populations lose cooperation"; Methods, "Failure route";
%            SI S4.1 (median f_D 0.36%) and Table S9, top rows (f_D, u, s_ret / u); prints them
%            itself, in the blocks labeled 'encounter rule' ([1.2] and [4])
% Inputs   : cfg (optional struct; defaults in fill_defaults; 'reportOnly', true reprints the
%              report from the saved .mat)
%            publication/<architecture>_100reps.mat   (written by run_signaling_architectures)
%            publication/member_share_by_generation_standard100.mat   (optional, for the check of block
%              [7]; written by share_of_members_by_generation)
% Outputs  : publication/defector_encounter_frequency.mat and .txt (the printed report)
% Run      : measure_defector_encounter_frequency()   (about 20 min on 6 workers; maxReps 3 on
%            3 workers took about 3 min; the default pool is 3 workers)
% Seed     : rng(20260926) (cfg.rngSeed, the session's generator type); job j reseeds with
%            20260926 + j before its random starts; membership_test seeds each wiring itself
%
% Each individual meets one partner drawn uniformly from the other N - 1, so f_D is computed
% over every possible partner from the count of each distinct wiring, averaged over members by
% individuals; medians are over successful replicates (final cooperation at least 0.95). 'encounter rule'
% is the runs' encounter rule (start U(0, 0.01), 300 steps, 10 encounters per ordered pair);
% 'det600' develops each pair as the membership test does (start 0.005, 600 steps); 'rand600'
% and 'det300' change one of the two at a time. u takes fixed values measured by other scripts
% (the local function retaliation_removal_rates at the end of this file).

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg = fill_defaults(cfg);
    matFile = fullfile(cfg.outDir, 'defector_encounter_frequency.mat');
    if cfg.reportOnly
        S = load(matFile, 'out', 'cfg'); out = S.out; cfgRun = S.cfg;
    else
        out = run_all(cfg); cfgRun = cfg;
        save(matFile, 'out', 'cfg');
    end
    report_file = fullfile(cfg.outDir, 'defector_encounter_frequency.txt');
    if exist(report_file, 'file'), delete(report_file); end
    diary(report_file);
    print_report(out, cfgRun);
    print_membership_test_check(out, cfgRun);
    diary off;
end


% =========================================================================
function out = run_all(cfg)
    rng(cfg.rngSeed);
    ARCH = {'graded separate',     'graded_separate_100reps.mat'
            'graded contrast',     'graded_contrast_100reps.mat'
            'saturating separate', 'saturating_separate_100reps.mat'
            'saturating contrast', 'saturating_contrast_100reps.mat'};
    [jobs, meta] = build_jobs(ARCH, cfg);
    fprintf('%d snapshot jobs\n', numel(jobs));
    if isempty(gcp('nocreate')), parpool('Processes', cfg.nWorkers); end
    res = cell(numel(jobs), 1);
    t0 = tic;
    parfor j = 1:numel(jobs)
        rng(cfg.rngSeed + j);
        res{j} = analyse_snapshot(jobs{j}, cfg.engineDraws);
    end
    fprintf('analysis done in %.1f min\n', toc(t0) / 60);
    out = assemble(res, jobs, meta, ARCH);
end


% =========================================================================
function print_membership_test_check(out, cfg)
% Block [7], a check of the membership test: the mean generation-500 share of members per
% architecture and, replicate by replicate, the largest difference from the shares stored in
% member_share_by_generation_standard100.mat (saturating contrast, the same membership test).
    fprintf('\n[7] MEMBERSHIP TEST CHECK: gen-500 share of members\n');
    for a = 1:numel(out)
        fprintf('  %-20s mean share of members %.1f%%\n', out(a).name, 100*mean([out(a).at500.occ]));
    end
    f = fullfile(cfg.dataDir, 'member_share_by_generation_standard100.mat');
    if exist(f, 'file')
        H = load(f, 'out'); ref = H.out(1).occ(:, end);
        mine = [out(4).at500.occ]'; n = min(numel(ref), numel(mine));
        fprintf('  saturating contrast vs %s: max |diff| over %d replicates = %.2e\n', ...
            'member_share_by_generation_standard100.mat', n, max(abs(ref(1:n) - mine(1:n))));
    end
end


% =========================================================================
function cfg = fill_defaults(cfg)
    def = struct('dataDir', fullfile(grn_data_file('publication')), ...
                 'outDir',  fullfile(grn_data_file('publication')), ...
                 'rngSeed', 20260926, 'engineDraws', 10, 'windowStart', 300, ...
                 'successThreshold', 0.95, 'nWorkers', 3, 'maxReps', Inf, 'reportOnly', false);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
    end
end


% =========================================================================
function [jobs, meta] = build_jobs(ARCH, cfg)
% One job per (replicate, snapshot): every replicate at its last stored generation (500;
% successful replicates are the estimate, failed ones the control), and successful
% replicates also at every stored generation in [windowStart, 500] (stationarity, block [5]).
% A snapshot is the population as stored before it plays that generation.
    jobs = {}; meta = struct([]);
    for a = 1:size(ARCH, 1)
        f = fullfile(cfg.dataDir, ARCH{a, 2});
        d = load(f, 'results', 'config');
        info = dir(f);
        c = d.config;
        sp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, ...
            'discreteOutput', c.discreteOutput, 'michaelisConstant', c.michaelisConstant, ...
            'degradationRate', c.degradationRate, 'timeStep', c.timeStep, ...
            'developmentSteps', 600, 'payoffDC', 0, 'sensorIntegrationRate', c.sensorIntegrationRate, ...
            'backend', 'simulate_pair');
        pay = [c.payoffDC c.payoffCC c.payoffDD c.payoffCD];   % [T R P S]
        nR = min(numel(d.results), cfg.maxReps);
        d.results = d.results(1:nR);
        fc = cellfun(@(r) double(r.historyCooperation(end)), d.results);
        sg = double(d.results{1}.historyGenerations(:))';
        meta(a).name = ARCH{a, 1}; meta(a).file = ARCH{a, 2}; meta(a).fileDate = info.date;
        meta(a).finalCoop = fc; meta(a).success = fc >= cfg.successThreshold;
        meta(a).gens = sg; meta(a).payoffs = pay; meta(a).nOpp = c.numOpponents;
        meta(a).payoffMode = c.payoffMode; meta(a).engineSteps = c.developmentSteps;
        meta(a).histCoop = cell2mat(cellfun(@(r) double(r.historyCooperation(:))', d.results, 'uni', 0));
        for r = 1:nR
            ks = numel(sg);
            if meta(a).success(r), ks = find(sg >= cfg.windowStart); end
            for k = ks
                jobs{end+1} = struct('arch', a, 'rep', r, 'gen', sg(k), ...
                    'M', d.results{r}.historyM(:,:,:,k), 'B', d.results{r}.historyB(:,:,:,k), ...
                    'sp', sp, 'pay', pay, 'encounter_rule', sg(k) == 500, ...
                    'engineSteps', c.developmentSteps, 'histCoop', double(d.results{r}.historyCooperation(sg(k)))); %#ok<AGROW>
            end
        end
    end
end


% =========================================================================
function r = analyse_snapshot(job, nDraws)
% One snapshot: its distinct wirings with their counts c, membership of each (membership_test), and
% the encounter tables, first as the membership test develops a pair (det600), then under the
% runs' encounter rule. Pm(u, v) is the chance that a focal of wiring u meets wiring v.
    N = size(job.M, 3);
    Wall = reshape(double(job.M) .* double(job.B), 16, N)';
    [Wu, ~, ic] = unique(Wall, 'rows');
    c = accumarray(ic, 1);
    nu = size(Wu, 1);
    mem = false(nu, 1);
    for u = 1:nu
        mem(u) = membership_test(reshape(Wu(u, :), 4, 4), job.sp).is_css;
    end
    [A1, A2, R2] = pair_table(Wu, job.sp, 1, []);
    Pm = (repmat(c', nu, 1) - eye(nu)) / (N - 1);   % P(focal u meets v)
    r = summarise(A1, A2, Pm, c, mem, job.pay);
    r.nUnique = nu;
    r.symMismatch = sum(sum(A2 ~= A1'));
    r.memTieShare = weighted_member_mean(Pm .* double(R2 == 0), c, mem);
    r.memResidueShare = weighted_member_mean(Pm .* double(abs(R2) < 1e-12 & ~A2), c, mem);
    r.memDefectorIsMember = defect_from_members(A2, Pm, c, mem);
    spE = job.sp; spE.developmentSteps = job.engineSteps;
    [E1, E2] = pair_table(Wu, spE, nDraws, 0.01);
    r.eng = summarise(E1, E2, Pm, c, mem, job.pay);
    r.eng.memDefectorIsMember = defect_from_members(E2, Pm, c, mem);
    r.rand600 = struct(); r.det300 = struct();
    if job.encounter_rule
        % At generation 500, the two differences (start, steps) one at a time
        [X1, X2] = pair_table(Wu, job.sp, nDraws, 0.01);
        r.rand600 = summarise(X1, X2, Pm, c, mem, job.pay);
        [Y1, Y2] = pair_table(Wu, spE, 1, []);
        r.det300 = summarise(Y1, Y2, Pm, c, mem, job.pay);
    end
end


% =========================================================================
function [A1, A2, R2] = pair_table(Wu, sp, nDraws, startScale)
% A1(u,v) / A2(u,v): share of nDraws encounters in which agent 1 (wiring u) / agent 2 (wiring v)
% cooperates when u meets v. startScale empty: simulate_pair's fixed start (0.005);
% otherwise each agent starts from U(0, startScale) in every draw. R2(u,v): agent 2's readout
% x_C - x_D (fixed start only).
    nu = size(Wu, 1); A1 = zeros(nu); A2 = zeros(nu); R2 = nan(nu);
    E = ones(4);
    for u = 1:nu
        Wa = reshape(Wu(u, :), 4, 4);
        for v = 1:nu
            Wb = reshape(Wu(v, :), 4, 4);
            for k = 1:nDraws
                if ~isempty(startScale)
                    sp.initialState1 = rand(4, 1) * startScale;
                    sp.initialState2 = rand(4, 1) * startScale;
                end
                [oc, fs] = simulate_pair(E, Wa, [], E, Wb, [], sp);
                A1(u, v) = A1(u, v) + (oc(1) == 'C') / nDraws;
                A2(u, v) = A2(u, v) + (oc(2) == 'C') / nDraws;
                if isempty(startScale), R2(u, v) = fs(5) - fs(6); end
            end
        end
    end
end


% =========================================================================
function r = summarise(A1, A2, Pm, c, mem, pay)
% Expected payoffs from the tables; f_D = clone-weighted member mean of P(partner defects);
% s = f_D (P - S) / wbar, with wbar the population's expected focal payoff (SI S4.1).
    T = pay(1); R = pay(2); P = pay(3); S = pay(4);
    D2 = 1 - A2;
    payoff = R * A1 .* A2 + S * A1 .* D2 + T * (1 - A1) .* A2 + P * (1 - A1) .* D2;
    N = sum(c);
    r.occ = sum(c(mem)) / N;
    r.nMem = sum(c(mem));
    r.fD = weighted_member_mean(Pm .* D2, c, mem);
    r.fExploited = weighted_member_mean(Pm .* A1 .* D2, c, mem);    % member C, partner D
    r.fRetaliated = weighted_member_mean(Pm .* (1 - A1) .* D2, c, mem);
    r.fD_all = sum(c .* sum(Pm .* D2, 2)) / N;
    r.coopPred = sum(c .* sum(Pm .* A1, 2)) / N;
    r.wbar = sum(c .* sum(Pm .* payoff, 2)) / N;
    r.wMem = weighted_member_mean(Pm .* payoff, c, mem);
    r.s = r.fD * (P - S) / r.wbar;
end


% =========================================================================
function m = weighted_member_mean(X, c, mem)
    if ~any(mem), m = NaN; return; end
    m = sum(c(mem) .* sum(X(mem, :), 2)) / sum(c(mem));
end


% =========================================================================
function m = defect_from_members(A2, Pm, c, mem)
% Share of members' defector encounters in which the defecting partner is itself a member.
    if ~any(mem), m = NaN; return; end
    X = Pm .* (1 - A2);
    tot = sum(c(mem) .* sum(X(mem, :), 2));
    m = sum(c(mem) .* sum(X(mem, mem), 2)) / max(tot, eps);
    if tot == 0, m = NaN; end
end


% =========================================================================
function out = assemble(res, jobs, meta, ARCH)
    out = meta;
    for a = 1:size(ARCH, 1)
        J = find(cellfun(@(j) j.arch == a, jobs));
        rep = cellfun(@(j) j.rep, jobs(J)); gen = cellfun(@(j) j.gen, jobs(J));
        hc = cellfun(@(j) j.histCoop, jobs(J));
        g500 = J(gen == 500);
        out(a).rep500 = rep(gen == 500);
        out(a).at500 = [res{g500}];
        out(a).histCoop500 = hc(gen == 500);
        out(a).window = struct('rep', num2cell(rep(:)), 'gen', num2cell(gen(:)), ...
            'fD', num2cell(cellfun(@(x) x.fD, res(J))), ...
            'fD_eng', num2cell(cellfun(@(x) x.eng.fD, res(J))));
    end
end


% =========================================================================
function rates = retaliation_removal_rates()
% The rate u = lambda x P(a mutation ends only retaliation) per member per generation (SI S4.1),
% from fixed values measured by other files.
%   lambda = 0.005 per individual per generation (Poisson; Table S1).
%   Pooled: the per-mutation census (tableS6_per_mutation_shares), share ending membership 0.360,
%     ending cooperation 0.125; no mutation ends cooperation without ending membership, so
%     the share ending only retaliation is 0.360 - 0.125.
%   Per architecture (two measurements combined): the snapshot share of mutations ending
%     membership, 0.591 / 0.415 / 0.421 / 0.381 (SI S3.2; share_ending_membership), times the
%     share of those that leave cooperation intact, 0.599 / 0.668 / 0.654 / 0.683
%     (membership_exits_split, 20 to 25 members per architecture).
    lambda = 0.005;
    rates.pooled = lambda * (0.360 - 0.125);
    rates.perArch = lambda * [0.591 0.415 0.421 0.381] .* [0.599 0.668 0.654 0.683];
end



% =========================================================================
function s = vlabel(v)
% The printed name of a variant or field; 'eng' and 'fD_eng' are the names stored in the saved .mat.
    s = strrep(strrep(v, 'fD_eng', 'fD, encounter rule'), 'eng', 'encounter rule');
end


% =========================================================================
function X = variant(X, v)
% 'det600': the top-level fields; otherwise the named sub-struct ('eng', 'rand600', 'det300').
    if ~strcmp(v, 'det600') && ~isempty(X), X = [X.(v)]; end
end


% =========================================================================
function print_report(out, cfg)
    rates = retaliation_removal_rates();
    fprintf('\nDEFECTOR-ENCOUNTER FREQUENCY OF MEMBERS  (run %s UTC, seed %d, %d draws/pair on random-start variants)\n', ...
        char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm')), cfg.rngSeed, cfg.engineDraws);
    fprintf('Reading: median member f_D over successful replicates < 0.1%% -> selection on retaliation relaxed; >= 1%% -> exercised; between -> ambiguous.\n');
    fprintf('Variants: det600 = simulate_pair 600 steps, x0 = 0.005 (the membership test''s dynamics);\n');
    fprintf('          encounter rule = the runs'' own rule, 300 steps, x0 ~ U(0, 0.01) per agent per encounter (the paper''s f_D);\n');
    fprintf('          rand600 / det300 = the two differences one at a time.\n');
    print_inputs(out);
    V = {'det600', 'eng', 'rand600', 'det300'};
    for i = 1:numel(V)
        fprintf('\n[1.%d] %s, gen 500, SUCCESSFUL replicates\n', i, vlabel(V{i}));
        print_fd_table(out, true, V{i});
    end
    for i = 1:2
        fprintf('\n[2.%d] CONTROL: %s, gen 500, FAILED replicates (n = replicates with >= 1 member)\n', i, vlabel(V{i}));
        print_fd_table(out, false, V{i});
    end
    print_crosscheck(out, 'det600');
    print_crosscheck(out, 'eng');
    print_selection(out, rates, 'det600');
    print_selection(out, rates, 'eng');
    print_window(out, cfg);
    print_diagnostics(out);
end


% =========================================================================
function print_inputs(out)
    fprintf('\nInputs (file date | nOpp | payoffMode | development steps | success, recomputed at >= 0.95):\n');
    for a = 1:numel(out)
        fprintf('  %-20s %s | %d | %s | %d | %d/%d\n', out(a).name, out(a).fileDate, out(a).nOpp, ...
            out(a).payoffMode, out(a).engineSteps, sum(out(a).success), numel(out(a).success));
    end
end


% =========================================================================
function print_fd_table(out, wantSuccess, v)
    fprintf('  %-20s %4s %5s | %9s %21s %21s | %9s %9s | %6s\n', 'architecture', 'n', 'noMem', ...
        'median fD', 'IQR', 'range', 'retal.', 'exploit.', 'mem sh');
    for a = 1:numel(out)
        X = variant(out(a).at500(out(a).success(out(a).rep500) == wantSuccess), v);
        if isempty(X), fprintf('  %-20s    0 (no replicates)\n', out(a).name); continue; end
        fD = [X.fD]; ok = ~isnan(fD);
        if ~any(ok), fprintf('  %-20s    0 %5d (no members)\n', out(a).name, sum(~ok)); continue; end
        ret = [X.fRetaliated]; exl = [X.fExploited];
        q = prctile(fD(ok), [25 75]);
        fprintf('  %-20s %4d %5d | %8.3f%% [%7.3f%%, %7.3f%%] [%7.3f%%, %7.3f%%] | %8.3f%% %8.3f%% | %5.1f%%\n', ...
            out(a).name, sum(ok), sum(~ok), 100*median(fD(ok)), 100*q, 100*min(fD(ok)), 100*max(fD(ok)), ...
            100*median(ret(ok)), 100*median(exl(ok)), 100*mean([X.occ]));
    end
    pooled = arrayfun(@(o) variant(o.at500(o.success(o.rep500) == wantSuccess), v), out, 'uni', 0);
    pooled = [pooled{~cellfun(@isempty, pooled)}];
    fAll = [pooled.fD]; fAll = fAll(~isnan(fAll));
    fprintf('  %-20s %4d       | %8.3f%% [%7.3f%%, %7.3f%%]   (four architectures pooled)\n', 'pooled', ...
        numel(fAll), 100*median(fAll), 100*prctile(fAll, [25 75]));
    fprintf('  (retal. / exploit. = member defects back / keeps cooperating, as shares of ALL member encounters, medians)\n');
    fprintf('  (granularity: one defecting individual among the N - 1 = 999 possible partners gives f_D = 0.100%%)\n');
end


% =========================================================================
function print_selection(out, rates, v)
    fprintf('\n[4:%s] IMPLIED SELECTION ON THE RETALIATION SIGN (successful, gen 500, medians over reps)\n', vlabel(v));
    fprintf('  s = f_D (P - S) / wbar;  u = lambda * (share of mutations that end only retaliation) per member per generation\n');
    fprintf('  %-20s %9s %8s %11s %11s %8s\n', 'architecture', 'f_D', 'wbar', 's', 'u (arch)', 's/u');
    for a = 1:numel(out)
        X = variant(out(a).at500(out(a).success(out(a).rep500)), v);
        if isempty(X), fprintf('  %-20s (no successful replicates)\n', out(a).name); continue; end
        s = median([X.s], 'omitnan');
        fprintf('  %-20s %8.3f%% %8.3f %11.2e %11.2e %8.2f\n', out(a).name, 100*median([X.fD], 'omitnan'), ...
            median([X.wbar]), s, rates.perArch(a), s / rates.perArch(a));
    end
    fprintf('  pooled u = 0.005 x (0.360 - 0.125) = %.2e (pooled per-mutation census)\n', rates.pooled);
end


% =========================================================================
function print_crosscheck(out, v)
    fprintf('\n[3:%s] CROSS-CHECK: predicted focal defection vs recorded 1 - historyCooperation(500) (one realized pairing, N = 1000)\n', vlabel(v));
    fprintf('  %-20s %4s %11s %11s %11s %11s\n', 'architecture', 'rep', '1-recorded', '1-predicted', 'pop f_D', 'member f_D');
    for a = 1:numel(out)
        X = variant(out(a).at500, v);
        idx = find(out(a).success(out(a).rep500));
        for i = idx(1:min(3, numel(idx)))'
            fprintf('  %-20s %4d %10.3f%% %10.3f%% %10.3f%% %10.3f%%\n', out(a).name, out(a).rep500(i), ...
                100*(1 - out(a).histCoop500(i)), 100*(1 - X(i).coopPred), 100*X(i).fD_all, 100*X(i).fD);
        end
        pred = 1 - [X.coopPred]; obs = 1 - out(a).histCoop500;
        z = (obs(:) - pred(:)) ./ sqrt(max(pred(:) .* (1 - pred(:)), 1 / 1000) / 1000);
        sel = out(a).success(out(a).rep500);
        fprintf('  %-20s all reps: mean |obs - pred| = %.2f pts; median |z| = %.2f (successful %.2f, failed %.2f); share |z| > 3 = %.2f\n', ...
            out(a).name, 100*mean(abs(obs(:) - pred(:))), median(abs(z)), median(abs(z(sel))), ...
            median(abs(z(~sel))), mean(abs(z) > 3));
    end
end


% =========================================================================
function print_window(out, cfg)
    fprintf('\n[5] STATIONARITY, successful replicates, snapshots in [%d, 500]\n', cfg.windowStart);
    for a = 1:numel(out)
        w = out(a).window; w = w(out(a).success([w.rep]));
        if isempty(w), fprintf('  %-20s (no successful replicates)\n', out(a).name); continue; end
        gens = unique([w.gen]); reps = unique([w.rep]);
        for fld = {'fD', 'fD_eng'}
            med = arrayfun(@(g) median([w([w.gen] == g).(fld{1})], 'omitnan'), gens);
            repMean = arrayfun(@(r) mean([w([w.rep] == r).(fld{1})], 'omitnan'), reps);
            repMed = arrayfun(@(r) median([w([w.rep] == r).(fld{1})], 'omitnan'), reps);
            fprintf('  %-20s %-6s median over reps of window-mean = %.3f%%, of window-median = %.3f%% | per-snapshot medians: %s\n', ...
                out(a).name, vlabel(fld{1}), 100*median(repMean, 'omitnan'), 100*median(repMed, 'omitnan'), ...
                sprintf('%d:%.3f%% ', [gens; 100*med]));
        end
        hc = out(a).histCoop(out(a).success, cfg.windowStart:500);
        fprintf('  %-20s population proxy 1 - historyCooperation(g), g = %d..500: median over reps of mean = %.3f%%\n', ...
            '', cfg.windowStart, 100*median(mean(1 - hc, 2)));
    end
end


% =========================================================================
function print_diagnostics(out)
    fprintf('\n[6] DIAGNOSTICS (gen 500, all replicates; det600 unless marked)\n');
    for a = 1:numel(out)
        X = out(a).at500; sel = out(a).success(out(a).rep500);
        E = [X.eng];
        fprintf(['  %-20s role-symmetry mismatches %d | unique wirings median %d | member tie share max %.2e | ' ...
            'residue-decided defections max %.2e\n'], out(a).name, sum([X.symMismatch]), ...
            median([X.nUnique]), max([X.memTieShare]), max([X.memResidueShare]));
        fprintf('  %-20s share of members'' defector encounters where the defector is itself a member (successful, median): det600 %.2f, encounter rule %.2f\n', ...
            '', median([X(sel).memDefectorIsMember], 'omitnan'), median([E(sel).memDefectorIsMember], 'omitnan'));
    end
end
